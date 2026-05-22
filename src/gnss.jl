abstract type AbstractGNSSConstants end
abstract type AbstractGNSSData end
abstract type AbstractGNSSCache end

"""
$(TYPEDEF)

Generic decoder state for GNSS signal decoding. This parametric struct holds all state
required for decoding navigation messages from GNSS satellites.

Soft symbols are buffered inside the mutable per-signal `cache` (a
`CircularDeque{Float32}` with the per-signal capacity of
`syncro_sequence_length + preamble_length`). The outer struct itself is
immutable; per-field reconstruction continues to work via the keyword
constructor.

# Type Parameters
- `D<:AbstractGNSSData`: The data type holding decoded navigation message fields
- `C<:AbstractGNSSConstants`: Constants specific to the GNSS system (e.g., preamble, timing)
- `CA<:AbstractGNSSCache`: Cache for intermediate decoding state (carries the soft-symbol buffer)

# Fields
$(TYPEDFIELDS)

# See Also
- [`GPSL1CADecoderState`](@ref): Constructor for GPS L1 C/A decoder state
- [`GalileoE1BDecoderState`](@ref): Constructor for Galileo E1B decoder state
- [`decode`](@ref): Main function to decode soft symbols using this state
- [`reset_decoder_state`](@ref): Reset decoder state after signal loss
"""
Base.@kwdef struct GNSSDecoderState{
    D<:AbstractGNSSData,
    C<:AbstractGNSSConstants,
    CA<:AbstractGNSSCache,
}
    "Pseudo-Random Noise code identifier for the satellite"
    prn::Int
    "Partially decoded navigation data (not yet validated)"
    raw_data::D
    "Validated navigation data ready for use"
    data::D
    "System-specific constants (preamble, timing parameters)"
    constants::C
    "Cache for intermediate decoding state (holds the soft-symbol `CircularDeque{Float32}`)"
    cache::CA
    "Number of symbols received after the last valid synchronization sequence, or `nothing` if not yet synchronized"
    num_bits_after_valid_syncro_sequence::Union{Nothing,Int} = 0
    "Whether the signal phase is inverted by 180 degrees"
    is_shifted_by_180_degrees::Bool = false
end

function GNSSDecoderState(
    state::GNSSDecoderState;
    raw_data = state.raw_data,
    data = state.data,
    cache = state.cache,
    num_bits_after_valid_syncro_sequence = state.num_bits_after_valid_syncro_sequence,
    is_shifted_by_180_degrees = state.is_shifted_by_180_degrees,
)
    GNSSDecoderState(
        state.prn,
        raw_data,
        data,
        state.constants,
        cache,
        num_bits_after_valid_syncro_sequence,
        is_shifted_by_180_degrees,
    )
end

# The default `==` for structs containing fields with mutable types (like the
# `CircularDeque{Float32}` soft-symbol buffer or the `Vector{GalileoAlmanac}`
# inside `GalileoE1BData`) falls back to `===`. Compare field-by-field so that
# two states with equal-but-not-identical contents are considered equal.
function Base.:(==)(a::GNSSDecoderState, b::GNSSDecoderState)
    typeof(a) === typeof(b) || return false
    for f in fieldnames(typeof(a))
        getfield(a, f) == getfield(b, f) || return false
    end
    return true
end

"""
    deques_equal(a::CircularDeque, b::CircularDeque)

Internal helper: structural equality on `CircularDeque`s of the same element
type. DataStructures.jl does not define `==` on `CircularDeque`, and we
deliberately avoid defining it ourselves (which would be type piracy and is
flagged by Aqua). Per-signal caches' `==` calls this directly.
"""
function deques_equal(a::CircularDeque{T}, b::CircularDeque{T}) where {T}
    length(a) == length(b) || return false
    capacity(a) == capacity(b) || return false
    for i in 1:length(a)
        a[i] == b[i] || return false
    end
    return true
end

"Soft-symbol buffer accessor — the per-signal cache stores it as `soft_buffer`."
soft_buffer(state::GNSSDecoderState) = state.cache.soft_buffer

"Number of soft symbols currently buffered."
num_bits_buffered(state::GNSSDecoderState) = length(soft_buffer(state))

"""
    push_soft_symbol!(state, sym)

Push one `Float32` soft symbol onto the per-signal circular deque. The deque
is sized at `syncro_sequence_length + preamble_length`; once full, the oldest
sample is overwritten via `popfirst!`.
"""
function push_soft_symbol!(state::GNSSDecoderState, sym::Real)
    deque = soft_buffer(state)
    if length(deque) >= capacity(deque)
        popfirst!(deque)
    end
    push!(deque, Float32(sym))
    return state
end

function is_enough_buffered_bits_to_decode(state::GNSSDecoderState)
    num_bits_buffered(state) >=
    state.constants.syncro_sequence_length + state.constants.preamble_length
end

"""
    hard_slice(soft_symbol) -> Bool

Convention: positive soft symbol ⇒ bit 0, negative ⇒ bit 1. (Matches AFF3CT's
LLR convention.) The returned `Bool` is `true` for bit 1.
"""
@inline hard_slice(sym::Real) = sym < zero(sym)

"""
    pack_soft_buffer(T, soft_buffer, total_bits)

Hard-slice the leading `total_bits` of `soft_buffer` (oldest first) into a
`T<:Unsigned` packed-bit buffer, MSB = oldest bit. Mirrors how the legacy
`push_bit` shifted bits into `raw_buffer`.
"""
function pack_soft_buffer(::Type{T}, deque::CircularDeque{Float32}, total_bits::Int) where {T<:Unsigned}
    word = T(0)
    @inbounds for i in 1:total_bits
        bit = hard_slice(deque[i]) ? T(1) : T(0)
        word = (word << 1) | bit
    end
    return word
end

calc_preamble_mask(state::GNSSDecoderState) =
    UInt(1) << UInt(state.constants.preamble_length) - UInt(1)

"""
    pack_buffer_into_cache!(state)

Hard-slice the leading `syncro_sequence_length + preamble_length` soft
symbols of `state.cache.soft_buffer` into the per-signal packed-bit buffer
held inside `state.cache.packed_buffer[]`. The packed buffer is the v1
`raw_buffer`: oldest bit at MSB, newest bit at LSB.
"""
function pack_buffer_into_cache!(state::GNSSDecoderState)
    cache = state.cache
    n = state.constants.syncro_sequence_length + state.constants.preamble_length
    T = eltype(cache.packed_buffer)
    cache.packed_buffer[] = pack_soft_buffer(T, soft_buffer(state), n)
    state
end

"""
    try_sync(state)

Default per-signal sync hook: hard-slice the deque tail into the per-signal
packed-bit buffer (via `pack_buffer_into_cache!`) and run the existing
`find_preamble` bit-pattern check (preamble visible at both ends of the
candidate syncro sequence, in either polarity). Returns the (possibly
modified) state and a `Bool` match flag.

Per-signal overrides (e.g. GPS L1C-D's TOI BCH match in a later slice)
override this method.
"""
function try_sync(state::GNSSDecoderState)
    pack_buffer_into_cache!(state)
    find_preamble(state)
end

"""
    find_preamble(state)

Bit-pattern preamble check on the per-signal packed-bit buffer. Mirrors the
v1 implementation: the preamble must be visible at *both* the oldest 8 bits
(start of this subframe) and the newest 8 bits (start of next subframe),
either both upright OR both inverted (180-degree polarity ambiguity).
"""
function find_preamble(state::GNSSDecoderState)
    raw = state.cache.packed_buffer[]
    state.cache.packed_buffer[] & calc_preamble_mask(state) == state.constants.preamble &&
        (raw >> state.constants.syncro_sequence_length) &
        calc_preamble_mask(state) == state.constants.preamble ||
        raw & calc_preamble_mask(state) ==
        ~state.constants.preamble & calc_preamble_mask(state) &&
            (raw >> state.constants.syncro_sequence_length) &
            calc_preamble_mask(state) ==
            ~state.constants.preamble & calc_preamble_mask(state)
end

"""
    complement_buffer_if_necessary(state)

If the newest preamble in the packed buffer is the *inverted* preamble,
record `is_shifted_by_180_degrees = true` and store the complemented buffer.
Otherwise store the packed buffer as-is. Mirrors v1 behaviour.
"""
function complement_buffer_if_necessary(state::GNSSDecoderState)
    raw = state.cache.packed_buffer[]
    if raw & calc_preamble_mask(state) ==
       ~state.constants.preamble & calc_preamble_mask(state)
        state.cache.complemented_buffer[] = ~raw
        return GNSSDecoderState(state; is_shifted_by_180_degrees = true)
    else
        state.cache.complemented_buffer[] = raw
        return GNSSDecoderState(state; is_shifted_by_180_degrees = false)
    end
end

"""
    drain_after_sync!(state)

Drop the consumed `syncro_sequence_length` oldest soft symbols from the
deque, keeping the trailing `preamble_length` symbols as the leading
preamble of the next subframe. Equivalent to v1's
`GNSSDecoderState(state; num_bits_buffered = preamble_length)`.
"""
function drain_after_sync!(state::GNSSDecoderState)
    deque = soft_buffer(state)
    n_drop = state.constants.syncro_sequence_length
    for _ in 1:n_drop
        popfirst!(deque)
    end
    state
end

"""
$(TYPEDSIGNATURES)

Decode GNSS navigation message soft symbols and update the decoder state.

Processes incoming soft symbols from a GNSS signal, detecting preambles and
decoding synchronization sequences to extract navigation data. The function
handles both normal and 180-degree phase-shifted signals automatically.

# Soft-symbol convention

`soft_symbols` is an `AbstractVector{<:Real}`; `Float32` is canonical. The
sign carries the bit decision and the magnitude carries confidence (matches
AFF3CT's LLR convention):

- **positive ⇒ transmitted bit 0**
- **negative ⇒ transmitted bit 1**
- magnitude ⇒ confidence (proportional to SNR × coherent integration)

Glue from `Tracking.jl` typically supplies `Float32(real(prompt))` after
bit-sync resolves polarity. See `CONTEXT.md` for the full glossary.

# Arguments
- `state::GNSSDecoderState`: Current decoder state
- `soft_symbols::AbstractVector{<:Real}`: Soft symbols to consume, oldest first
- `num_symbols::Int`: Number of leading entries of `soft_symbols` to process

# Keywords
- `decode_once::Bool=false`: If `true`, stops once all required positioning
  data has been validated (subframes 1-3 for GPS L1 C/A; word types 1-5 for
  Galileo E1B)

# Returns
- `GNSSDecoderState`: Updated decoder state with newly decoded data

# Example
```julia
state = GPSL1CADecoderState(1)            # PRN 1
state = decode(state, Float32[+1, -1, ...], 8)
```

# See Also
- [`GNSSDecoderState`](@ref): The state structure being updated
- [`is_sat_healthy`](@ref): Check satellite health after decoding
"""
function decode(
    state::GNSSDecoderState,
    soft_symbols::AbstractVector{<:Real},
    num_symbols::Int;
    decode_once::Bool = false,
)
    num_symbols <= length(soft_symbols) ||
        throw(ArgumentError("num_symbols exceeds length(soft_symbols)"))
    for i in 1:num_symbols
        sym = soft_symbols[i]
        state = push_soft_symbol!(state, sym)
        if !isnothing(state.num_bits_after_valid_syncro_sequence)
            state = GNSSDecoderState(
                state;
                num_bits_after_valid_syncro_sequence = state.num_bits_after_valid_syncro_sequence +
                                                       1,
            )
        end

        if is_enough_buffered_bits_to_decode(state) && try_sync(state)
            state = decode_syncro_sequence(complement_buffer_if_necessary(state))
            if !decode_once || !is_decoding_completed_for_positioning(state.data)
                state = validate_data(state)
            end
            state = drain_after_sync!(state)
        end
    end
    return state
end
