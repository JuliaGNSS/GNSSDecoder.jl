abstract type AbstractGNSSConstants end
abstract type AbstractGNSSData end
abstract type AbstractGNSSCache end

"""
    AbstractGPSData <: AbstractGNSSData

Abstract supertype for the decoded navigation data of a signal transmitted by
the GPS constellation, e.g. `GPSL1CAData`, `GPSCNAVData`.

Its purpose is to carry the constellation-level facts every GPS signal's data
shares, so they can be stated once (on the supertype, via subtype dispatch)
instead of once per signal. Constellation membership is encoded at the struct
definition site — the `<: AbstractGPSData` line written anyway — so a new GPS
signal inherits the shared behaviour with nothing to remember. Genuinely
per-signal facts (the subframe/message-type completeness checks, the health-bit
selection in `is_sat_healthy`) stay defined on the concrete data types.
"""
abstract type AbstractGPSData <: AbstractGNSSData end

"""
    AbstractGalileoData <: AbstractGNSSData

Abstract supertype for the decoded navigation data of a signal transmitted by
the Galileo constellation, e.g. `GalileoE1BData`, `GalileoE5aData`.

The Galileo counterpart to [`AbstractGPSData`](@ref). It carries the facts every
Galileo signal's data shares: `is_ephemeris_decoded` and
`is_clock_correction_decoded` check the same orbital and clock fields for I/NAV
(E1B) and F/NAV (E5a), so they are defined once on this supertype (see
`src/galileo/galileo.jl`) instead of once per signal. The health-status and
positioning-readiness checks genuinely differ per signal and stay on the
concrete data types.
"""
abstract type AbstractGalileoData <: AbstractGNSSData end

# Physical constants common to every GNSS handled here. Each per-signal
# `*Constants` struct exposes these as fields (so the orbit/clock math reads
# `state.constants.PI` etc.); the defaults are sourced from here to keep a single
# source of truth. Constellation-specific values that genuinely differ by
# reference frame — the Earth gravitational parameter μ and the relativistic
# correction F — are *not* here: GPS uses WGS-84 values and Galileo the GTRF
# values, defined alongside each constellation.
#
# `GNSS_PI` is deliberately the truncated value the ICDs fix for the
# semicircle→radian scaling (IS-GPS-200 Table 20-IV; Galileo OS SIS ICD Table 68),
# *not* `Base.π` — broadcast angular quantities must be scaled with exactly this
# value to reproduce the transmitted numbers bit-for-bit.
const GNSS_PI = 3.1415926535898
const SPEED_OF_LIGHT = 2.99792458e8        # m/s
const EARTH_ROTATION_RATE = 7.2921151467e-5  # rad/s (WGS-84 and GTRF agree)

"""
$(TYPEDEF)

Generic decoder state for GNSS signal decoding. This parametric struct holds all state
required for decoding navigation messages from GNSS satellites.

The struct itself is immutable; per-field reconstruction works via the keyword
constructor, and the per-signal constants and decoded data carry value
semantics. The one piece of intentionally-mutable state is the soft-symbol
buffer inside the `cache`: a `CircularDeque{Float32}` of capacity
`syncro_sequence_length + preamble_length` that accumulates incoming symbols
across successive [`decode`](@ref) calls. It is a mutable container shared by
reference between an input state and the state `decode` returns — fully
immutable threading would copy the whole buffer on every symbol, which is the
wrong trade for a streaming decoder. Treat the value returned by `decode` as
*the* live state and do not keep mutating an earlier snapshot in parallel. The
transient packed-bit buffer used for preamble matching is **not** stored here;
it is computed as a local value at sync time and threaded through the sync
path (see `pack_buffer` / `try_sync`).

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

calc_preamble_mask(constants::AbstractGNSSConstants) =
    UInt(1) << UInt(constants.preamble_length) - UInt(1)

"""
    pack_buffer(state) -> Unsigned

Hard-slice the leading `syncro_sequence_length + preamble_length` soft
symbols of the per-signal soft buffer into a packed-bit value (the v1
`raw_buffer`: oldest bit at MSB, newest bit at LSB). The concrete unsigned
type is signal-specific and supplied by `packed_buffer_type`. The
result is a plain value — it is threaded through the sync path rather than
stashed in mutable cache state.
"""
function pack_buffer(state::GNSSDecoderState)
    n = state.constants.syncro_sequence_length + state.constants.preamble_length
    pack_soft_buffer(packed_buffer_type(state), soft_buffer(state), n)
end

"""
    try_sync(state) -> Union{Nothing,Unsigned}

Default per-signal sync hook: hard-slice the deque into a packed-bit buffer
(via `pack_buffer`) and run the `find_preamble` bit-pattern
check (preamble visible at both ends of the candidate syncro sequence, in
either polarity). Returns the packed buffer on a match, or `nothing` if there
is no sync. Returning the buffer lets the caller reuse it without recomputing
and keeps it out of mutable cache state.

Per-signal overrides (e.g. GPS L1C-D's TOI BCH match in a later slice)
override this method.
"""
function try_sync(state::GNSSDecoderState)
    buffer = pack_buffer(state)
    find_preamble(buffer, state.constants) ? buffer : nothing
end

"""
    find_preamble(buffer, constants) -> Bool

Bit-pattern preamble check on a packed-bit buffer. Mirrors the v1
implementation: the preamble must be visible at *both* the oldest 8 bits
(start of this subframe) and the newest 8 bits (start of next subframe),
either both upright OR both inverted (180-degree polarity ambiguity).
"""
function find_preamble(buffer, constants::AbstractGNSSConstants)
    mask = calc_preamble_mask(constants)
    buffer & mask == constants.preamble &&
        (buffer >> constants.syncro_sequence_length) & mask == constants.preamble ||
        buffer & mask == ~constants.preamble & mask &&
            (buffer >> constants.syncro_sequence_length) & mask ==
            ~constants.preamble & mask
end

"""
    complement_buffer_if_necessary(state, buffer) -> (state, resolved_buffer)

If the newest preamble in `buffer` is the *inverted* preamble, return the
state flagged `is_shifted_by_180_degrees = true` together with the
complemented buffer; otherwise return the state flagged `false` with `buffer`
unchanged. The polarity-resolved buffer is returned as a value (the v1
`buffer`) for the per-signal `decode_syncro_sequence` to consume. Mirrors v1
behaviour.
"""
function complement_buffer_if_necessary(state::GNSSDecoderState, buffer)
    mask = calc_preamble_mask(state.constants)
    if buffer & mask == ~state.constants.preamble & mask
        return GNSSDecoderState(state; is_shifted_by_180_degrees = true), ~buffer
    else
        return GNSSDecoderState(state; is_shifted_by_180_degrees = false), buffer
    end
end

"""
    drain_after_sync!(state)

Drop the consumed `syncro_sequence_length` oldest soft symbols from the
deque, keeping the trailing `preamble_length` symbols as the leading
preamble of the next subframe. Equivalent to v1's
`GNSSDecoderState(state; num_bits_buffered = preamble_length)`.

Drops at most `length(deque)` symbols: a `decode_syncro_sequence` hook may
reset the decoder mid-frame (e.g. GPS L1C-D on a TOI discontinuity), which
empties the buffer. Without the clamp the unconditional drain in `decode`
would `popfirst!` an empty `CircularDeque` and throw.
"""
function drain_after_sync!(state::GNSSDecoderState)
    deque = soft_buffer(state)
    n_drop = min(state.constants.syncro_sequence_length, length(deque))
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

`soft_symbols` is an `AbstractVector{<:Real}`; `Float32` is canonical. The sign
carries the bit decision and the magnitude carries confidence (standard LLR
convention):

- **positive ⇒ bit 0**, **negative ⇒ bit 1** — but treat this as a *convention*,
  not a hard input requirement. The absolute polarity of a Costas-tracked signal
  is inherently 180°-ambiguous, so the decoder does not depend on it: it matches
  the preamble in either polarity and flips internally (recording the result in
  `is_shifted_by_180_degrees`). Feeding the opposite sign decodes the same data;
  only the reported polarity flag differs. (Note: `Tracking.jl`'s
  `get_soft_bits` happens to use the opposite sign — positive ⇒ bit 1 — which is
  harmless for exactly this reason.)
- magnitude ⇒ confidence. **No normalization is required; values need not lie in
  `[-1, 1]`.** GPS L1 C/A (hard-slice + parity) and Galileo E1B (Viterbi, whose
  ML path is invariant to a global scale) use the sign and are indifferent to
  the magnitude scale. GPS L1C-D's LDPC decode is flooding sum-product, which
  *is* scale-sensitive, so there the magnitudes should be confidence-weighted on
  a roughly LLR-like scale (`≈ 2·r/σ²`) for best performance at marginal SNR —
  but still need not be normalized to a fixed range.

Glue from `Tracking.jl`: feed `get_soft_bits` (polarity-corrected,
amplitude-weighted soft bits) for every signal. See `CONTEXT.md` for the full
glossary.

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

        if is_enough_buffered_bits_to_decode(state)
            buffer = try_sync(state)
            if !isnothing(buffer)
                state, resolved_buffer = complement_buffer_if_necessary(state, buffer)
                state = decode_syncro_sequence(state, resolved_buffer)
                if !decode_once || !is_decoding_completed_for_positioning(state.data)
                    state = validate_data(state)
                end
                state = drain_after_sync!(state)
            end
        end
    end
    return state
end

# ---- Shared GPS decoder primitives ------------------------------------------
#
# Signal-agnostic primitives used by more than one GPS signal decoder. They
# live here (a shared file included before every signal) rather than in a
# per-signal file so that no signal decoder has to be included after another
# just to borrow them.

# Packed-word integer type shared across the GPS decoders: a GPS L1 C/A
# subframe, a GPS L1C-D subframe, and a GPS CNAV message (GPS L5I / L2C) each
# pack into a single `UInt320` — 300 data bits plus up to 8 trailing sync bits.
# `BitIntegers.@define_integers` also defines the signed companion `Int320`.
BitIntegers.@define_integers 320

"""
Insert/overwrite `value` keyed by `key` in a (possibly `nothing`) `Dictionary`, returning the updated copy.
"""
function _merge_keyed(dict::Union{Nothing,Dictionary{Int,V}}, key::Int, value::V) where {V}
    out = isnothing(dict) ? Dictionary{Int,V}() : copy(dict)
    set!(out, key, value)
    return out
end
