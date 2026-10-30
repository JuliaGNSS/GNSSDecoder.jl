"""
Soft-symbol buffer accessor — the per-signal cache stores it as `soft_buffer`.
"""
soft_buffer(state::GNSSDecoderState) = state.cache.soft_buffer

"""
Number of soft symbols currently buffered.
"""
num_symbols_buffered(state::GNSSDecoderState) = length(soft_buffer(state))

"""
    push_soft_symbol!(state, sym)

Push one `Float32` soft symbol onto the per-signal circular deque. The deque
is sized at `sync_unit_length + preamble_length`; once full, the oldest
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

function is_enough_buffered_symbols_to_decode(state::GNSSDecoderState)
    num_symbols_buffered(state) >=
    state.constants.sync_unit_length + state.constants.preamble_length
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
function pack_soft_buffer(
    ::Type{T},
    deque::CircularDeque{Float32},
    total_bits::Int,
) where {T<:Unsigned}
    word = T(0)
    @inbounds for i = 1:total_bits
        bit = hard_slice(deque[i]) ? T(1) : T(0)
        word = (word << 1) | bit
    end
    return word
end

calc_preamble_mask(constants::AbstractGNSSConstants) =
    UInt(1) << UInt(constants.preamble_length) - UInt(1)

"""
    pack_bits_msb_first(deque, offset, len) -> UInt64

Hard-slice `len` soft symbols of `deque` starting at 1-based `offset` into a
`UInt64`, oldest symbol at the MSB — the same bit order [`pack_buffer`](@ref)
produces, but for a slice rather than the whole window.

Sync only ever inspects a handful of known positions (the preamble at either
end, and any unencoded header field), so a decoder whose window is long and
whose symbol rate is high can read exactly those bits instead of repacking
the entire window on every symbol.

`len` must be at most 64, and `offset + len - 1` must be within the deque: the
read is `@inbounds`, so the caller owns both bounds. Every caller is a
`try_sync` hook, which the decode loop runs only behind
`is_enough_buffered_symbols_to_decode`.
"""
function pack_bits_msb_first(deque::CircularDeque{Float32}, offset::Int, len::Int)
    word = UInt64(0)
    @inbounds for i = 0:(len-1)
        word = (word << 1) | (hard_slice(deque[offset+i]) ? UInt64(1) : UInt64(0))
    end
    return word
end

"""
    pack_bits_lsb_first(deque, offset, len) -> UInt64

Hard-slice `len` soft symbols of `deque` starting at 1-based `offset` into a
`UInt64` with **bit 0 = the first symbol** — the opposite end from
[`pack_bits_msb_first`](@ref), and the order the codeword tables use
(`BCH_TOI_CODEWORDS`, `b1c_prn_codeword`, `b1c_soh_codeword`; cf.
[`pack_hard_codeword`](@ref), which is this over an iterable of hard decisions).

The two BCH-synchronised signals — GPS L1C-D on its TOI codeword and BeiDou B1C
on its PRN + SOH pair — compare against those tables once per symbol, at both
ends of the window, so this reads the deque in place rather than materialising a
slice.

`len` must be at most 64, and `offset + len - 1` must be within the deque: the
read is `@inbounds`, so the caller owns both bounds — the length check that
`pack_hard_codeword` performs on a materialised slice is not available here.
Every caller is a `try_sync` hook, which the decode loop runs only behind
`is_enough_buffered_symbols_to_decode`.
"""
function pack_bits_lsb_first(deque::CircularDeque{Float32}, offset::Int, len::Int)
    word = UInt64(0)
    @inbounds for i = 0:(len-1)
        if hard_slice(deque[offset+i])
            word |= UInt64(1) << i
        end
    end
    return word
end

"""
    pack_buffer(state) -> Unsigned

Hard-slice the leading `sync_unit_length + preamble_length` soft
symbols of the per-signal soft buffer into a packed-bit value (oldest bit at
MSB, newest bit at LSB). The concrete unsigned
type is signal-specific and supplied by `packed_buffer_type`. The
result is a plain value — it is threaded through the sync path rather than
stashed in mutable cache state.
"""
function pack_buffer(state::GNSSDecoderState)
    n = state.constants.sync_unit_length + state.constants.preamble_length
    pack_soft_buffer(packed_buffer_type(state), soft_buffer(state), n)
end

"""
    try_sync(state) -> Union{Nothing,Unsigned}

Default per-signal sync hook: hard-slice the deque into a packed-bit buffer
(via `pack_buffer`) on a preamble match (preamble visible at both ends of the
candidate sync unit, in either polarity; see `find_preamble_in_deque`). Returns the packed buffer on a match, or `nothing` if there
is no sync. Returning the buffer lets the caller reuse it without recomputing
and keeps it out of mutable cache state.

Per-signal overrides (e.g. GPS L1C-D's TOI BCH match in a later slice)
override this method.
"""
function try_sync(state::GNSSDecoderState)
    # Match the preamble on the few symbols at either end first, and pack the
    # whole window only on a match: packing it (hundreds of symbols) for every
    # incoming symbol is otherwise most of the per-symbol cost.
    constants = state.constants
    isnothing(
        find_preamble_in_deque(
            soft_buffer(state),
            constants.preamble,
            constants.preamble_length,
            constants.sync_unit_length,
        ),
    ) && return nothing
    return pack_buffer(state)
end

"""
    find_preamble_in_deque(deque, preamble, preamble_length, sync_unit_length)
        -> Union{Nothing,Bool}

Match the preamble at both ends of the candidate window — the start of this sync
unit and the start of the next — by slicing the soft-symbol deque directly,
rather than hard-slicing the whole window into a packed integer first.

Returns `nothing` when there is no sync, and otherwise the resolved polarity —
`true` when both ends carry the *inverted* preamble. Both ends must match, and in
a common polarity, which is what resolves the 180-degree carrier ambiguity.

This is the form every FEC-bearing decoder wants. The default `pack_buffer` path
would shift a window-wide integer once per symbol to expose a preamble sitting at
two known offsets, and those decoders consume their payload as soft symbols
anyway, so the packed window is never needed at all — hence no
`packed_buffer_type` method for Galileo I/NAV, E5a or E6-B, nor for BeiDou B2b.
Measured per `try_sync` call on a full window: 195 ns → 18 ns for I/NAV's
260-symbol window and 521 ns → 18 ns for E5a's 512-symbol one, which is 0.26 →
0.09 and 0.57 → 0.07 µs per symbol end to end. For the two 1000 sps decoders it
is a 1016-bit shift per symbol that never happens at all.

The head is tested before the tail is packed: on noise the head matches with
probability `2^-(preamble_length-1)`, so packing the tail first is wasted work on
essentially every symbol.

!!! note "Pass the lengths, not the state"

    This runs once per *symbol* — the hottest path in the package. The two
    lengths are taken as arguments rather than read from `state.constants` so
    that callers can hand over their own module-level constants: those fold at
    compile time and let `pack_bits_msb_first` unroll, which reading the equivalent
    struct fields does not. Measured at 1000 sps, sourcing them from `constants`
    instead costs about 14 % of the whole per-symbol budget.
"""
@inline function find_preamble_in_deque(
    deque::CircularDeque{Float32},
    preamble::Unsigned,
    preamble_length::Int,
    sync_unit_length::Int,
)
    pattern = UInt64(preamble)
    inverted = pattern ⊻ ((UInt64(1) << preamble_length) - UInt64(1))
    head = pack_bits_msb_first(deque, 1, preamble_length)
    upright = head == pattern
    (upright || head == inverted) || return nothing
    tail = pack_bits_msb_first(deque, sync_unit_length + 1, preamble_length)
    tail == (upright ? pattern : inverted) || return nothing
    return !upright
end

"""
    copy_soft_window!(dst, deque, offset, num_symbols, polarity_flipped) -> dst

Copy `num_symbols` soft symbols out of `deque`, starting at 1-based
`offset + 1`, resolving the 180-degree polarity as it goes: an inverted symbol
stream is a negated one, so the LLR magnitudes (and hence the confidences the
FEC decoders weigh) survive unchanged.

`dst` is the caller's long-lived scratch buffer from its cache, so the copy does
not allocate. Every FEC-bearing decoder here needs exactly this step between
"sync found" and "hand the payload to the decoder".
"""
function copy_soft_window!(
    dst::AbstractVector{Float32},
    deque::CircularDeque{Float32},
    offset::Int,
    num_symbols::Int,
    polarity_flipped::Bool,
)
    sign = polarity_flipped ? -1.0f0 : 1.0f0
    @inbounds for i = 1:num_symbols
        dst[i] = sign * deque[offset+i]
    end
    return dst
end

"""
    complement_buffer_if_necessary(state, buffer) -> (state, resolved_buffer)

If the newest preamble in `buffer` is the *inverted* preamble, return the
state flagged `is_shifted_by_180_degrees = true` together with the
complemented buffer; otherwise return the state flagged `false` with `buffer`
unchanged. The polarity-resolved buffer is returned as a value for the
per-signal `decode_sync_unit` to consume.
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

Drop the consumed `sync_unit_length` oldest soft symbols from the
deque, keeping the trailing `preamble_length` symbols as the leading
preamble of the next subframe.

Drops at most `length(deque)` symbols: a `decode_sync_unit` hook may
reset the decoder mid-frame (e.g. GPS L1C-D on a TOI discontinuity), which
empties the buffer. Without the clamp the unconditional drain in `decode`
would `popfirst!` an empty `CircularDeque` and throw.
"""
function drain_after_sync!(state::GNSSDecoderState)
    deque = soft_buffer(state)
    n_drop = min(state.constants.sync_unit_length, length(deque))
    for _ = 1:n_drop
        popfirst!(deque)
    end
    state
end

"""
$(TYPEDSIGNATURES)

Decode the soft bits of a GNSS navigation message, **overwriting** the storage `state`
holds, and return the updated state. Allocates nothing.

`decode!` writes straight into the buffers that `state` was constructed with —
the soft-symbol buffer, the FEC scratch, the voting tallies in `state.cache`,
and every container referenced from `state.raw_data` and `state.data` (almanac
stores, health tables, HAS masks, ...). Those are all sized once, when the
decoder state is constructed, so a streaming receiver pays no allocation per
symbol:

```julia
state = GPSL1CADecoderState(25)
for chunk in soft_symbol_chunks
    state = decode!(state, chunk, length(chunk))
end
```

!!! warning "Earlier states are overwritten"

    The returned state shares its storage with `state`, and later calls keep
    overwriting it. Treat the return value as *the* live decoder and do not use
    `state` (or anything read out of its containers, such as `state.data.almanacs`)
    after the call. Take a [`copy`](@ref Base.copy(::GNSSDecoderState)) first if
    a snapshot is needed.

Processes incoming soft bits from a GNSS signal, detecting preambles and
decoding sync units to extract navigation data. The function
handles both normal and 180-degree phase-shifted signals automatically.

# Soft-bit convention

`soft_bits` are what a tracking loop hands over, such as TrackingLoops'
`get_soft_bits`: one per navigation bit as tracked, that is per channel symbol,
before any error-correction decoding — the bits GNSSSignals' `get_bit_rate` counts.
On a signal with forward error correction they are not the message's bits: the
decoded message bits exist only inside this package.

`soft_bits` is an `AbstractVector{<:Real}`; `Float32` is canonical. The sign
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
    the magnitude scale. The LDPC decodes (GPS L1C-D and the BeiDou B-CNAV
    family: B1C, B2a, B2b) are flooding sum-product, which *is* scale-sensitive,
    so there the magnitudes should be confidence-weighted on a roughly LLR-like
    scale (`≈ 2·r/σ²`) for best performance at marginal SNR — but still need not
    be normalized to a fixed range.

Glue from a tracking loop: feed its `get_soft_bits` (TrackingLoops' or
`Tracking.jl`'s polarity-corrected, amplitude-weighted soft bits) for every
signal. See `CONTEXT.md` for the full glossary.

# Arguments

  - `state::GNSSDecoderState`: Current decoder state
  - `soft_bits::AbstractVector{<:Real}`: Soft bits to consume, oldest first
  - `num_bits::Int=length(soft_bits)`: Number of leading symbols to process;
    pass an explicit count when the vector contains unused trailing entries.

# Returns

  - `GNSSDecoderState`: Updated decoder state with newly decoded data, sharing
    (and having overwritten) `state`'s storage

# See Also

  - [`GNSSDecoderState`](@ref): The state structure being updated
  - [`reset_decoder_state!`](@ref): Reset the decoder after a signal loss
  - [`is_sat_healthy`](@ref): Check satellite health after decoding
"""
function decode!(
    state::GNSSDecoderState,
    soft_bits::AbstractVector{<:Real},
    num_bits::Int = length(soft_bits),
)
    num_bits >= 0 || throw(ArgumentError("num_bits must be nonnegative"))
    num_bits <= length(soft_bits) ||
        throw(ArgumentError("num_bits exceeds length(soft_bits)"))
    num_bits_since_time_reference = state.num_bits_since_time_reference
    next = 1
    while next <= num_bits
        next, buffer, num_bits_since_time_reference = push_until_sync!(
            state,
            soft_bits,
            next,
            num_bits,
            num_bits_since_time_reference,
        )
        isnothing(buffer) && break
        state = with_symbol_counter(state, num_bits_since_time_reference)
        state, resolved_buffer = complement_buffer_if_necessary(state, buffer)
        state = decode_sync_unit(state, resolved_buffer)
        state = validate_data(state)
        state = drain_after_sync!(state)
        num_bits_since_time_reference = state.num_bits_since_time_reference
    end
    return with_symbol_counter(state, num_bits_since_time_reference)
end

"""
    push_until_sync!(state, soft_symbols, first, last, num_bits_since_time_reference)
        -> (next, buffer, num_bits_since_time_reference)

The per-symbol half of [`decode!`](@ref): push `soft_symbols[first:last]` onto
`state`'s soft-symbol buffer (overwriting its oldest entries) one at a time,
counting them in `num_bits_since_time_reference` (the state's running counter, or
`nothing` before the first sync), until `try_sync` finds a sync. Returns the
index of the next unconsumed symbol, the sync buffer (`nothing` if the symbols
ran out first) and the updated count.

It never rebuilds `state`, which is what keeps it cheap: `state` holds the whole
decoded data inline, and a loop that reassigns it carries a copy of it through
every iteration. No `try_sync` hook reads the symbol counter, so it is written
back into the state only once a sync is to be decoded.
"""
function push_until_sync!(
    state::GNSSDecoderState,
    soft_symbols::AbstractVector{<:Real},
    first::Int,
    last::Int,
    num_bits_since_time_reference::Union{Nothing,Int},
)
    for i = first:last
        push_soft_symbol!(state, soft_symbols[i])
        if num_bits_since_time_reference !== nothing
            num_bits_since_time_reference += 1
        end
        if is_enough_buffered_symbols_to_decode(state)
            buffer = try_sync(state)
            isnothing(buffer) || return (i + 1, buffer, num_bits_since_time_reference)
        end
    end
    return (last + 1, nothing, num_bits_since_time_reference)
end

# Split on the counter being set, so each rebuild has a concretely typed keyword
# (a `Union` keyword value takes the allocating keyword path on Julia 1.10).
with_symbol_counter(state::GNSSDecoderState, num_bits_since_time_reference) =
    num_bits_since_time_reference === nothing ?
    GNSSDecoderState(state; num_bits_since_time_reference = nothing) :
    GNSSDecoderState(state; num_bits_since_time_reference)

"""
    copy(state::GNSSDecoderState) -> GNSSDecoderState

Independent copy of a decoder state: its soft-symbol buffer, cache and every
container referenced from `raw_data` and `data` are copied (see `duplicate`),
so [`decode!`](@ref) on the copy overwrites nothing `state` can see, and vice
versa. The copy keeps every buffer's capacity, so it decodes without allocating
as well.
"""
Base.copy(state::GNSSDecoderState) = GNSSDecoderState(
    state.prn,
    duplicate(state.raw_data),
    duplicate(state.data),
    state.constants,
    duplicate(state.cache),
    state.num_bits_since_time_reference,
    state.is_shifted_by_180_degrees,
)

"""
    reset_decoder_state!(state::GNSSDecoderState) -> GNSSDecoderState

Reset a decoder after a signal loss or reacquisition, **overwriting** the
buffers `state` holds (the soft-symbol buffer is emptied), and return the reset
state. Allocates nothing. What is reset (the time of week, the validated `data`)
and what is kept for a fast recovery is per signal, and documented on each
signal's method. As with [`decode!`](@ref), treat the return value as the live
decoder and do not use `state` afterwards.
"""
function reset_decoder_state! end

"""
    reset_decoder_state!(state::GNSSDecoderState, prn::Integer) -> GNSSDecoderState

Reset a decoder for another satellite, `prn`, **overwriting** the buffers
`state` holds, and return the reset state. Beyond what the one-argument
`reset_decoder_state!` resets, it clears everything that keeps for a fast
recovery of the same satellite: the decoded `raw_data` and `data`, the
repetition-voting tallies, and any message assembly in flight. The polarity is
unknown again. The preallocated buffers and the decoding scratch are kept, so
the result decodes `prn` as a newly constructed decoder of the same signal
would, and the reset allocates nothing.

This is how a receiver moves a channel's decoder to another satellite instead
of constructing a new state, whose buffers would be allocated anew. Everything a
decoder does differently per satellite follows `state.prn`, so it changes with
it: BeiDou B1I and B3I switch between the D1 and D2 formats and their symbol
rates, and B1C syncs on the new satellite's code word. Galileo E6-B retains
its HAS assembly, masks and corrections: these belong to the shared service,
whose messages can be assembled from pages of any satellite. It clears only
the channel's soft-symbol buffer, time counter and polarity.

As with [`decode!`](@ref), treat the return value as the live decoder and do not
use `state` afterwards.
"""
function reset_decoder_state!(state::GNSSDecoderState{D}, prn::Integer) where {D}
    state = reset_decoder_state!(state)
    GNSSDecoderState(
        state;
        prn,
        raw_data = D(),
        data = D(),
        cache = clear_satellite_cache(state.cache),
        num_bits_since_time_reference = nothing,
        is_shifted_by_180_degrees = false,
    )
end

# The part of a decoder's cache that belongs to the satellite it decoded and that the
# one-argument `reset_decoder_state!` keeps for a fast recovery: voting tallies and
# message assembly in flight. Every cache states its own method, beside the cache,
# returning it unchanged where it holds none; there is no fallback, so a new cache
# without one fails loudly rather than carrying its satellite's state over.
function clear_satellite_cache end
