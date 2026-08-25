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
the Galileo constellation, e.g. `GalileoINAVData` (E1-B, E5b-I), `GalileoE5aData`,
`GalileoE6BData`.

The Galileo counterpart to [`AbstractGPSData`](@ref): it carries the facts every
Galileo signal's data shares, whatever that signal broadcasts. The
health-status and positioning-readiness checks genuinely differ per signal and
stay on the concrete data types.

Note that "every Galileo signal" is not "every Galileo signal carries an
ephemeris" — `GalileoE6BData` holds HAS corrections *to* another signal's
navigation data and has no orbital fields at all. The ephemeris and clock
completeness checks therefore live one level down, on
[`AbstractGalileoEphemerisData`](@ref).
"""
abstract type AbstractGalileoData <: AbstractGNSSData end

"""
    AbstractGalileoEphemerisData <: AbstractGalileoData

Abstract supertype for the Galileo signals whose navigation message carries an
ephemeris and clock of its own: `GalileoINAVData` (E1-B, E5b-I) and
`GalileoE5aData` (E5a-I).

This exists so `is_ephemeris_decoded` and `is_clock_correction_decoded` — which
check the same orbital and clock fields for I/NAV and F/NAV, and so are stated
once (see `src/galileo/galileo.jl`) — are dispatched on a type that actually has
those fields. Declaring them on `AbstractGalileoData` instead would put a method
on the supertype that raises a `FieldError` for `GalileoE6BData`, whose C/NAV
message broadcasts corrections rather than an ephemeris; a future
corrections-only or almanac-only Galileo signal would inherit the same trap.
"""
abstract type AbstractGalileoEphemerisData <: AbstractGalileoData end

"""
    AbstractBeiDouData <: AbstractGNSSData

Abstract supertype for the decoded navigation data of a signal transmitted by
the BeiDou constellation, e.g. `BeiDouDNAVData` (B1I/B3I), `BeiDouB1CData`,
`BeiDouB2aData`, `BeiDouB2bData`.

The BeiDou counterpart to [`AbstractGPSData`](@ref) and
[`AbstractGalileoData`](@ref): it carries the constellation-level facts every
BeiDou signal's data shares, so they can be stated once via subtype dispatch.
Genuinely per-signal facts (the message-set completeness checks, the
health-flag selection in `is_sat_healthy`) stay defined on the concrete data
types.
"""
abstract type AbstractBeiDouData <: AbstractGNSSData end

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

# Every constellation here counts time as a week number plus a time of week,
# so the week length is shared rather than owned by whichever signal happened
# to need it first.
const SECONDS_PER_WEEK = 604_800

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
    """
    Pseudo-Random Noise code identifier for the satellite
    """
    prn::Int
    """
    Partially decoded navigation data (not yet validated)
    """
    raw_data::D
    """
    Validated navigation data ready for use
    """
    data::D
    """
    System-specific constants (preamble, timing parameters)
    """
    constants::C
    """
    Cache for intermediate decoding state (holds the soft-symbol `CircularDeque{Float32}`)
    """
    cache::CA
    """
    Number of symbols received after the last valid synchronization sequence, or `nothing` if not yet synchronized
    """
    num_bits_after_valid_syncro_sequence::Union{Nothing,Int} = 0
    """
    Whether the signal phase is inverted by 180 degrees
    """
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
# inside `GalileoINAVData`) falls back to `===`. Compare field-by-field so that
# two states with equal-but-not-identical contents are considered equal.
function Base.:(==)(a::GNSSDecoderState, b::GNSSDecoderState)
    typeof(a) === typeof(b) || return false
    return fields_equal(a, b)
end

"""
    fields_equal(a::T, b::T) -> Bool

Internal helper: structural equality field by field.

Julia's default `==` on a struct falls back to `===`, which is *reference*
equality for any field of mutable type — a `Vector`, a `Dictionary`, a
`CircularDeque`. Every container in this package that holds one therefore has to
define `==` explicitly, and this is that definition, written once.
"""
function fields_equal(a::T, b::T) where {T}
    for f in fieldnames(T)
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
    for i = 1:length(a)
        a[i] == b[i] || return false
    end
    return true
end

"""
$(TYPEDSIGNATURES)

Report whether a decoder has recovered the minimum navigation data a
positioning engine needs from this satellite: a time of week, a full ephemeris
(orbit) set, and the SV clock-correction polynomial (plus, on the signals that
carry it in the same required set, the broadcast week number and single-band
group delay). Dispatches on the validated [`data`](@ref GNSSDecoderState) field,
so it only becomes `true` once the required message set has passed CRC/parity
**and** the cross-subframe issue-of-data consistency check that promotes
`raw_data` to `data`.

This is the readiness gate a receiver (e.g. `PositionVelocityTime.jl`) should
pair with [`is_sat_healthy`](@ref): whenever this returns `true`, the health
field `is_sat_healthy` inspects is guaranteed to have been decoded, so the two
can be checked together without a separate `nothing` guard.

!!! note "What this deliberately does *not* gate on"

    A `true` here means the *data set* is complete and self-consistent — it is a
    necessary condition for using the SV in a fix, not a blanket guarantee that
    no further judgement is required:

      - **Ephemeris freshness.** Only presence is checked, not age. The decoder
        has no notion of "now", so the consumer must still reject ephemerides
        outside their fit interval (`fit_interval` / `t_oe` age).
      - **Second-order corrections.** Group delay / inter-signal corrections
        (`T_GD`, `ISC_*`) beyond the single required band, Klobuchar
        ionosphere, and UTC parameters are intentionally excluded because they
        are broadcast far less often; apply them when present and treat
        `nothing` as zero rather than waiting for them.
      - **Alert flag.** `is_sat_healthy` reflects the broadcast health bits
        only; a receiver that wants to honour the L1 C/A alert flag (or
        equivalent) must check it separately.
"""
is_decoding_completed_for_positioning(state::GNSSDecoderState) =
    is_decoding_completed_for_positioning(state.data)

# ---- Signal metadata --------------------------------------------------------
#
# A decoder state is 1:1 with the signal whose navigation message it
# demodulates, so everything GNSSSignals already knows about that signal — what
# it is called, which constellation and band it belongs to, how fast its
# navigation message runs, which time scale its week numbers and times of week
# are counted in — is knowable from the state alone. Restating any of it here
# would be a second source of truth to keep in sync, so instead each signal file
# states the *one* fact that is genuinely ours (the constants → signal mapping,
# `get_signal_type`) and every accessor below forwards through it.

"""
    get_signal_type(state::GNSSDecoderState) -> Type{<:AbstractGNSSSignal}
    get_signal_type(constants::AbstractGNSSConstants) -> Type{<:AbstractGNSSSignal}

Get the GNSSSignals signal type whose navigation message a decoder demodulates,
e.g. `GPSL1CA` for a [`GPSL1CADecoderState`](@ref). This is the single mapping
from this package's decoders back into GNSSSignals; every signal accessor listed
under "Signal metadata" in the docs is forwarded through it, and it is the entry
point for anything not forwarded (`get_signal_type(state)` accepts every
GNSSSignals accessor that takes a signal type).

Returns the *type*, not an instance: constructing a signal builds its spreading
code matrix and SIMD lookup table, which a decoder — working purely in the
symbol domain — never needs. Being a type, it also folds to a compile-time
constant, so the forwarded accessors cost nothing at run time.

Stated per signal file, dispatched on the constants type rather than the data
type, because the constants are what tell apart decoders that share a data
container: GPS L5-I and L2C-M both decode into a `GPSCNAVData` but are distinct
signals (`GPSCNAVConstants{:GPSL5I}` vs `GPSCNAVConstants{:GPSL2CM}`).

The mapping names the data-bearing component of a signal pair, since that is
what carries the navigation message: `GPSL2CM` (not the `GPSL2CL` pilot) and
`GalileoE5aI` (not the `GalileoE5aQ` pilot). A decoder built from an
approximation of a signal reports the signal it approximates — the E1B BOC(1,1)
decoder state is a `GalileoE1B`, since the approximation is a
tracking/acquisition concern and the I/NAV stream it decodes is identical.

# Examples

```julia
using GNSSDecoder, GNSSSignals
get_signal_type(GPSL1CADecoderState(1))            # GPSL1CA
get_signal_type(GPSL2CMDecoderState(1))            # GPSL2CM
get_code_length(get_signal_type(GPSL1CADecoderState(1)))  # 1023
```

# See Also

  - [`GNSSDecoderState`](@ref): The state this is queried from
"""
@inline get_signal_type(state::GNSSDecoderState) = get_signal_type(state.constants)

# The GNSSSignals accessors a decoder state can answer for itself, forwarded
# through `get_signal_type`: identity (`get_signal_id` / `get_signal_name`,
# `get_constellation_id` / `get_constellation_name`, `get_band` / `get_band_id` /
# `get_band_name`), the navigation-message symbol rate (`get_data_frequency`),
# and the time scale the decoded week numbers and times of week are referenced
# to (`get_time_system` and friends — the epoch `get_system_start_time` and the
# TAI offset `get_tai_offset` are what turn a decoded WN/TOW pair into an
# absolute instant, and the Galileo epoch is 13 s shy of a UTC midnight, so
# sourcing it rather than restating it matters).
#
# These extend GNSSSignals' own functions rather than introducing GNSSDecoder
# names, so a caller already doing `using GNSSSignals` needs no new import and
# can query a decoder state exactly as it queries a signal. Accessors that
# describe the *code* (`get_code_length`, `get_code_frequency`,
# `get_carrier_phase_offset`, …) are deliberately not forwarded: they belong to
# acquisition and tracking, not to a symbol-domain decoder, and remain one
# `get_signal_type(state)` away.
for accessor in (
    :get_signal_id,
    :get_signal_name,
    :get_constellation_id,
    :get_constellation_name,
    :get_band,
    :get_band_id,
    :get_band_name,
    :get_data_frequency,
    :get_time_system,
    :get_time_system_id,
    :get_time_system_name,
    :get_system_start_time,
    :get_tai_offset,
)
    @eval @inline GNSSSignals.$accessor(state::GNSSDecoderState) =
        GNSSSignals.$accessor(get_signal_type(state))
end

"""
Soft-symbol buffer accessor — the per-signal cache stores it as `soft_buffer`.
"""
soft_buffer(state::GNSSDecoderState) = state.cache.soft_buffer

"""
Number of soft symbols currently buffered.
"""
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
    pack_soft_bits(deque, offset, len) -> UInt64

Hard-slice `len` soft symbols of `deque` starting at 1-based `offset` into a
`UInt64`, oldest symbol at the MSB — the same bit order [`pack_buffer`](@ref)
produces, but for a slice rather than the whole window.

Sync only ever inspects a handful of known positions (the preamble at either
end, and any unencoded header field), so a decoder whose window is long and
whose symbol rate is high can read exactly those bits instead of repacking
the entire window on every symbol. `len` must be at most 64.
"""
function pack_soft_bits(deque::CircularDeque{Float32}, offset::Int, len::Int)
    word = UInt64(0)
    @inbounds for i = 0:(len-1)
        word = (word << 1) | (hard_slice(deque[offset+i]) ? UInt64(1) : UInt64(0))
    end
    return word
end

"""
    pack_soft_codeword(deque, offset, len) -> UInt64

Hard-slice `len` soft symbols of `deque` starting at 1-based `offset` into a
`UInt64` with **bit 0 = the first symbol** — the opposite end from
[`pack_soft_bits`](@ref), and the order the codeword tables use
(`BCH_TOI_CODEWORDS`, `b1c_prn_codeword`, `b1c_soh_codeword`; cf.
[`soft_to_hard_codeword`](@ref), which is this over an iterable).

The two BCH-synchronised signals — GPS L1C-D on its TOI codeword and BeiDou B1C
on its PRN + SOH pair — compare against those tables once per symbol, at both
ends of the window, so this reads the deque in place rather than materialising a
slice. `len` must be at most 64.
"""
function pack_soft_codeword(deque::CircularDeque{Float32}, offset::Int, len::Int)
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
        (buffer >> constants.syncro_sequence_length) & mask == ~constants.preamble & mask
end

"""
    find_preamble_in_deque(deque, preamble, preamble_length, syncro_sequence_length)
        -> Union{Nothing,Bool}

Soft-domain counterpart of [`find_preamble`](@ref): match the preamble at both
ends of the candidate window by slicing the soft-symbol deque directly, rather
than hard-slicing the whole window into a packed integer first.

Returns `nothing` when there is no sync, and otherwise the resolved polarity —
`true` when both ends carry the *inverted* preamble. The rule is `find_preamble`'s:
both ends must match, and in a common polarity, which is what resolves the
180-degree carrier ambiguity.

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
    compile time and let `pack_soft_bits` unroll, which reading the equivalent
    struct fields does not. Measured at 1000 sps, sourcing them from `constants`
    instead costs about 14 % of the whole per-symbol budget.
"""
@inline function find_preamble_in_deque(
    deque::CircularDeque{Float32},
    preamble::Unsigned,
    preamble_length::Int,
    syncro_sequence_length::Int,
)
    pattern = UInt64(preamble)
    inverted = pattern ⊻ ((UInt64(1) << preamble_length) - UInt64(1))
    head = pack_soft_bits(deque, 1, preamble_length)
    upright = head == pattern
    (upright || head == inverted) || return nothing
    tail = pack_soft_bits(deque, syncro_sequence_length + 1, preamble_length)
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
    for _ = 1:n_drop
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
    the magnitude scale. The LDPC decodes (GPS L1C-D and the BeiDou B-CNAV
    family: B1C, B2a, B2b) are flooding sum-product, which *is* scale-sensitive,
    so there the magnitudes should be confidence-weighted on a roughly LLR-like
    scale (`≈ 2·r/σ²`) for best performance at marginal SNR — but still need not
    be normalized to a fixed range.

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
state = decode(state, Float32[+1, -1, +1, +1, -1, -1, -1, -1], 8)
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
    for i = 1:num_symbols
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

# ---- Shared decoder primitives ----------------------------------------------
#
# Signal-agnostic primitives used by more than one signal decoder. They live
# here (a shared file included before every signal) rather than in a per-signal
# file so that no signal decoder has to be included after another just to borrow
# them.

# Packed-word integer types shared across signal decoders. `BitIntegers` widths
# are global once defined, so the ones more than one constellation reaches for
# are stated here rather than in whichever signal file happened to need them
# first — otherwise the include order silently becomes load-bearing.
#
#   - `UInt320`: a GPS L1 C/A subframe, a GPS L1C-D subframe, a GPS CNAV message
#     (L5I / L2C) and a BeiDou D1/D2 subframe — 300 data bits plus up to 8
#     trailing sync bits.
#   - `UInt512`: the 486-bit information word of a Galileo E6-B C/NAV page and
#     of a BeiDou B2b message.
#
# `BitIntegers.@define_integers` also defines the signed companions `Int320` and
# `Int512`.
BitIntegers.@define_integers 320
BitIntegers.@define_integers 512

"""
Insert/overwrite `value` keyed by `key` in a (possibly `nothing`) `Dictionary`, returning the updated copy.
"""
function _merge_keyed(dict::Union{Nothing,Dictionary{Int,V}}, key::Int, value::V) where {V}
    out = isnothing(dict) ? Dictionary{Int,V}() : copy(dict)
    set!(out, key, value)
    return out
end
