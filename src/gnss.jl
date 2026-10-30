abstract type AbstractGNSSConstants end
abstract type AbstractGNSSData end
abstract type AbstractGNSSCache end

"""
    AbstractGPSData <: AbstractGNSSData

Abstract supertype for the decoded navigation data of a signal transmitted by
the GPS constellation, e.g. `GPSLNAVData`, `GPSCNAVData`.

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
    AbstractGPSCNAVData <: AbstractGPSData

Abstract supertype for the modern GPS civil navigation messages: CNAV
(`GPSCNAVData`, on L2C and L5) and CNAV-2 (`GPSCNAV2Data`, on L1C).

What the two messages share — and legacy LNAV (`GPSLNAVData`) does not — is a
full-width 13-bit week number and the quasi-Keplerian ephemeris broadcast as
deltas off the fixed reference values `A_REF` and `Ω̇_REF` (IS-GPS-200N
§30.3.3.1.1, IS-GPS-800J §3.5.3), which [`get_broadcast_ephemeris`](@ref)
applies for both messages through this type. The GPS counterpart of
[`AbstractBeiDouCNAVData`](@ref).
"""
abstract type AbstractGPSCNAVData <: AbstractGPSData end

"""
    AbstractGalileoData <: AbstractGNSSData

Abstract supertype for the decoded navigation data of a signal transmitted by
the Galileo constellation, e.g. `GalileoINAVData` (E1-B, E5b-I), `GalileoFNAVData`,
`GalileoE6BData`.

The Galileo counterpart to [`AbstractGPSData`](@ref): it carries the facts every
Galileo signal's data shares, whatever that signal broadcasts. The
health-status and positioning-readiness checks genuinely differ per signal and
stay on the concrete data types.

Note that "every Galileo signal" is not "every Galileo signal carries an
ephemeris" — `GalileoE6BData` holds HAS corrections *to* another signal's
navigation data and has no orbital fields at all. The ephemeris and clock
completeness checks therefore live one level down, on the internal
`AbstractGalileoEphemerisData`.
"""
abstract type AbstractGalileoData <: AbstractGNSSData end

"""
    AbstractGalileoEphemerisData <: AbstractGalileoData

Abstract supertype for the Galileo signals whose navigation message carries an
ephemeris and clock of its own: `GalileoINAVData` (E1-B, E5b-I) and
`GalileoFNAVData` (E5a-I).

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
the BeiDou constellation, e.g. `BeiDouDNAVData` (B1I/B3I), `BeiDouCNAV1Data`,
`BeiDouCNAV2Data`, `BeiDouCNAV3Data`.

The BeiDou counterpart to [`AbstractGPSData`](@ref) and
[`AbstractGalileoData`](@ref): it carries the constellation-level facts every
BeiDou signal's data shares, so they can be stated once via subtype dispatch.
Genuinely per-signal facts (the message-set completeness checks, the
health-flag selection in `is_sat_healthy`) stay defined on the concrete data
types.
"""
abstract type AbstractBeiDouData <: AbstractGNSSData end

"""
    AbstractBeiDouCNAVData <: AbstractBeiDouData

Abstract supertype for the BDS-3 B-CNAV navigation messages: B-CNAV1
(`BeiDouCNAV1Data`, on B1C), B-CNAV2 (`BeiDouCNAV2Data`, on B2a) and B-CNAV3
(`BeiDouCNAV3Data`, on B2b).

What the three messages share — and the legacy D1/D2 message (`BeiDouDNAVData`,
on B1I/B3I) does not — is the quasi-Keplerian ephemeris broadcast as `ΔA` off a
`sat_type`-selected reference semi-major axis with an outright `Ω̇`, the
broadcast `sat_type` orbit class itself ([`get_orbit_class`](@ref) dispatches on
this type), and the BDGIM ionospheric coefficient set. The BeiDou counterpart
of [`AbstractGPSCNAVData`](@ref).
"""
abstract type AbstractBeiDouCNAVData <: AbstractBeiDouData end

# Physical constants common to every GNSS handled here. The constants each ICD
# fixes for its own ephemeris algorithm — the Earth gravitational parameter μ,
# the relativistic correction F that follows from it and, for BeiDou, the Earth
# rotation rate — are *not* here: GPS uses the WGS-84 values, BeiDou the BDCS
# values and Galileo its own, defined alongside each constellation.
#
# `GNSS_PI` is deliberately the truncated value the ICDs fix for the
# semicircle→radian scaling, *not* `Base.π`: every ICD here states the same
# 3.1415926535898 (IS-GPS-200 Table 20-IV; Galileo OS SIS ICD Table 68;
# BDS-SIS-ICD-B1I-3.0 and B3I-1.0 Table 5-11, B1C-1.0 and B2a-1.0 Table 7-9,
# B2b-1.0 Table 7-7) — broadcast angular quantities must be scaled with exactly this
# value to reproduce the transmitted numbers bit-for-bit.
const GNSS_PI = 3.1415926535898
const EARTH_ROTATION_RATE = 7.2921151467e-5  # rad/s (WGS-84 and the Galileo ICD agree)

"""
$(TYPEDEF)

Generic decoder state for GNSS signal decoding. This parametric struct holds all state
required for decoding navigation messages from GNSS satellites.

The struct itself is immutable; per-field reconstruction works via the keyword
constructor. Every buffer a decoder needs is allocated once, by the per-signal
constructor, and owned by the state: the soft-symbol `CircularDeque{Float32}`
(capacity `sync_unit_length + preamble_length`), the FEC scratch and
voting tallies in the `cache`, and the containers (almanac stores, health
tables, ...) that `raw_data` and `data` reference, preallocated in the cache's
`DataStorage`.

Those buffers are **overwritten** by [`decode!`](@ref), which is what makes it
allocation-free: the state it returns shares them with the state passed in, so
treat the returned value as *the* live state and do not use the earlier one;
take a [`copy`](@ref Base.copy(::GNSSDecoderState)) if a snapshot is needed. The transient packed-bit buffer used for
preamble matching is **not** stored here; it is computed as a local value at
sync time and threaded through the sync path (see `pack_buffer` / `try_sync`).

Of the fields, `prn`, `raw_data`, `data`, `num_bits_since_time_reference` and
`is_shifted_by_180_degrees` are part of the API: a receiver reads them, and
`num_bits_since_time_reference` together with
[`get_broadcast_time_of_week`](@ref) and
`get_bit_rate` dates the latest symbol (see `get_broadcast_time_of_week`). A tracking
loop that counts symbols itself may move the counter with the keyword
constructor, `GNSSDecoderState(state; num_bits_since_time_reference = n)`. The
`constants` and the `cache` are internal.

# Construction

`GNSSDecoderState(signal, prn)` creates the decoder for satellite `prn` on a
GNSSSignals signal: `GNSSDecoderState(GPSL1CA(), 1)` is `GPSL1CADecoderState(1)`.
Every decoded signal has both forms: GPS L1 C/A, L1C-D, L2C-M and L5-I; Galileo E1-B
(both the default BOC(1,1) and full CBOC replicas), E5a-I, E5b-I and E6-B; BeiDou B1I, B3I, B1C-D,
B2a-I and B2b-I. A signal without a decoder, such as a dataless pilot, is a
`MethodError`.

# Type Parameters

  - `D<:AbstractGNSSData`: The data type holding decoded navigation message fields
  - `C<:AbstractGNSSConstants`: Constants specific to the GNSS system (e.g., preamble, timing)
  - `CA<:AbstractGNSSCache`: Cache for intermediate decoding state (carries the soft-symbol buffer)

# Fields

$(TYPEDFIELDS)

# See Also

  - [`GPSL1CADecoderState`](@ref): Constructor for GPS L1 C/A decoder state
  - [`GalileoE1BDecoderState`](@ref): Constructor for Galileo E1B decoder state
  - [`decode!`](@ref): Decode soft symbols into this state's buffers, allocation-free
  - [`reset_decoder_state!`](@ref): Reset decoder state after signal loss
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
    Number of channel symbols since the epoch stamped by `data`, or `nothing` without a time anchor
    """
    num_bits_since_time_reference::Union{Nothing,Int} = 0
    """
    Whether the signal phase is inverted by 180 degrees
    """
    is_shifted_by_180_degrees::Bool = false
end

@inline function GNSSDecoderState(
    state::GNSSDecoderState;
    prn = state.prn,
    raw_data = state.raw_data,
    data = state.data,
    cache = state.cache,
    num_bits_since_time_reference = state.num_bits_since_time_reference,
    is_shifted_by_180_degrees = state.is_shifted_by_180_degrees,
)
    GNSSDecoderState(
        Int(prn),
        raw_data,
        data,
        state.constants,
        cache,
        num_bits_since_time_reference,
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
@generated function fields_equal(a::T, b::T) where {T}
    # Unrolled at compile time: a runtime loop over `fieldnames(T)` reads each
    # field through a non-constant name, so every `==` would dispatch
    # dynamically on the union of all field types — which `juliac --trim`
    # rejects.
    comparisons =
        [:(getfield(a, $i) == getfield(b, $i) || return false) for i = 1:fieldcount(T)]
    return quote
        $(comparisons...)
        return true
    end
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
        has no notion of "now", so the consumer must still reject stale
        ephemerides. [`get_broadcast_ephemeris_validity`](@ref) reports the
        absolute reference epochs and the fit duration where broadcast.

      - **Second-order corrections.** Klobuchar / BDGIM ionosphere and UTC
        parameters are intentionally excluded: they are broadcast far less
        often, so apply them when present and treat `nothing` as zero rather
        than waiting for them. Group delay and inter-signal corrections
        (`T_GD`, `ISC_*`) split on one question — *does this signal's own
        correction ride in the block the gate already waits for?*

        Where it does, it is required, because it is a metre-scale bias rather
        than a refinement and costs nothing to wait for. BeiDou B1C's
        `T_GD_B1Cp` and `ISC_B1Cd` are in the same CRC-protected subframe 2 as
        the ephemeris and clock, B2b's `T_GD_B2bI` is in the same MT30 as the
        clock, and B1I/B3I's `T_GD1` is in the same D1 subframe 1 / D2 page 1
        as the clock, so all are gated on here.

        Where it does not, it is excluded, because requiring it would mean
        waiting for a message the ICD declines to schedule: BeiDou B2a's
        `T_GD_B2ap` / `ISC_B2ad` are in MT30 alone while the clock is in
        MT30-34, and the GPS CNAV `T_GD` and `ISC_*` (message type 30) and the
        CNAV-2 ISCs of the other civil signals (subframe 3 page 1) likewise
        arrive on their own cadence. CNAV-2's own `T_GD`, `ISC_L1CD` and
        `ISC_L1CP` ride in subframe 2 with the ephemeris.

        Corrections for a band this decoder is not on — B1C's `T_GD_B2ap`,
        B1I/B3I's `T_GD2` (there is no B2I signal here to range on, and none
        planned), the GPS `ISC_L5Q5` on an L1 fix — are never required either
        way. Such a field is still decoded and published; it is only this
        readiness gate that ignores it.

      - **Alert flag.** `is_sat_healthy` reflects the broadcast health bits
        only; a receiver that wants to honour the L1 C/A alert flag (or
        equivalent) must check it separately.
"""
has_positioning_data(state::GNSSDecoderState) = has_positioning_data(state.data)

# ---- Shared decoder primitives ----------------------------------------------
#
# Signal-agnostic primitives used by more than one signal decoder.

# Packed-word integer types used by more than one signal decoder (`BitIntegers`
# widths are global once defined, so each is defined exactly once, here).
#
#   - `UInt288`: the 220 bits a Galileo I/NAV page pair runs the CRC-24Q over,
#     the 274-bit GPS L1C-D subframe 3 and the 264-bit BeiDou B1C subframe 3.
#   - `UInt320`: a GPS L1 C/A subframe, a GPS L1C-D subframe, a GPS CNAV message
#     (L5I / L2C) and a BeiDou D1/D2 subframe — 300 data bits plus up to 8
#     trailing sync bits.
#   - `UInt512`: the 486-bit information word of a Galileo E6-B C/NAV page and
#     of a BeiDou B2b message.
#   - `UInt600`: the 600-bit subframe 2 of GPS L1C-D and of BeiDou B1C.
#
# `BitIntegers.@define_integers` also defines the signed companions `Int288`,
# `Int320`, `Int512` and `Int600`.
BitIntegers.@define_integers 288
BitIntegers.@define_integers 320
BitIntegers.@define_integers 512
BitIntegers.@define_integers 600

"""
    increment_voting(old_vote, max_vote) -> Int

Raise a broadcast-repetition vote by one, saturating at `max_vote`.

Two decoders promote `raw_data` to `data` by repetition voting rather than by an
issue-of-data match — GPS L1 C/A and the BeiDou D1/D2 legacy message, which has
no IODs at all — and both count with this.
"""
increment_voting(old_vote, max_vote) = min(max_vote, old_vote + 1)

"""
    elapsed_symbols(anchor, now) -> Union{Nothing,Int64}

Symbols counted between the decoder's symbol-counter reading `anchor` and its
current reading `now` (both `num_bits_since_time_reference` values), or
`nothing` if either reading is missing.
"""
elapsed_symbols(anchor, now) = isnothing(anchor) || isnothing(now) ? nothing : now - anchor
