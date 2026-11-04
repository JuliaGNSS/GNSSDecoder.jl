abstract type AbstractGNSSConstants end
abstract type AbstractGNSSData end
abstract type AbstractGNSSCache end

"""
    AbstractGPSData <: AbstractGNSSData

Abstract supertype for the decoded navigation data of a GPS signal, e.g. `GPSLNAVData`,
`GPSCNAVData`. Carries the constellation-level facts every GPS signal shares, stated once
via subtype dispatch. Per-signal facts (subframe/message-type completeness checks, the
health-bit selection in `is_sat_healthy`) stay on the concrete data types.
"""
abstract type AbstractGPSData <: AbstractGNSSData end

"""
    AbstractGPSCNAVData <: AbstractGPSData

Abstract supertype for the modern GPS civil messages: CNAV (`GPSCNAVData`, on L2C and L5)
and CNAV-2 (`GPSCNAV2Data`, on L1C). The GPS counterpart of
[`AbstractBeiDouCNAVData`](@ref).

Unlike legacy LNAV (`GPSLNAVData`), both carry a full-width 13-bit week number and the
quasi-Keplerian ephemeris as deltas off the fixed references `A_REF` and `Ω̇_REF`
(IS-GPS-200N §30.3.3.1.1, IS-GPS-800J §3.5.3), so a consumer's propagator and week
handling dispatch on this type rather than enumerating the two messages.
"""
abstract type AbstractGPSCNAVData <: AbstractGPSData end

"""
    AbstractGalileoData <: AbstractGNSSData

Abstract supertype for the decoded navigation data of a Galileo signal, e.g.
`GalileoINAVData` (E1-B, E5b-I), `GalileoFNAVData`, `GalileoE6BData`. The Galileo
counterpart to [`AbstractGPSData`](@ref); health-status and positioning-readiness checks
differ per signal and stay on the concrete types.

Not every Galileo signal carries an ephemeris: `GalileoE6BData` holds HAS corrections *to*
another signal's navigation data and has no orbital fields. The ephemeris and clock
completeness checks therefore live one level down, on
[`AbstractGalileoEphemerisData`](@ref).
"""
abstract type AbstractGalileoData <: AbstractGNSSData end

"""
    AbstractGalileoEphemerisData <: AbstractGalileoData

Abstract supertype for the Galileo signals whose message carries its own ephemeris and
clock: `GalileoINAVData` (E1-B, E5b-I) and `GalileoFNAVData` (E5a-I).

`is_ephemeris_decoded` and `is_clock_correction_decoded` (shared by I/NAV and F/NAV, see
`src/galileo/galileo.jl`) dispatch on this type because it has those fields. On
`AbstractGalileoData` they would raise a `FieldError` for `GalileoE6BData` (and any future
corrections-only or almanac-only signal).
"""
abstract type AbstractGalileoEphemerisData <: AbstractGalileoData end

"""
    AbstractBeiDouData <: AbstractGNSSData

Abstract supertype for the decoded navigation data of a BeiDou signal, e.g.
`BeiDouDNAVData` (B1I/B3I), `BeiDouCNAV1Data`, `BeiDouCNAV2Data`, `BeiDouCNAV3Data`. The
BeiDou counterpart to [`AbstractGPSData`](@ref) and [`AbstractGalileoData`](@ref).
Per-signal facts (message-set completeness, health-flag selection in `is_sat_healthy`)
stay on the concrete types.
"""
abstract type AbstractBeiDouData <: AbstractGNSSData end

"""
    AbstractBeiDouCNAVData <: AbstractBeiDouData

Abstract supertype for the BDS-3 B-CNAV messages: B-CNAV1 (`BeiDouCNAV1Data`, on B1C),
B-CNAV2 (`BeiDouCNAV2Data`, on B2a) and B-CNAV3 (`BeiDouCNAV3Data`, on B2b). The BeiDou
counterpart of [`AbstractGPSCNAVData`](@ref).

Unlike legacy D1/D2 (`BeiDouDNAVData`, on B1I/B3I), all three carry the quasi-Keplerian
ephemeris as `ΔA` off a `sat_type`-selected reference semi-major axis with an outright
`Ω̇`, the `sat_type` orbit class itself ([`get_orbit_class`](@ref) dispatches on this
type), and the BDGIM ionospheric coefficients.
"""
abstract type AbstractBeiDouCNAVData <: AbstractBeiDouData end

# Physical constants common to every GNSS handled here. The per-ICD ephemeris constants
# (gravitational parameter μ, relativistic correction F, BeiDou Earth rotation rate) are
# defined alongside each constellation: WGS-84 for GPS, BDCS for BeiDou, Galileo's own.
#
# `GNSS_PI` is deliberately the truncated value the ICDs fix for the semicircle→radian
# scaling, *not* `Base.π`: every ICD states 3.1415926535898 (IS-GPS-200 Table 20-IV;
# Galileo OS SIS ICD Table 68; BDS-SIS-ICD-B1I-3.0 and B3I-1.0 Table 5-11, B1C-1.0 and
# B2a-1.0 Table 7-9, B2b-1.0 Table 7-7), and angles must be scaled with exactly it to
# reproduce the transmitted numbers bit-for-bit.
const GNSS_PI = 3.1415926535898
const EARTH_ROTATION_RATE = 7.2921151467e-5  # rad/s (WGS-84 and the Galileo ICD agree)

"""
$(TYPEDEF)

Generic decoder state for GNSS navigation-message decoding.

The struct is immutable; per-field reconstruction works via the keyword constructor. The
per-signal constructor allocates every buffer once and the state owns them: the
soft-symbol `CircularDeque{Float32}` (capacity `sync_unit_length + preamble_length`), the
FEC scratch and voting tallies in the `cache`, and the containers (almanac stores, health
tables, ...) that `raw_data` and `data` reference, preallocated in the cache's
`DataStorage`.

[`decode!`](@ref) **overwrites** those buffers, which makes it allocation-free: the state
it returns shares them with the state passed in, so use the returned value as *the* live
state and drop the earlier one. Take a [`copy`](@ref Base.copy(::GNSSDecoderState)) for a
snapshot.

# Construction

`GNSSDecoderState(signal, prn)` creates the decoder for satellite `prn` on a GNSSSignals
signal: `GNSSDecoderState(GPSL1CA(), 1)` is `GPSL1CADecoderState(1)`. Every decoded signal
has both forms: GPS L1 C/A, L1C-D, L2C-M and L5-I; Galileo E1-B (default BOC(1,1) and full
CBOC replicas), E5a-I, E5b-I and E6-B; BeiDou B1I, B3I, B1C-D, B2a-I and B2b-I. A signal
without a decoder, such as a dataless pilot, is a `MethodError`.

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
    is_polarity_inverted::Bool = false
end

@inline function GNSSDecoderState(
    state::GNSSDecoderState;
    raw_data = state.raw_data,
    data = state.data,
    cache = state.cache,
    num_bits_since_time_reference = state.num_bits_since_time_reference,
    is_polarity_inverted = state.is_polarity_inverted,
)
    GNSSDecoderState(
        state.prn,
        raw_data,
        data,
        state.constants,
        cache,
        num_bits_since_time_reference,
        is_polarity_inverted,
    )
end

# The default `==` falls back to `===` for mutable fields (the soft-symbol deque, the
# almanac vectors); compare field by field so equal-content states are equal.
function Base.:(==)(a::GNSSDecoderState, b::GNSSDecoderState)
    typeof(a) === typeof(b) || return false
    return fields_equal(a, b)
end

"""
    fields_equal(a::T, b::T) -> Bool

Internal helper: structural equality field by field. Julia's default `==` on a struct is
*reference* equality for mutable fields (`Vector`, `Dictionary`, `CircularDeque`), so every
container here that holds one defines `==` through this.
"""
@generated function fields_equal(a::T, b::T) where {T}
    # Unrolled: a runtime loop over `fieldnames(T)` would dispatch dynamically on the union
    # of all field types, which `juliac --trim` rejects.
    comparisons =
        [:(getfield(a, $i) == getfield(b, $i) || return false) for i = 1:fieldcount(T)]
    return quote
        $(comparisons...)
        return true
    end
end

"""
    deques_equal(a::CircularDeque, b::CircularDeque)

Internal helper: structural equality on `CircularDeque`s of the same element type.
DataStructures.jl defines no `==` for it, and defining one here would be type piracy
(flagged by Aqua). Per-signal caches' `==` calls this.
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

Whether a decoder has recovered the minimum navigation data a positioning engine needs
from this satellite: a time of week, a full ephemeris set and the SV clock-correction
polynomial (plus, on signals that carry them in the same required set, the broadcast week
number and single-band group delay). Dispatches on the validated
[`data`](@ref GNSSDecoderState), so it becomes `true` only once the required message set
has passed CRC/parity **and** the cross-subframe issue-of-data check that promotes
`raw_data` to `data`.

This is the readiness gate a receiver (e.g. `PositionVelocityTime.jl`) should pair with
[`is_sat_healthy`](@ref): when it is `true`, the health field `is_sat_healthy` inspects
has been decoded, so no separate `nothing` guard is needed.

!!! note "What this deliberately does *not* gate on"

    A `true` here means the *data set* is complete and self-consistent: necessary for using
    the SV in a fix, not sufficient.

      - **Ephemeris freshness.** Only presence is checked, not age. The decoder has no
        notion of "now", so the consumer must reject ephemerides outside their fit
        interval (`fit_interval` / `t_0e` age).

      - **Second-order corrections.** Klobuchar / BDGIM ionosphere and UTC parameters are
        excluded: they are broadcast far less often, so apply them when present and treat
        `nothing` as zero rather than waiting. Group delay and inter-signal corrections
        (`T_GD`, `ISC_*`) split on one question: *does this signal's own correction ride
        in the block the gate already waits for?*

        If so, it is required, being a metre-scale bias that costs nothing to wait for.
        BeiDou B1C's `T_GD_B1Cp` and `ISC_B1Cd` are in the same CRC-protected subframe 2 as
        the ephemeris and clock, B2b's `T_GD_B2bI` in the same MT30 as the clock, and
        B1I/B3I's `T_GD1` in the same D1 subframe 1 / D2 page 1 as the clock.

        If not, it is excluded, because requiring it would mean waiting for a message the
        ICD does not schedule: BeiDou B2a's `T_GD_B2ap` / `ISC_B2ad` are in MT30 alone
        while the clock is in MT30-34, and the GPS CNAV `T_GD` and `ISC_*` arrive on their
        own cadence.

        Corrections for a band this decoder is not on (B1C's `T_GD_B2ap`, B1I/B3I's
        `T_GD2` since there is no B2I signal here, the GPS `ISC_L5Q5` on an L1 fix) are
        never required. They are still decoded and published; only this gate ignores them.

      - **Alert flag.** `is_sat_healthy` reflects the broadcast health bits only; check the
        L1 C/A alert flag (or equivalent) separately.
"""
is_decoding_completed_for_positioning(state::GNSSDecoderState) =
    is_decoding_completed_for_positioning(state.data)

# ---- Signal metadata --------------------------------------------------------
#
# A decoder state is 1:1 with the signal it demodulates, so everything GNSSSignals knows
# about that signal (name, constellation, band, message rate, time system of week numbers)
# follows from the state alone. To avoid a second source of truth, each signal file states
# only the one fact that is ours (the constants -> signal mapping, `get_signal_type`) and
# every accessor below forwards through it.

"""
    get_signal_type(state::GNSSDecoderState) -> Type{<:AbstractGNSSSignal}
    get_signal_type(constants::AbstractGNSSConstants) -> Type{<:AbstractGNSSSignal}

The GNSSSignals signal type whose navigation message a decoder demodulates, e.g.
`GPSL1CA` for a [`GPSL1CADecoderState`](@ref). The single mapping from this package's
decoders back into GNSSSignals: every accessor under "Signal metadata" in the docs forwards
through it, and it is the entry point for anything not forwarded (`get_signal_type(state)`
works wherever a GNSSSignals accessor takes a signal type).

Returns the *type*, not an instance: constructing a signal builds its spreading code matrix
and SIMD lookup table, which a symbol-domain decoder never needs. As a type it folds to a
compile-time constant, so forwarded accessors cost nothing at run time.

Stated per signal file and dispatched on the constants type rather than the data type,
since only the constants tell apart decoders sharing a data container: GPS L5-I and L2C-M
both decode into `GPSCNAVData` but are distinct signals (`GPSCNAVConstants{:GPSL5I}` vs
`GPSCNAVConstants{:GPSL2CM}`).

The mapping names the data-bearing component of a signal pair: `GPSL2CM` (not the `GPSL2CL`
pilot), `GalileoE5aI` (not `GalileoE5aQ`). Replica variants map to the same canonical type:
`GalileoE1B()` (default BOC(1,1)) and `GalileoE1B(CBOC)` both construct a decoder reporting
`GalileoE1B`, as the decoder does not retain which replica supplied its symbols. Forwarded
`get_signal_id` and `get_signal_name` therefore report `:GalileoE1B` and `"Galileo E1B"` for
both. Replica identity belongs to the tracking signal (GNSSSignals' `get_modulation`);
modulation is not forwarded for a decoder.

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
# `get_band_name`),
# and the time system the decoded week numbers and times of week are referenced
# to (`get_time_system` and friends — the Galileo epoch
# is 13 s
# shy of a UTC midnight *and* anchored to GPS Time rather than to UTC, so
# sourcing them rather than restating them matters).
#
# These extend GNSSSignals' own functions rather than introducing GNSSDecoder
# names, so a caller already doing `using GNSSSignals` needs no new import and
# can query a decoder state exactly as it queries a signal. Accessors that
# describe the *code* (`get_code_length`, `get_chip_rate`,
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
    :get_time_system,
    :get_time_system_id,
    :get_time_system_name,
)
    @eval @inline GNSSSignals.$accessor(state::GNSSDecoderState) =
        GNSSSignals.$accessor(get_signal_type(state))
end

# The symbol rate is the satellite's: a BeiDou B1I or B3I GEO satellite broadcasts its
# D2 message at 500 Hz, where the signal's rate is D1's 50 Hz, so the PRN goes along.
@inline GNSSSignals.get_bit_rate(state::GNSSDecoderState) =
    GNSSSignals.get_bit_rate(get_signal_type(state), state.prn)

# ---- Satellite and time metadata --------------------------------------------
#
# Three questions a positioning engine asks of every decoder: the satellite's orbit, the
# time of week its message stamps, and how its system time relates to another's. Each is
# answered from broadcast data by some signals and from constellation structure by others.
# Deriving a broadcast quantity is decoding work: it belongs here, next to the ICD
# knowledge, stated once per signal rather than once per consumer.

"""
    OrbitClass

Orbit class of a navigation satellite, reported by [`get_orbit_class`](@ref). The values
are deliberately *not* exported, like the health enums' values: `medium_earth_orbit` is a
name a receiver or orbit-propagation package may define, and consumers compare against
them rarely enough that `GNSSDecoder.geostationary_orbit` costs nothing.

# Values

  - `geostationary_orbit`: GEO
  - `inclined_geosynchronous_orbit`: IGSO
  - `medium_earth_orbit`: MEO

The integer values are this package's own, not any ICD's encoding; BeiDou's 2-bit
`sat_type` (1 GEO / 2 IGSO / 3 MEO) is mapped explicitly.
"""
@enum OrbitClass begin
    geostationary_orbit
    inclined_geosynchronous_orbit
    medium_earth_orbit
end

"""
    get_orbit_class(state::GNSSDecoderState) -> Union{Nothing,OrbitClass}

Orbit class of the satellite this decoder is tracking, or `nothing` when the signal does
not make it knowable.

  - **GPS and Galileo**: [`medium_earth_orbit`](@ref OrbitClass) for every satellite (both
    are MEO-only by design), so no decoded data is needed.
  - **BeiDou B1C, B2a, B2b**: from the broadcast `sat_type` field of the satellite's
    Ephemeris I block (1 GEO / 2 IGSO / 3 MEO, reserved 0; B1C and B2a Table 7-8, B2b
    Table 7-6). `nothing` until that block is decoded and promoted, and for the reserved
    code (the screen `is_known_sat_type` applies to positioning).
  - **BeiDou B1I and B3I**: [`geostationary_orbit`](@ref OrbitClass) for the GEO PRNs,
    `nothing` otherwise. D1/D2 NAV broadcasts no orbit-type field and the ICD's PRN
    partition (BDS-SIS-ICD-B1I-3.0 Table 4-1) separates only GEO from non-GEO, so IGSO
    and MEO cannot be told apart.

!!! warning "`nothing` is not `medium_earth_orbit`"

    `nothing` means "this signal cannot say" (on B1C, B2a and B2b also "not decoded
    yet"), not "not GEO" and not "MEO". Testing `=== GNSSDecoder.geostationary_orbit`
    is exact once a class is reported. Testing `=== GNSSDecoder.medium_earth_orbit` is
    *not*: a B1I MEO satellite reports `nothing`.
    Treat `nothing` as unknown and branch on the GEO case.

# Example

```julia
# Which reference does this BeiDou satellite's ephemeris use?
if get_orbit_class(state) === GNSSDecoder.geostationary_orbit
    # GEO: use the GEO/IGSO semi-major-axis reference
end
```

# See Also

  - [`OrbitClass`](@ref): the returned enumeration
"""
function get_orbit_class end

# GPS and Galileo are MEO-only, a constellation fact stated on the supertype. There is
# deliberately no `AbstractGNSSData` fallback: a constellation with GEO or IGSO members
# (QZSS, NavIC) must state its own answer rather than silently inherit "MEO".
get_orbit_class(::GNSSDecoderState{<:AbstractGPSData}) = medium_earth_orbit
get_orbit_class(::GNSSDecoderState{<:AbstractGalileoData}) = medium_earth_orbit

"""
$(TYPEDSIGNATURES)

Seconds of week that the decoder's validated navigation data stamps, or `nothing` before
a time has been decoded. It is the time the message broadcasts for the epoch it stamps, not
the current one (see below).

Every signal except Galileo E6-B broadcasts one, but GPS L1C-D counts two-hour intervals
plus 18-second TOI steps and BeiDou B1C counts hours plus 18-second SOH steps; the rest
carry a plain `TOW` or `SOW`. The ICD's reconstruction is done here rather than in every
consumer:

| Signal                                 | Broadcast fields | Seconds of week                                                                               |
|:-------------------------------------- |:---------------- |:--------------------------------------------------------------------------------------------- |
| GPS L1 C/A, GPS CNAV (L2C/L5)          | `TOW`            | `TOW`                                                                                         |
| GPS L1C-D                              | `ITOW`, `TOI`    | `ITOW·7200 + TOI·18`, `TOI = 0` carrying into the next interval (IS-GPS-800 §3.5.2, §3.5.3.2) |
| Galileo I/NAV (E1-B, E5b-I), E5a F/NAV | `TOW`            | `TOW`                                                                                         |
| BeiDou B1I, B3I (D1/D2)                | `SOW`            | `SOW`                                                                                         |
| BeiDou B1C                             | `HOW`, `SOH`     | `HOW·3600 + SOH·18` (BDS-SIS-ICD-B1C §7.3)                                                    |
| BeiDou B2a, B2b                        | `SOW`            | `SOW`                                                                                         |

Galileo E6-B always answers `nothing`: C/NAV carries HAS corrections stamped with a time
*of hour* per correction block, which stay on the blocks that own them (see
[`GalileoHASCorrectionBlock`](@ref)).

# The epoch it denotes, and reading the current time

Every signal's value is the time of the epoch that
[`num_bits_since_time_reference`](@ref GNSSDecoderState) is anchored to; the anchor is set
per signal at promotion time to the frame the broadcast field stamps (the *next* subframe
for GPS LNAV, CNAV and L1C-D; the *current* frame for BeiDou). So the current time of week
is the same expression on every signal:

```julia
using Unitful: Hz, ustrip
tow = get_broadcast_time_of_week(state)
now = tow + state.num_bits_since_time_reference / ustrip(Hz, get_bit_rate(state))
```

It reads [`data`](@ref GNSSDecoderState) rather than `raw_data`, so it returns a time only
once the message set has passed CRC/parity and the cross-subframe consistency check (the
gate [`is_decoding_completed_for_positioning`](@ref) reports). Call it on `state.raw_data`
for an unvalidated time.

# See Also

  - [`get_time_offset`](@ref): relating that time to another constellation's
"""
get_broadcast_time_of_week(state::GNSSDecoderState) = get_broadcast_time_of_week(state.data)

"""
$(TYPEDEF)

One broadcast time offset between two GNSS time systems, normalised across the five shapes
the ICDs state it in. [`get_time_offset`](@ref) returns it.

The offset at time `t` (seconds of week) in week `WN` of the *broadcasting* satellite's own
time system is

```julia
Δτ = t - t_0 + 604800 * (WN - WN_0)
Δt = A_0 + A_1 * Δτ + A_2 * Δτ^2
```

and `Δt` is **the broadcasting system's time minus the target system's**,
`t_own - t_target`, the direction every ICD here uses (Galileo OS SIS ICD §5.1.8;
IS-GPS-800/200 GGTO; BDS-SIS-ICD-B1C §7.13.2 Eq. 7-30). To convert a seconds-of-week
reading, *subtract*: `t_target = t_own - Δt`.

!!! note "`A_0` is not the broadcast coefficient alone"

    Two systems differ by a *defined* whole-second offset plus a broadcast steering
    residual; only the residual is on the air. The bias coefficient is 16 bits at 2⁻³⁵ s
    (±0.95 µs), so it cannot express a whole-second offset: real BeiDou satellites
    broadcast tens of nanoseconds against a BDT-to-GPST offset of 14 seconds.

    `A_0` is the sum of both, so the subtraction above holds as written. The defined part
    is `get_time_system_offset(target, get_time_system(state))`, each scale being
    a fixed offset from TAI with no leap seconds of its own; the broadcast part stays raw
    on the decoded data (`data.A_0BGTO`, `data.A_0G`, `data.A_0GGTO`).

    A Galileo decoder never shows it: GST and GPST are both TAI − 19 s, so the defined part
    is exactly zero. It is −14 s on every BeiDou signal (19 − 33: a BDT reading is 14 s
    *behind*, so converting away from BDT subtracts a negative), which is why the two
    cases must not be told apart by testing one. IS-GPS-800J §3.5.4.2.1.1 states the
    split: the GGTO parameters provide the offset "modulo one second", and "users must
    also apply any integer seconds difference between the systems using definitions of
    each system time scale"; the integer part is what this field folds in.

!!! note "`WN_0` is resolved, not as broadcast"

    The reference week is sent truncated, not always to the width of the week it is
    subtracted from (Galileo's `WN_0G` is 6 bits against a 12-bit `WN`). `WN_0` here is
    lifted into the decoder's week numbering, so `WN - WN_0` is a real week count on every
    signal. The raw field stays on the decoded data (`data.WN_0G`, `data.WN_0BGTO`,
    `data.WN_GGTO`).

!!! warning "Seconds of week, not week numbers"

    This converts a time *of week*, not which week. The scales do not share a week origin
    (BDT week 0 is GPS week 1356, exactly 9492 days later), so a full `(WN, TOW)`
    conversion also needs the epoch difference from `get_week_start_epoch` and must carry
    the rollover when adding `Δt` crosses a week boundary.

Signals that broadcast fewer terms are widened rather than given a separate type:
Galileo's GGTO has no quadratic term and reports `A_2 == 0`, and BeiDou D1/D2's two-term
offsets carry no reference epoch and report `t_0 === nothing` / `WN_0 === nothing`. Both
are exact: a zero coefficient contributes nothing, and an absent epoch is reported as
absent. It has a defined meaning: the D1/D2 polynomial's argument is the time of week
itself (BDS-SIS-ICD-B1I-3.0 §5.2.4.19, `Δt_GPS = A_0GPS + A_1GPS · t_E`), so
`t_0 === nothing` evaluates as `Δτ = t`, equivalent to `t_0 = 0` with `WN_0` the current
week, which is why such an offset needs no week number.

# Fields

$(TYPEDFIELDS)

# See Also

  - [`get_time_offset`](@ref): the accessor that produces one
"""
struct GNSSTimeOffset{T<:TimeSystem}
    """
    Time scale the offset is *to*; the satellite's own scale is `get_time_system(state)`
    """
    target::T
    """
    Bias (s): the broadcast steering residual *plus* the defined whole-second offset between
    the two scales (see above)
    """
    A_0::Float64
    """
    Drift (s/s)
    """
    A_1::Float64
    """
    Drift rate (s/s²); zero where only two terms are broadcast
    """
    A_2::Float64
    """
    Reference time of week (s); `nothing` on BeiDou D1/D2
    """
    t_0::Union{Nothing,Int}
    """
    Reference week, lifted from its truncated broadcast field into the decoder's week
    numbering (see `resolve_reference_week`), so `WN - WN_0` is a real week count;
    `nothing` on BeiDou D1/D2. Valid only in that difference: across a week-number rollover
    it can fall outside the field's range, even below zero, so it is no absolute week to
    place an epoch with
    """
    WN_0::Union{Nothing,Int}
end

"""
    get_time_offset(state::GNSSDecoderState, target::TimeSystem) -> Union{Nothing,GNSSTimeOffset}

The broadcast offset from this satellite's own time system to `target`, as a
[`GNSSTimeOffset`](@ref), or `nothing` when this decoder has no such offset.

`target` is a `GNSSSignals.TimeSystem`: `GPST()`, `GST()` or `BDT()`. The satellite's own
scale is `get_time_system(state)`; asking for it returns `nothing`, since the offset from a
scale to itself is not a broadcast quantity.

`nothing` also does not distinguish the three ways a signal can lack one:

  - it broadcasts no inter-system offset at all (GPS L1 C/A, Galileo E6-B),
  - it has not decoded one yet,
  - it has one, but for a different system than `target`: the GPS, B2a and B2b messages
    carry a single offset set tagged with its system (the satellite's choice), Galileo's
    is GPS-only by definition, and BeiDou B1C keys a set per system.

The ICDs' "not available" sentinels are screened, so a returned offset is real: a
GPS/BeiDou `GNSS_ID` of 0 means "no offset in this message" (IS-GPS-800 Table 6.2-7,
BDS-SIS-ICD-B1C §7.13.1), and Galileo marks an absent GGTO with all four fields all ones
(§5.1.8, handled at decode by `galileo_ggto`).

GLONASS offsets are decoded and published as fields but cannot be asked for here:
`GNSSSignals` defines no GLONASS `TimeSystem` and this package decodes no GLONASS signal.
Read `data.GNSS_ID` and the coefficient fields directly for those.

# Example

```julia
# Steer a BeiDou observable onto GPS time.
offset = get_time_offset(state, GPST())
if !isnothing(offset)
    # D1/D2 broadcasts no reference epoch (`t_0 === nothing`); its polynomial's
    # argument is the time of week itself (BDS-SIS-ICD-B1I-3.0 §5.2.4.19).
    Δτ = isnothing(offset.t_0) ? tow : tow - offset.t_0 + 604800 * (wn - offset.WN_0)
    t_gps = t_bdt - (offset.A_0 + offset.A_1 * Δτ + offset.A_2 * Δτ^2)
end
```

# See Also

  - [`GNSSTimeOffset`](@ref): the returned record, and the sign convention
  - [`get_broadcast_time_of_week`](@ref): the `t` the polynomial is evaluated at
"""
function get_time_offset end

# Every signal states its own method in its own file. No blanket `nothing` fallback, as for
# `get_orbit_class`: a future signal broadcasting an offset would silently inherit "none"
# and its correction would go unapplied.

"""
    broadcast_time_offset(state, target, A_0, A_1, A_2,
                          t_0, WN_0, WN, WN_0_modulus) -> Union{Nothing,GNSSTimeOffset}

Assemble a [`GNSSTimeOffset`](@ref) from one message's broadcast coefficients, adding the
defined whole-second offset between this decoder's time system and `target` to the
broadcast bias.

Every `get_time_offset` method routes through this instead of the constructor, because the
added term is the one a call site is liable to omit: it is exactly zero on the Galileo and
GPS paths (GST and GPST are both TAI − 19 s) and 14 s on every BeiDou one, so an omission
is invisible to any test that does not run a BeiDou decoder. Adding it where the record is
built gets it right in one place.

The scales are all fixed offsets from TAI without leap seconds of their own, so the term is
constant per pair and `get_time_system_offset` is its whole source; no leap-second table
is needed.

Positional arguments on purpose: `t_0`, `WN_0` and `WN` are `Union{Nothing,Int}` fields at
every call site, and as keywords they would form a non-concrete `NamedTuple`, heap-allocated
per call, on a path (`get_time_offset` is read once per satellite per epoch) that must stay
allocation-free.
"""
function broadcast_time_offset(
    state::GNSSDecoderState,
    target::TimeSystem,
    A_0::Real,
    A_1::Real,
    A_2::Real,
    t_0::Union{Nothing,Integer},
    WN_0::Union{Nothing,Integer},
    WN::Union{Nothing,Integer},
    WN_0_modulus::Int,
)
    reference_week = resolve_reference_week(WN_0, WN, WN_0_modulus)
    # An epoch that cannot be placed in the decoder's week numbering cannot be evaluated
    # against (and `Δτ` needs `WN` anyway): no usable offset yet.
    isnothing(reference_week) && !isnothing(WN_0) && return nothing
    defined = get_time_system_offset(target, get_time_system(state))
    GNSSTimeOffset(
        target,
        A_0 + defined,
        A_1,
        A_2,
        isnothing(t_0) ? nothing : Int(t_0),
        reference_week,
    )
end

"""
    resolve_reference_week(WN_0, WN, modulus) -> Union{Nothing,Int}

Lift a truncated broadcast reference week `WN_0` into the numbering of the decoder's own
week `WN`, picking the representative nearest to `WN`.

Every message sends the reference week truncated, not always to the width of the week
number it is subtracted from:

| field             | width   | against a `WN` of |
|:----------------- |:------- |:----------------- |
| Galileo `WN_0G`   | 6 bits  | 12 bits           |
| GPS `WN_GGTO`     | 13 bits | 13 bits           |
| BeiDou `WN_0BGTO` | 13 bits | 13 bits           |

So `WN - WN_0` as broadcast is meaningful on GPS and BeiDou but not on Galileo, where it
subtracts a 0..63 count from a 0..4095 one. Unresolved it fails silently: `Δτ` is inflated
by up to 4032 weeks, which a full-scale `A_1G` turns into 0.7 ms (200 km of range) on a
quantity whose true value is tens of nanoseconds.

The *nearest* representative, not the most recent, matters at a week boundary: a reference
epoch a few hours ahead of the current week is legal, and rounding down would place it
`modulus - 1` weeks in the past. On the two 13-bit fields the result is `WN_0` unchanged
except across a week-number rollover, where it keeps the difference small.

`nothing` when either week is unknown.
"""
function resolve_reference_week(
    WN_0::Union{Nothing,Integer},
    WN::Union{Nothing,Integer},
    modulus::Int,
)
    (isnothing(WN_0) || isnothing(WN)) && return nothing
    half = modulus ÷ 2
    return Int(WN) - (mod(Int(WN) - Int(WN_0) + half, modulus) - half)
end

"""
Soft-symbol buffer accessor; the per-signal cache stores it as `soft_buffer`.
"""
soft_buffer(state::GNSSDecoderState) = state.cache.soft_buffer

"""
Number of soft symbols currently buffered.
"""
num_symbols_buffered(state::GNSSDecoderState) = length(soft_buffer(state))

"""
    push_soft_symbol!(state, sym)

Push one `Float32` soft symbol onto the per-signal circular deque (capacity
`sync_unit_length + preamble_length`); once full, the oldest sample is dropped via
`popfirst!`.
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

Positive soft symbol ⇒ bit 0, negative ⇒ bit 1 (AFF3CT's LLR convention). Returns `true`
for bit 1.
"""
@inline hard_slice(sym::Real) = sym < zero(sym)

"""
    pack_soft_buffer(T, soft_buffer, total_bits)

Hard-slice the leading `total_bits` of `soft_buffer` (oldest first) into a `T<:Unsigned`
packed-bit buffer, MSB = oldest bit.
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

Hard-slice `len` soft symbols of `deque` starting at 1-based `offset` into a `UInt64`,
oldest symbol at the MSB (the bit order of [`pack_buffer`](@ref), for a slice).

Sync inspects only a few known positions (the preamble at either end, any unencoded header
field), so long, high-rate windows read just those bits instead of repacking the window on
every symbol.

`len` must be at most 64 and `offset + len - 1` within the deque: the read is `@inbounds`,
so the caller owns both bounds. Every caller is a `try_sync` hook, which the decode loop
runs only behind `is_enough_buffered_symbols_to_decode`.
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

Hard-slice `len` soft symbols of `deque` starting at 1-based `offset` into a `UInt64`
with **bit 0 = the first symbol**: the opposite end from [`pack_bits_msb_first`](@ref) and
the order of the codeword tables (`BCH_TOI_CODEWORDS`, `b1c_prn_codeword`,
`b1c_soh_codeword`; cf. [`pack_hard_codeword`](@ref), the same over an iterable of hard
decisions).

The BCH-synchronised signals (GPS L1C-D on its TOI codeword, BeiDou B1C on its PRN + SOH
pair) compare against those tables once per symbol at both ends of the window, so this
reads the deque in place instead of materialising a slice.

`len` must be at most 64 and `offset + len - 1` within the deque: the read is `@inbounds`,
so the caller owns both bounds (the length check of `pack_hard_codeword` is not available
here). Every caller is a `try_sync` hook, run only behind
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

Hard-slice the leading `sync_unit_length + preamble_length` soft symbols of the soft
buffer into a packed-bit value (oldest bit at MSB). The unsigned type is signal-specific,
from `packed_buffer_type`. The result is a plain value threaded through the sync path, not
stashed in mutable cache state.
"""
function pack_buffer(state::GNSSDecoderState)
    n = state.constants.sync_unit_length + state.constants.preamble_length
    pack_soft_buffer(packed_buffer_type(state), soft_buffer(state), n)
end

"""
    try_sync(state) -> Union{Nothing,Unsigned}

Default per-signal sync hook: on a preamble match (visible at both ends of the candidate
sync unit, in either polarity; see `detect_preamble_pair`), hard-slice the deque into a
packed-bit buffer via `pack_buffer` and return it, else `nothing`. Returning the buffer
lets the caller reuse it and keeps it out of mutable cache state.

Signals override it where they sync differently (e.g. GPS L1C-D's TOI BCH match).
"""
function try_sync(state::GNSSDecoderState)
    # Match the preamble at either end first and pack the whole window only on a match:
    # packing hundreds of symbols per incoming symbol is otherwise most of the cost.
    constants = state.constants
    isnothing(
        detect_preamble_pair(
            soft_buffer(state),
            constants.preamble,
            constants.preamble_length,
            constants.sync_unit_length,
        ),
    ) && return nothing
    return pack_buffer(state)
end

"""
    detect_preamble_pair(deque, preamble, preamble_length, sync_unit_length)
        -> Union{Nothing,Bool}

Match the preamble at both ends of the candidate window (the start of this sync unit and
of the next) by slicing the soft-symbol deque directly, instead of hard-slicing the whole
window into a packed integer first.

Returns `nothing` when there is no sync, otherwise the resolved polarity: `true` when both
ends carry the *inverted* preamble. Both ends must match in a common polarity, which
resolves the 180-degree carrier ambiguity.

The form every FEC-bearing decoder wants: the default `pack_buffer` path shifts a
window-wide integer per symbol to expose a preamble at two known offsets, and these
decoders consume their payload as soft symbols anyway. Hence no `packed_buffer_type` method
for Galileo I/NAV, E5a or E6-B, nor BeiDou B2b. Measured per `try_sync` call on a full
window: 195 ns → 18 ns for I/NAV's 260-symbol window and 521 ns → 18 ns for E5a's
512-symbol one (0.26 → 0.09 and 0.57 → 0.07 µs per symbol end to end). The two 1000 sps
decoders avoid a 1016-bit shift per symbol.

The head is tested before the tail is packed: on noise the head matches with probability
`2^-(preamble_length-1)`, so packing the tail first is almost always wasted.

!!! note "Pass the lengths, not the state"

    This runs once per *symbol*, the hottest path in the package. The lengths are
    arguments, not read from `state.constants`, so callers can pass module-level constants
    that fold at compile time and let `pack_bits_msb_first` unroll; struct fields do not.
    At 1000 sps, reading them from `constants` costs about 14 % of the per-symbol budget.
"""
@inline function detect_preamble_pair(
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
    copy_soft_window!(dst, deque, offset, num_symbols, is_polarity_inverted) -> dst

Copy `num_symbols` soft symbols out of `deque`, starting at 1-based `offset + 1`,
resolving the 180-degree polarity on the way: an inverted stream is a negated one, so the
LLR magnitudes (the FEC confidences) are unchanged.

`dst` is the caller's long-lived scratch buffer from its cache; the copy does not allocate.
Every FEC-bearing decoder needs this step between "sync found" and the FEC decode.
"""
function copy_soft_window!(
    dst::AbstractVector{Float32},
    deque::CircularDeque{Float32},
    offset::Int,
    num_symbols::Int,
    is_polarity_inverted::Bool,
)
    sign = is_polarity_inverted ? -1.0f0 : 1.0f0
    @inbounds for i = 1:num_symbols
        dst[i] = sign * deque[offset+i]
    end
    return dst
end

"""
    complement_buffer_if_necessary(state, buffer) -> (state, resolved_buffer)

If the newest preamble in `buffer` is the *inverted* preamble, return the state flagged
`is_polarity_inverted = true` with the complemented buffer; otherwise the state
flagged `false` with `buffer` unchanged. The polarity-resolved buffer is for the
per-signal `decode_sync_unit` to consume.
"""
function complement_buffer_if_necessary(state::GNSSDecoderState, buffer)
    mask = calc_preamble_mask(state.constants)
    if buffer & mask == ~state.constants.preamble & mask
        return GNSSDecoderState(state; is_polarity_inverted = true), ~buffer
    else
        return GNSSDecoderState(state; is_polarity_inverted = false), buffer
    end
end

"""
    drain_after_sync!(state)

Drop the consumed `sync_unit_length` oldest soft symbols from the deque, keeping the
trailing `preamble_length` symbols as the leading preamble of the next subframe.

Drops at most `length(deque)` symbols: a `decode_sync_unit` hook may reset the decoder
mid-frame (e.g. GPS L1C-D on a TOI discontinuity), emptying the buffer, and the
unconditional drain would then `popfirst!` an empty `CircularDeque` and throw.
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
holds, and return the updated state. Detects preambles in the incoming soft bits and
decodes sync units into navigation data, handling normal and 180-degree phase-shifted
signals automatically. Allocates nothing.

`decode!` writes straight into the buffers `state` was constructed with: the soft-symbol
buffer, the FEC scratch, the voting tallies in `state.cache`, and every container
referenced from `state.raw_data` and `state.data` (almanac stores, health tables, HAS
masks, ...). All are sized once at construction:

```julia
state = GPSL1CADecoderState(25)
for chunk in soft_symbol_chunks
    state = decode!(state, chunk, length(chunk))
end
```

!!! warning "Earlier states are overwritten"

    The returned state shares its storage with `state`, and later calls keep overwriting
    it. Treat the return value as *the* live decoder and do not use `state` (or anything
    read out of its containers, such as `state.data.almanacs`) afterwards. Take a
    [`copy`](@ref Base.copy(::GNSSDecoderState)) first for a snapshot.

# Soft-bit convention

`soft_bits` are what a tracking loop hands over, such as TrackingLoops' `get_soft_bits`:
one per channel symbol as tracked, before any error-correction decoding (the bits
GNSSSignals' `get_bit_rate` counts). With forward error correction they are not the
message's bits; the decoded message bits exist only inside this package.

`soft_bits` is an `AbstractVector{<:Real}`; `Float32` is canonical. The sign carries the
bit decision and the magnitude the confidence (standard LLR convention):

  - **positive ⇒ bit 0**, **negative ⇒ bit 1**, but only as a *convention*. The absolute
    polarity of a Costas-tracked signal is 180°-ambiguous, so the decoder matches the
    preamble in either polarity and flips internally (recording it in
    `is_polarity_inverted`). The opposite sign decodes the same data; only the
    polarity flag differs. A tracking loop's soft bits may come in either polarity, which
    is harmless for this reason.
  - magnitude ⇒ confidence. **No normalization is required; values need not lie in
    `[-1, 1]`.** GPS L1 C/A (hard-slice + parity) and Galileo E1B (Viterbi, whose ML path
    is invariant to a global scale) use only the sign. The LDPC decodes (GPS L1C-D and the
    BeiDou B-CNAV family: B1C, B2a, B2b) are flooding sum-product, which *is*
    scale-sensitive: give them confidence-weighted magnitudes on a roughly LLR-like scale
    (`≈ 2·r/σ²`) for best performance at marginal SNR, but no fixed range is needed.

# Arguments

  - `state::GNSSDecoderState`: Current decoder state
  - `soft_bits::AbstractVector{<:Real}`: Soft bits to consume, oldest first
  - `num_bits::Int`: Number of leading entries of `soft_bits` to process

# Keywords

  - `decode_once::Bool=false`: If `true`, stops once all required positioning
    data has been validated (subframes 1-3 for GPS L1 C/A; word types 1-5 for
    Galileo E1B)

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
    num_bits::Int;
    decode_once::Bool = false,
)
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
        state = with_num_bits_since_time_reference(state, num_bits_since_time_reference)
        state, resolved_buffer = complement_buffer_if_necessary(state, buffer)
        state = decode_sync_unit(state, resolved_buffer)
        if !decode_once || !is_decoding_completed_for_positioning(state.data)
            state = validate_data(state)
        end
        state = drain_after_sync!(state)
        num_bits_since_time_reference = state.num_bits_since_time_reference
    end
    return with_num_bits_since_time_reference(state, num_bits_since_time_reference)
end

"""
    push_until_sync!(state, soft_symbols, first, last, num_bits_since_time_reference)
        -> (next, buffer, num_bits_since_time_reference)

The per-symbol half of [`decode!`](@ref): push `soft_symbols[first:last]` onto `state`'s
soft-symbol buffer one at a time, counting them in `num_bits_since_time_reference` (the
running counter, or `nothing` before the first sync), until `try_sync` finds a sync.
Returns the index of the next unconsumed symbol, the sync buffer (`nothing` if the symbols
ran out) and the updated count.

It never rebuilds `state`, which keeps it cheap: `state` holds the decoded data inline, and
a loop reassigning it copies it every iteration. No `try_sync` hook reads the counter, so
it is written back only once a sync is to be decoded.
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

# Split on the counter being set so each rebuild has a concretely typed keyword (a `Union`
# keyword value takes the allocating path on Julia 1.10).
with_num_bits_since_time_reference(state::GNSSDecoderState, num_bits_since_time_reference) =
    num_bits_since_time_reference === nothing ?
    GNSSDecoderState(state; num_bits_since_time_reference = nothing) :
    GNSSDecoderState(state; num_bits_since_time_reference)

"""
    copy(state::GNSSDecoderState) -> GNSSDecoderState

Independent copy of a decoder state: soft-symbol buffer, cache and every container
referenced from `raw_data` and `data` are copied (see `duplicate`), so [`decode!`](@ref) on
the copy overwrites nothing `state` can see, and vice versa. Buffer capacities are kept, so
the copy decodes without allocating too.
"""
Base.copy(state::GNSSDecoderState) = GNSSDecoderState(
    state.prn,
    duplicate(state.raw_data),
    duplicate(state.data),
    state.constants,
    duplicate(state.cache),
    state.num_bits_since_time_reference,
    state.is_polarity_inverted,
)

"""
    reset_decoder_state!(state::GNSSDecoderState) -> GNSSDecoderState

Reset a decoder after a signal loss or reacquisition, **overwriting** the buffers `state`
holds (the soft-symbol buffer is emptied), and return the reset state. Allocates nothing.
What is reset (the time of week, the validated `data`) and what is kept for fast recovery
is per signal, documented on each signal's method. As with [`decode!`](@ref), use the
return value as the live decoder and not `state`.
"""
function reset_decoder_state! end

# ---- Shared decoder primitives ----------------------------------------------
#
# Signal-agnostic primitives used by more than one signal decoder.

# Packed-word integer types used by more than one decoder (`BitIntegers` widths are
# global, so each is defined once, here). `@define_integers` also defines the signed
# `Int288`, `Int320`, `Int512` and `Int600`.
#
#   - `UInt288`: the 220 bits a Galileo I/NAV page pair runs the CRC-24Q over, the 274-bit
#     GPS L1C-D subframe 3 and the 264-bit BeiDou B1C subframe 3.
#   - `UInt320`: a GPS L1 C/A subframe, a GPS L1C-D subframe, a GPS CNAV message (L5I /
#     L2C) and a BeiDou D1/D2 subframe: 300 data bits plus up to 8 trailing sync bits.
#   - `UInt512`: the 486-bit information word of a Galileo E6-B C/NAV page and of a BeiDou
#     B2b message.
#   - `UInt600`: the 600-bit subframe 2 of GPS L1C-D and of BeiDou B1C.
BitIntegers.@define_integers 288
BitIntegers.@define_integers 320
BitIntegers.@define_integers 512
BitIntegers.@define_integers 600

"""
    increment_voting(old_vote, max_vote) -> Int

Raise a broadcast-repetition vote by one, saturating at `max_vote`. Used by the two
decoders that promote `raw_data` to `data` by repetition voting instead of an
issue-of-data match: GPS L1 C/A and BeiDou D1/D2 (no IODs at all).
"""
increment_voting(old_vote, max_vote) = min(max_vote, old_vote + 1)

"""
    elapsed_symbols(anchor, now) -> Union{Nothing,Int64}

Symbols counted between the time-reference counter readings `anchor` and `now` (both
`num_bits_since_time_reference` values), or `nothing` if either is missing.
"""
elapsed_symbols(anchor, now) = isnothing(anchor) || isnothing(now) ? nothing : now - anchor
