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
(IS-GPS-200N §30.3.3.1.1, IS-GPS-800J §3.5.3), which [`get_broadcast_ephemeris`](@ref)
applies through this type.
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

# `has_positioning_messages(data)`, one method per message container, checks the required
# message set alone; publishing validated data and its accessors gate on that, not on the
# time anchor.
"""
$(TYPEDSIGNATURES)

Whether a decoder has recovered the minimum navigation data a positioning engine needs
from this satellite: a time of week, a full ephemeris set and the SV clock-correction
polynomial (plus, on signals that carry them in the same required set, the broadcast week
number and single-band group delay). Reads the validated [`data`](@ref GNSSDecoderState),
so it becomes `true` only once the required message set has passed CRC/parity **and** the
cross-subframe issue-of-data check that promotes `raw_data` to `data`. It also requires the
time anchor [`num_bits_since_time_reference`](@ref GNSSDecoderState), so whenever it is
`true`, [`get_decoded_time_of_week`](@ref) returns a number. Always `false` on Galileo
E6-B, whose HAS corrections apply to another signal's data and form no positioning set.

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
        while the clock is in MT30-34, and the GPS CNAV `T_GD` and `ISC_*` (message
        type 30) and the CNAV-2 ISCs of the other civil signals (subframe 3 page 1) arrive
        on their own cadence. CNAV-2's own `T_GD`, `ISC_L1CD` and `ISC_L1CP` ride in
        subframe 2 with the ephemeris.

        Corrections for a band this decoder is not on (B1C's `T_GD_B2ap`, B1I/B3I's
        `T_GD2` since there is no B2I signal here, the GPS `ISC_L5Q5` on an L1 fix) are
        never required. They are still decoded and published; only this gate ignores them.

      - **Alert flag.** `is_sat_healthy` reflects the broadcast health bits only; check the
        L1 C/A alert flag (or equivalent) separately.
"""
has_positioning_data(state::GNSSDecoderState) =
    has_positioning_messages(state.data) && !isnothing(state.num_bits_since_time_reference)

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
# to (`get_time_system` and friends: that system's epoch on the continuous TAI scale,
# `get_week_start_epoch(state, 0)`; the Galileo
# epoch is 13 s shy of a UTC midnight *and* anchored to GPS Time rather than to
# UTC, so sourcing them rather than restating them matters).
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

# The TAI epoch a week of the decoder's time system begins at, forwarded like the
# accessors above, so that a decoded week and time of week date an instant:
# `get_week_start_epoch(state, week) + time_of_week`.
@inline GNSSSignals.get_week_start_epoch(state::GNSSDecoderState, week::Integer) =
    GNSSSignals.get_week_start_epoch(get_signal_type(state), week)

# ---- Satellite and time metadata --------------------------------------------
#
# Questions a positioning engine asks of every decoder: the satellite's orbit, the week and
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

The value is the time of the epoch that
[`num_bits_since_time_reference`](@ref GNSSDecoderState) is anchored to; the anchor is set
per signal at promotion time to the frame the broadcast field stamps (the *next* subframe
for GPS LNAV, CNAV and L1C-D; the *current* frame for BeiDou). The current time of week is
available through [`get_decoded_time_of_week`](@ref):

```julia
time_of_week = get_decoded_time_of_week(state)
```

It reads [`data`](@ref GNSSDecoderState) rather than `raw_data`, so it returns a time only
once the message set has passed CRC/parity and the cross-subframe consistency check (the
data gate of [`has_positioning_data`](@ref)). Call it on `state.raw_data` for an unvalidated
time.

# See Also

  - [`get_broadcast_inter_system_time_offset`](@ref): relating that time to another
    constellation's
"""
get_broadcast_time_of_week(state::GNSSDecoderState) = get_broadcast_time_of_week(state.data)

"""
    get_decoded_time_of_week(state::GNSSDecoderState) -> Union{Nothing,Float64}

Time at the last channel symbol consumed by [`decode!`](@ref), in seconds relative to the
week of [`get_broadcast_time_of_week`](@ref): the validated broadcast stamp advanced by the
channel-symbol count at the decoder's symbol rate. `nothing` when the stamp or its time
anchor is unavailable, including on Galileo E6-B.

Deliberately not reduced modulo one week: it can exceed 604800 when symbols arrive past a
week boundary before a new stamp is validated. Keep it paired with the broadcast week;
wrapping it alone would lose the week carry. Like the broadcast week, it may inherit a
temporarily lagging week field.

This is the uncorrected satellite time. A ranging consumer adds the residual code/carrier
phase and applies satellite clock and group-delay corrections. Vector tracking may adjust
the time-reference counter to align it with its code phase; the counter does not report the
age of the last preamble.
"""
function get_decoded_time_of_week(state::GNSSDecoderState)
    stamp = get_broadcast_time_of_week(state)
    count = state.num_bits_since_time_reference
    (isnothing(stamp) || isnothing(count)) && return nothing
    Float64(stamp) + count / ustrip(Hz, get_bit_rate(state))
end

"""
    get_broadcast_week(state::GNSSDecoderState; approximate_year::Integer) -> Union{Nothing,Int}
    get_broadcast_week(data::AbstractGNSSData; approximate_year::Integer) -> Union{Nothing,Int}

Absolute week number of the decoder's validated navigation data, counted from the origin
of its own time system, or `nothing` before a week number has been decoded.

It is the week the message broadcasts, not the current one: it can lag across a week
crossover, by a subframe, until the message carrying `WN` is decoded again (see the
signals' `reset_decoder_state!` notes). GNSSSignals' `get_week_start_epoch` answers for a
decoder state, so the broadcast week and time of week date an instant on TAI:

```julia
week = get_broadcast_week(state; approximate_year = 2026)
time_of_week = get_broadcast_time_of_week(state)
if !isnothing(week) && !isnothing(time_of_week)
    epoch = get_week_start_epoch(state, week) + time_of_week  # a GNSSSignals.TAITime
end
```

Every message here broadcasts its week number wide enough to be absolute for decades
(Galileo I/NAV and F/NAV 12 bits; GPS CNAV, CNAV-2, BeiDou D1/D2 and B-CNAV 13 bits),
except GPS L1 C/A: LNAV sends only the week modulo 1024 (IS-GPS-200 §20.3.3.3.1.1), which
rolls over every 19.6 years (1999-08-22, 2019-04-07, 2038-11-21, ...). That rollover is
resolved against `approximate_year`: the 1024-week cycle that puts the broadcast week
nearest the middle of that year. Any year within ±9 years of the observation selects the
right cycle (the current year for a live signal, the recording year for a replay).

`approximate_year` is required on every signal, though only GPS L1 C/A reads it, so a call
reads the same on any decoder. There is no default: the decoder has no notion of the
current time.

Galileo E6-B always answers `nothing`: C/NAV broadcasts no week number (see
[`get_broadcast_time_of_week`](@ref)).

The state form reads [`data`](@ref GNSSDecoderState) rather than `raw_data`, so it answers
only once the message set has been validated; call it on `state.raw_data` for an
unvalidated week.

# See Also

  - [`get_broadcast_time_of_week`](@ref):
    the time within that week
"""
get_broadcast_week(state::GNSSDecoderState; approximate_year::Integer) =
    get_broadcast_week(state.data; approximate_year)

# Every signal states its own `get_broadcast_week(data; approximate_year)` beside its
# `get_broadcast_time_of_week`. There is no fallback reading `data.WN`: a future signal with
# a truncated week number must resolve it, not silently inherit the truncated value.

"""
    get_broadcast_inter_system_time_offset(state::GNSSDecoderState, target::TimeSystem,
        time_of_week) -> Union{Nothing,Float64}

The broadcast offset from this satellite's own time system to `target` (s), at
`time_of_week` (s) on the satellite's own scale, or `nothing` when this decoder has no such
offset.

`target` is a `GNSSSignals.TimeSystem`: `GPST()`, `GST()` or `BDT()`. The satellite's own
scale is `get_time_system(state)`; asking for it returns `nothing`, since the offset from a
scale to itself is not a broadcast quantity.

The offset `Δt` is **the broadcasting system's time minus the target system's**,
`t_own - t_target`, the direction every ICD here uses (Galileo OS SIS ICD §5.1.8;
IS-GPS-800/200 GGTO; BDS-SIS-ICD-B1C §7.13.2 Eq. 7-30). To convert a seconds-of-week
reading, *subtract* it and add the leap seconds between the scales, GNSSSignals'
`get_time_system_offset(get_time_system(state), target)` (see the note below):
`t_target = t_own - Δt + leap_seconds`.

It is the broadcast polynomial

```julia
Δτ = time_of_week - t_0 + 604800 * (WN - WN_0)
Δt = A_0 + A_1 * Δτ + A_2 * Δτ^2
```

evaluated against the decoder's own week number `WN`. The reference week `WN_0` is sent
truncated, not always to the width of `WN` (Galileo's `WN_0G` is 6 bits against a 12-bit
`WN`), so it is lifted into the decoder's week numbering first and the caller needs only the
time of week. `WN` is the decoded one, so right after a week crossover, until the message
carrying it is decoded again, it can still be the previous week's. Galileo's GGTO has no
quadratic term. BeiDou D1/D2's two-term offsets carry no reference epoch: their argument is
the time of week itself (BDS-SIS-ICD-B1I-3.0 §5.2.4.19, `Δt_GPS = A_0GPS + A_1GPS · t_E`).

!!! note "The steering residual alone"

    Two systems differ by the leap seconds between them plus a broadcast steering residual;
    only the residual is on the air. The bias coefficient is 16 bits at 2⁻³⁵ s (±0.95 µs),
    or on BeiDou D1/D2 14 bits at 0.1 ns (±0.82 µs), so it cannot express a whole-second
    offset: real BeiDou satellites broadcast tens of nanoseconds against a BDT-to-GPST
    offset of 14 seconds. IS-GPS-800J §3.5.4.2.1.1 states the split: the GGTO parameters
    provide the offset "modulo one second", and "users must also apply any integer seconds
    difference between the systems using definitions of each system time scale".

    `Δt` is evaluated from the coefficients as broadcast (as RINEX navigation files and
    other receivers report them), so it is the steering residual alone. Each of these
    scales is a fixed offset from TAI, set to TAI − UTC when defined and never adjusted
    (19 s for GPST, which GST adopted, and 33 s for BDT), so two of them differ by the leap
    seconds UTC gained in between: zero between GPST and GST, +14 s from BDT toward either
    (a BDT reading is 14 s *behind*). A conversion only ever tried on Galileo never notices
    them missing.

!!! warning "Seconds of week, not week numbers"

    This converts a time *of week*, not which week. The scales do not share a week origin
    (BDT week 0 begins 14 s into GPS week 1356), so a full `(WN, TOW)` conversion goes
    through TAI, which carries the leap seconds, epoch offset and week over in one step:

    ```julia
    epoch = get_week_start_epoch(BDT(), WN_bdt) + t_bdt - Δt  # a GNSSSignals.TAITime
    WN_gps, t_gps = get_week(GPST(), epoch), get_time_of_week(GPST(), epoch)
    ```

    GNSSSignals takes the whole weeks off in integer seconds, so the nanoseconds of `Δt`
    survive the round trip.

`nothing` does not distinguish the three ways a signal can lack an offset:

  - it broadcasts no inter-system offset at all (GPS L1 C/A, Galileo E6-B),
  - it has not decoded one yet,
  - it has one, but for a different system than `target`: the GPS, B2a and B2b messages
    carry a single offset set tagged with its system (the satellite's choice), Galileo's
    is GPS-only by definition, and BeiDou B1C keys a set per system.

The ICDs' "not available" sentinels are screened, so a returned offset is real: a
GPS/BeiDou `GNSS_ID` of 0 means "no offset in this message" (IS-GPS-800 Table 6.2-7,
BDS-SIS-ICD-B1C §7.13.1), and Galileo marks an absent GGTO with all four fields all ones
(§5.1.8, handled at decode by `galileo_ggto`).

[`get_inter_system_time_offset_parameters`](@ref) returns the polynomial itself, with its
reference week resolved, for a RINEX writer. The coefficients as broadcast stay on the
decoded data (`data.A_0G`, `data.A_0BGTO`, `data.A_0GGTO`, reference weeks `data.WN_0G`,
`data.WN_0BGTO`, `data.WN_GGTO`, ...). GLONASS offsets are decoded and published there too
but cannot be asked for here: `GNSSSignals` defines no GLONASS `TimeSystem`.

# Example

```julia
# Steer a BeiDou observable onto GPS time.
Δt = get_broadcast_inter_system_time_offset(state, GPST(), t_bdt)
if !isnothing(Δt)
    # The leap seconds between the scales are not broadcast: +14 s from BDT to GPST.
    t_gps = t_bdt - Δt + get_time_system_offset(BDT(), GPST())
end
```

# See Also

  - [`get_broadcast_time_of_week`](@ref):
    the time of week the offset is evaluated at
  - [`get_broadcast_inter_system_time_offset_epoch`](@ref): the reference epoch it is
    evaluated from, to judge its age by
  - [`get_inter_system_time_offset_parameters`](@ref): the broadcast polynomial itself
"""
function get_broadcast_inter_system_time_offset(
    state::GNSSDecoderState,
    target::TimeSystem,
    time_of_week::Real,
)
    offset = inter_system_time_offset_parameters(state, target)
    isnothing(offset) && return nothing
    (; A_0, A_1, A_2, t_0, WN_0) = offset
    # A resolved `WN_0` implies a decoded `WN`: the record is only built with one.
    Δτ =
        isnothing(t_0) || isnothing(WN_0) ? Float64(time_of_week) :
        time_of_week - t_0 + SECONDS_PER_WEEK * (something(state.data.WN) - WN_0)
    return A_0 + A_1 * Δτ + A_2 * Δτ^2
end

"""
    get_broadcast_inter_system_time_offset_epoch(state::GNSSDecoderState,
        target::TimeSystem) -> Union{Nothing,TAITime}

The reference epoch of the broadcast offset to `target` that
[`get_broadcast_inter_system_time_offset`](@ref) evaluates: its reference week and time of
week, dated on TAI as GNSSSignals' `get_week_start_epoch` dates the satellite's own time
system. `nothing` when there is no such offset or it carries no reference epoch.

The polynomial is extrapolated from this epoch, and a decoder keeps an offset set until the
satellite broadcasts a newer one (BeiDou B1C keeps one per system, which a satellite may
stop refreshing). How old is too old is the caller's decision, and where several
satellites broadcast the same offset, the latest epoch names the newest set:

```julia
epoch = get_week_start_epoch(BDT(), WN_bdt) + t_bdt
reference = get_broadcast_inter_system_time_offset_epoch(state, GPST())
age = isnothing(reference) ? nothing : epoch - reference  # s
```

The same epoch as [`get_inter_system_time_offset_parameters`](@ref)' `(WN_0, t_0)`, but a
`TAITime` rather than the record, so a navigation solver reading it per satellite per epoch
stays allocation-free on Julia 1.10 too. BeiDou D1/D2 answers `nothing` even where it has an
offset: its two-term offsets carry no reference epoch.
"""
function get_broadcast_inter_system_time_offset_epoch(
    state::GNSSDecoderState,
    target::TimeSystem,
)
    offset = inter_system_time_offset_parameters(state, target)
    isnothing(offset) && return nothing
    (; t_0, WN_0) = offset
    (isnothing(t_0) || isnothing(WN_0)) && return nothing
    return GNSSSignals.get_week_start_epoch(get_time_system(state), WN_0) + t_0
end

"""
$(TYPEDEF)

One broadcast offset from the satellite's own time system to `target`: the polynomial
`A_0 + A_1·Δτ + A_2·Δτ²` (s) at `Δτ` seconds from its reference epoch `(WN_0, t_0)`, the
form [`get_broadcast_inter_system_time_offset`](@ref) evaluates.
[`get_inter_system_time_offset_parameters`](@ref) returns it.

It normalises the five forms the ICDs broadcast: Galileo's GGTO toward GPS time (OS SIS ICD
§5.1.8), GPS CNAV and CNAV-2's GGTO toward the system their `GNSS_ID` names
(IS-GPS-705/800), BeiDou B-CNAV's BGTO likewise (BDS-SIS-ICD-B1C-1.0 §7.13), and BeiDou
D1/D2's offsets toward GPS and Galileo time (BDS-SIS-ICD-B1I-3.0 §5.2.4.19). The
coefficients are the steering residual as broadcast, without the leap seconds between the
scales. `A_2` is `0.0` where only two terms are broadcast (Galileo, D1/D2).

`WN_0` is a full week number on the decoder's own count, not the truncated field the
message carries. BeiDou D1/D2 broadcasts no reference epoch: `t_0` and `WN_0` are
`nothing` and the polynomial's argument is the BDT time of week itself.

[`get_broadcast_inter_system_time_offset_epoch`](@ref) dates the reference epoch on TAI, to
judge the offset's age by.

# Fields

$(TYPEDFIELDS)
"""
struct InterSystemTimeOffsetParameters{T<:TimeSystem}
    """
    Time scale the offset is *to*; the satellite's own scale is `get_time_system(state)`
    """
    target::T
    """
    Bias (s)
    """
    A_0::Float64
    """
    Drift (s/s)
    """
    A_1::Float64
    """
    Drift rate (s/s²); `0.0` where not broadcast
    """
    A_2::Float64
    """
    Reference time of week (s); `nothing` on BeiDou D1/D2
    """
    t_0::Union{Nothing,Int}
    """
    Reference week, on the decoder's own week count; `nothing` on BeiDou D1/D2
    """
    WN_0::Union{Nothing,Int}
end

# The broadcast offset from the satellite's own scale to `target` as a
# `InterSystemTimeOffsetParameters`, or `nothing` in the cases
# `get_broadcast_inter_system_time_offset` lists. Every signal states its own method in its
# own file. No blanket `nothing` fallback, as for `get_orbit_class`: a future signal
# broadcasting an offset would silently inherit "none" and its correction would go
# unapplied.
function inter_system_time_offset_parameters end

"""
    assemble_inter_system_time_offset_parameters(target, A_0, A_1, A_2, t_0, WN_0, WN,
                                                 WN_0_modulus)
        -> Union{Nothing,InterSystemTimeOffsetParameters}

Assemble an `InterSystemTimeOffsetParameters` from one message's broadcast coefficients,
resolving the truncated reference week `WN_0` against the decoder's week `WN` (see
`resolve_reference_week`). Every `inter_system_time_offset_parameters` method routes through
this instead of the constructor, so the resolution is right in one place.

Positional arguments on purpose: `t_0`, `WN_0` and `WN` are `Union{Nothing,Int}` fields at
every call site, and as keywords they would form a non-concrete `NamedTuple`, heap-allocated
per call, on a path (read once per satellite per epoch) that must stay allocation-free.
"""
function assemble_inter_system_time_offset_parameters(
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
    InterSystemTimeOffsetParameters(
        target,
        A_0,
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
    get_group_delay(state::GNSSDecoderState, signal::AbstractGNSSSignal) -> Union{Nothing,Float64}

Group delay (s) of the ranging `signal` against the reference the broadcast satellite
clock is stated for, from the decoder's validated data. Subtract it from the clock
polynomial's offset `Δt_SV` to get the clock offset on `signal`:
`Δt_SV(signal) = Δt_SV - get_group_delay(state, signal)`. `nothing` until the validated
data holds the set [`has_positioning_data`](@ref) requires and every term of the signal's
correction has been decoded (see below).

!!! note "Combining signals"

    This is each signal's own (single-frequency) correction. A combination of ranges takes
    the same combination of group delays: with `γ = (f₁/f₂)²`, the ionosphere-free pair's
    is `(g₂ - γ·g₁) / (1 - γ)` for the group delays `g₁`, `g₂` returned here, as the ICDs
    state their dual-frequency corrections (BDS-SIS-ICD-B2a-1.0 Eqs. 7-18 to 7-21). It
    vanishes for the pair a clock is stated for (GPS L1/L2 P(Y), Galileo E1/E5b on I/NAV
    and E1/E5a on F/NAV).

`signal` is the GNSSSignals signal instance the range was measured on, as the tracking
loop holds it. It may be a dataless pilot whose navigation data this decoder takes from the
band's data component (GPS L5-Q and L1C-P, Galileo E1-C, E5a-Q and E5b-Q, BeiDou B1C and
B2a pilots); the two need not share a correction, so the signal is asked for rather than
read off the decoder:

| Message                   | Ranging signal       | Group delay                                   |
|:------------------------- |:-------------------- |:--------------------------------------------- |
| GPS LNAV (L1 C/A)         | L1 C/A               | `T_GD`                                        |
| GPS CNAV (L5, L2C)        | L5-I, L5-Q, L2 CM/CL | `T_GD - ISC_L5I5`, `- ISC_L5Q5`, `- ISC_L2C`  |
| GPS CNAV-2 (L1C)          | L1C-D, L1C-P, L1 C/A | `T_GD - ISC_L1CD`, `- ISC_L1CP`, `- ISC_L1CA` |
| Galileo I/NAV (E1-B, E5b) | E1-B/C, E5b-I/Q      | `(f_E1 / f)² · BGD_E1_E5b`                    |
| Galileo F/NAV (E5a)       | E5a-I/Q              | `(f_E1 / f)² · BGD_E1_E5a`                    |
| BeiDou D1/D2 (B1I, B3I)   | B1I, B3I             | `T_GD1`, `0`                                  |
| BeiDou B-CNAV1 (B1C)      | B1C data, pilot      | `T_GD_B1Cp + ISC_B1Cd`, `T_GD_B1Cp`           |
| BeiDou B-CNAV2 (B2a)      | B2a data, pilot      | `T_GD_B2ap + ISC_B2ad`, `T_GD_B2ap`           |
| BeiDou B-CNAV3 (B2b)      | B2b-I                | `T_GD_B2bI`                                   |

GPS states its clock for the L1 P(Y)/L2 P(Y) ionosphere-free pair, with an inter-signal
correction per civil signal (IS-GPS-200 §20.3.3.3.3.2 and §30.3.3.3.1.1, IS-GPS-705
§20.3.3.3.1.2, IS-GPS-800 §3.5.3.9.1). Galileo states it for an ionosphere-free pair (E1/E5b
for I/NAV, E1/E5a for F/NAV); a single-frequency user on the pair's second band scales the
BGD by the squared frequency ratio, 1.79 on E5a and 1.70 on E5b, exactly 1 on E1 (OS SIS
ICD §5.1.5). BeiDou states it for B3I, which needs none, and publishes a group delay per
signal from there, with the data components adding an ISC to their band's pilot
(BDS-SIS-ICD-B1I-3.0 §5.2.4.10, BDS-SIS-ICD-B1C-1.0 §7.6 and the B2a/B2b counterparts).

Most terms ride with the clock and are present once the gate passes, including CNAV-2's
`T_GD`, `ISC_L1CD` and `ISC_L1CP` (subframe 2, with the ephemeris). GPS CNAV's `T_GD` and
ISCs (message type 30), CNAV-2's ISCs of the other civil signals (`ISC_L1CA` and the L2/L5
ones, subframe 3 page 1) and B2a's `T_GD_B2ap` and `ISC_B2ad` (message type 30) arrive on
their own cadence; the gate does not wait for them (see [`has_positioning_data`](@ref)), so
this returns `nothing` until they have. What to do meanwhile is the caller's choice: a
single-frequency fix taking the correction as zero,
`something(get_group_delay(state, signal), 0.0)`, is off by about a metre, and an
ionosphere-free combination scales that by the missing range's coefficient (1.3 to 2.6 for
the GPS and BeiDou pairs), so it may rather wait.

Each message has methods for exactly the ranging signals it broadcasts a correction for
(the table rows), whichever of its signals the decoder is on: an I/NAV decoder on E1-B
answers for E5b ranges too, a CNAV-2 decoder for L1 C/A. Any other pairing, such as an E5a
range on an I/NAV decoder, or any range on Galileo E6-B (no group delay), throws a
`MethodError` rather than returning a plausible number.
"""
function get_group_delay end

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
        state = with_num_bits_since_time_reference(state, num_bits_since_time_reference)
        state, resolved_buffer = complement_buffer_if_necessary(state, buffer)
        state = decode_sync_unit(state, resolved_buffer)
        state = validate_data(state)
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
