# API Reference

## Decoder State

```@docs
GNSSDecoderState
```

## Constructors

Construct a decoder from its GNSSSignals signal and satellite PRN:

```julia
using GNSSDecoder, GNSSSignals
state = GNSSDecoderState(GPSL1CA(), 25)
state = decode!(state, soft_bits)
```

The signal-specific constructors below are conveniences for the same states.

```@docs
GPSL1CADecoderState
GalileoE1BDecoderState
GalileoE5bDecoderState
GalileoE5aDecoderState
GalileoE6BDecoderState
GPSL1C_DDecoderState
GPSL5IDecoderState
GPSL2CMDecoderState
BeiDouB1IDecoderState
BeiDouB3IDecoderState
BeiDouB1CDecoderState
BeiDouB2aDecoderState
BeiDouB2bDecoderState
```

## Decoding

`decode!` overwrites the buffers the decoder state was constructed with and
allocates nothing. Continue with the state it returns; take a `copy` first if
you need a snapshot of the earlier one.

```@docs
decode!
Base.copy(::GNSSDecoderState)
```

## Preallocated Storage

The keyed stores in the decoded data are sized once, when the decoder state is
constructed, so that `decode!` can overwrite them instead of allocating. The
broadcast text messages are inline `CStaticString`s from StaticStrings.jl for
the same reason.

```@docs
GNSSDecoder.SlotDictionary
```

## State Management

```@docs
reset_decoder_state!
```

## Health Status

```@docs
is_sat_healthy
```

## Positioning Readiness

Pair [`has_positioning_data`](@ref) with [`is_sat_healthy`](@ref)
to gate use of a satellite in a fix: the first confirms the required navigation
data set has been decoded and validated, the second that the satellite is
broadcasting healthy. See the docstring for what it deliberately does *not*
gate on (ephemeris freshness, second-order corrections, the alert flag).

```@docs
has_positioning_data
```

## Signal Metadata

A decoder state knows which signal it demodulates, so GNSSSignals' signal
accessors are extended for [`GNSSDecoderState`](@ref) and answer directly from
the state — no need to carry the signal alongside the decoder just to ask what
it is. Each forwards to the corresponding signal in GNSSSignals through
[`get_signal_type`](@ref), so every value stays single-sourced, and each folds
to a compile-time constant.

| Accessor | Answers | Example (`GPSL2CMDecoderState(1)`) |
|:---|:---|:---|
| `get_signal_id` / `get_signal_name` | which signal | `:GPSL2CM` / `"GPS L2CM"` |
| `get_constellation_id` / `get_constellation_name` | which constellation | `:GPS` / `"GPS"` |
| `get_band` / `get_band_id` / `get_band_name` | which RF band | `L2()` / `:L2` / `"L2"` |
| `get_bit_rate` | bit rate as tracked: the channel symbol rate, before error-correction decoding, of the satellite, a BeiDou B1I/B3I GEO satellite's D2 rate where it differs from the signal's | `50 Hz` |
| `get_time_system` / `get_time_system_id` / `get_time_system_name` | time scale the decoded week numbers and times of week are counted in | `GPST()` / `:GPST` / `"GPS Time"` |
| `get_system_start_time` / `get_tai_system_start_time` / `get_tai_offset` | that scale's epoch (as a UTC label, and on the continuous TAI scale) and offset from TAI | `1980-01-06T00:00:00` / `1980-01-06T00:00:19` / `19 s` |
| `get_week_start_epoch(state, week)` | the TAI epoch a week of that scale begins at, a GNSSSignals `TAITime`: plus a decoded time of week, the instant (see [`get_broadcast_week`](@ref)) | `1980-01-06T00:00:19 TAI` for week 0 |

Dispatch is on the constants type, which keeps decoders that share a data
container distinct: GPS L5-I and L2C-M both decode into a `GPSCNAVData`, and
Galileo E1-B and E5b-I both into a `GalileoINAVData`, but each reports its own
band, ids and symbol rates. Tracking replica variants share one decoder:
BOC(1,1) and full CBOC E1B both report the canonical `GalileoE1B` type. Forwarded
signal names and ids identify the broadcast component (`"Galileo E1B"` and
`:GalileoE1B`), independently of the tracking replica. Query GNSSSignals'
`get_replica_id` or `get_replica_name` on the tracking signal to distinguish
replicas; a decoder does not retain or expose replica identity.

```julia
using GNSSDecoder, GNSSSignals
get_bit_rate(GPSL5IDecoderState(1))     # 100 Hz
get_bit_rate(GPSL2CMDecoderState(1))    #  50 Hz
get_signal_name(GalileoE5aDecoderState(1))    # "Galileo E5a-I"
get_signal_name(GalileoE5bDecoderState(1))    # "Galileo E5b-I"
get_band_id(GalileoE6BDecoderState(1))        # :E6
get_band_id(GalileoE1BDecoderState(1))        # :L1 — bands are identified by RF
                                              #       frequency, not ICD label
get_system_start_time(GalileoE1BDecoderState(1))  # 1999-08-21T23:59:47
```

Accessors describing the spreading code (`get_code_length`,
`get_chip_rate`, `get_carrier_phase_offset`, …) are not forwarded — they
belong to acquisition and tracking rather than to a symbol-domain decoder — but
remain one step away via [`get_signal_type`](@ref):

```julia
get_code_length(get_signal_type(GPSL1CADecoderState(1)))  # 1023
```

```@docs
get_signal_type
```

## Satellite and Time Metadata

The questions a positioning engine asks of every decoder, whatever signal it
is on. Each is answered from broadcast data by some signals and from
constellation structure by others, and each is reconstruction work the ICD
defines — so it is done here, once per signal, rather than in every consumer.

| Accessor | Answers | `nothing` when |
|:---|:---|:---|
| [`get_orbit_class`](@ref) | GEO / IGSO / MEO, as an [`OrbitClass`](@ref) | the signal cannot say — see the warning in its docstring |
| [`get_broadcast_time_of_week`](@ref) | seconds of week at the epoch the message stamps | no time has been decoded and validated yet |
| [`get_broadcast_week`](@ref) | absolute week number, the GPS L1 C/A rollover resolved | no week number has been decoded and validated yet |
| [`get_broadcast_time_offset`](@ref) | offset to another constellation's time scale at a time of week, in seconds | this satellite is not broadcasting one for that target |
| [`get_broadcast_time_offset_epoch`](@ref) | the reference epoch that offset is extrapolated from, to judge its age by | no such offset, or one without a reference epoch (BeiDou D1/D2) |

`get_broadcast_time_of_week` folds away the four spellings the ICDs use — a plain `TOW`
or `SOW` on most signals, `ITOW·7200 + toi·18` on GPS L1C-D, and
`HOW·3600 + soh·18` on BeiDou B1C — and pairs with the state's symbol counter in `get_decoded_time_of_week`, which
returns the uncorrected time at the last consumed symbol:

```julia
using GNSSDecoder, GNSSSignals
using Unitful: Hz, ustrip
time_of_week = get_decoded_time_of_week(state)
```

`get_broadcast_time_offset` folds away the five shapes the same quantity is
broadcast in — Galileo's two-term GPS-only GGTO, the GPS CNAV and CNAV-2 GGTOs tagged by a
`GNSS_ID`/`GGTO_ID` code, BeiDou B2a/B2b's single tagged BGTO set, B1C's set
keyed per system, and D1/D2's three epochless pairs named per system — including
each ICD's "not available" sentinel:

```julia
Δt = get_broadcast_time_offset(state, GPST(), t_bdt)  # `GPST()`, `GST()` or `BDT()`
scale_offset = ustrip(s, get_time_scale_offset(BDT(), GPST()))  # 14 s
t_gpst = t_bdt - Δt + scale_offset
```

It evaluates the broadcast polynomial at the given time of week against the
decoder's own week number, which the truncated reference week is resolved
against, so a caller needs only the time of week.

`Δt` comes from the coefficients as broadcast. Two time scales differ by the leap
seconds between them plus a steering residual, and only the residual is on the
air — the bias coefficient is 16 bits at 2⁻³⁵ s, a range of ±0.95 µs (on BeiDou
D1/D2 14 bits at 0.1 ns), which cannot express a whole second. The leap seconds follow from the scales'
definitions, GNSSSignals' `get_time_scale_offset`: zero between GPST and GST
(both TAI − 19 s) and 14 s from BDT to either, so a conversion that is only
ever tried on Galileo never notices them missing.

```@docs
OrbitClass
get_orbit_class
get_broadcast_time_of_week
get_decoded_time_of_week
get_broadcast_week
get_broadcast_time_offset
BroadcastTimeOffsetParameters
get_broadcast_time_offset_parameters
get_broadcast_time_offset_epoch
```

## Broadcast Ephemeris

The decoded data keeps every ephemeris field under its ICD name, which differs
between the messages that state the orbit in Keplerian form (GPS LNAV, Galileo
I/NAV and F/NAV, BeiDou D1/D2) and those that state it quasi-Keplerian, as
deltas off reference values (GPS CNAV and CNAV-2, BeiDou B-CNAV).
[`get_broadcast_ephemeris`](@ref) folds both into one record of plain numbers,
with the reference values applied and the system's constants alongside, so one
propagator covers every signal:

```julia
using GNSSDecoder
ephemeris = get_broadcast_ephemeris(state)  # `nothing` until decoding is complete
```

```@docs
BroadcastEphemeris
get_broadcast_ephemeris
BroadcastEphemerisValidity
get_broadcast_ephemeris_validity
```

## Group Delay

The broadcast satellite clock is stated for a reference signal or signal pair,
and a range on any other signal needs that signal's group delay subtracted from
the clock offset. Which broadcast terms apply, with which sign and scale,
depends on the message and on the ranging signal, which may be a dataless pilot
whose navigation data the decoder took from the band's data component — so the
signal instance the range was measured on is passed in:

```julia
using GNSSDecoder, GNSSSignals
signal = GalileoE5aQ()  # as the tracking loop holds it
Δt_SV_E5a = Δt_SV - get_group_delay(state, signal)  # (f_E1/f_E5a)² · BGD_E1_E5a
```

A combination of ranges, such as the ionosphere-free one, takes the same
combination of their group delays. The terms GPS CNAV, BeiDou B2a and, for the
other civil signals, CNAV-2 broadcast apart from the clock can still be missing
once a fix is possible;
until they arrive this returns `nothing`, and the caller chooses whether to
take the correction as zero or to wait.

```@docs
get_group_delay
```

## Range Accuracy

Every message broadcasts the accuracy of the satellite's ranging signal as an
index — GPS's user range accuracy (URA), Galileo's signal-in-space accuracy
(SISA), BeiDou's URA — and each ICD states the metres an index stands for.
[`get_range_accuracy`](@ref) returns those metres, `Inf` where the satellite
broadcasts no accuracy prediction, so a positioning engine can weight its
measurements by it and a navigation file can record it.

```@docs
get_range_accuracy
```

## Ionospheric Parameters

Each constellation broadcasts the coefficients of its own ionospheric model:
Klobuchar on GPS (every civil message), BeiDou's Klobuchar variant on D1/D2,
NeQuick-G on Galileo I/NAV and F/NAV, and BDGIM on BeiDou B-CNAV. Each accessor
answers on every decoder, with `nothing` where the message carries no such set
or has not decoded it yet. The sets are not part of the positioning readiness
gate; the models that evaluate them belong to the positioning engine.

```@docs
KlobucharParameters
BeiDouKlobucharParameters
NeQuickGParameters
BDGIMParameters
get_klobuchar_parameters
get_nequick_g_parameters
get_bdgim_parameters
```

## UTC, Leap Seconds and GGTO

The time parameters a navigation file records beside the ephemeris: the
broadcast offset of the system's time from UTC, its leap seconds, and Galileo's
offset from GPS time. Every message truncates their reference weeks — to 6, 8 or
13 bits — so each accessor lifts them to the week nearest the decoder's own
[`get_broadcast_week`](@ref), on the same count,
and a consumer needs no field widths. The UTC and leap-second accessors take
`get_broadcast_week`'s `approximate_year` for that reason; Galileo's weeks are full, so the
GGTO accessor needs none.

| Accessor | Answers on | `nothing` when |
|:---|:---|:---|
| [`get_utc_parameters`](@ref) | GPS LNAV, CNAV, CNAV-2; Galileo I/NAV, F/NAV; BeiDou D1/D2 (without a reference epoch) and B-CNAV | not yet decoded, or on Galileo E6-B |
| [`get_leap_second_parameters`](@ref) | the same | not yet decoded, or on Galileo E6-B |
| [`get_ggto_parameters`](@ref) | Galileo I/NAV and F/NAV | not yet decoded, marked unavailable, or on any other decoder |

```julia
using GNSSDecoder
utc = get_utc_parameters(state; approximate_year = 2026)
if !isnothing(utc)
    utc.WN_0t  # a full week number, not the 8 or 13 bits broadcast
end
```

```@docs
UTCParameters
LeapSecondParameters
GGTOParameters
get_utc_parameters
get_leap_second_parameters
get_ggto_parameters
```

## Data Types

Every concrete message data type subtypes the abstract supertype of its
constellation, so a caller can dispatch on the constellation of a decoder's data
without enumerating its signals.

```@docs
GNSSDecoder.AbstractGPSData
GNSSDecoder.AbstractGalileoData
GNSSDecoder.AbstractBeiDouData
```

### GPS L1 C/A

```@docs
GPSLNAVData
GNSSDecoder.GPSL1CAAlmanac
```

### Galileo I/NAV (shared by E1-B and E5b)

Galileo E1-B and E5b-I carry the identical I/NAV message — the OS SIS ICD states
they use "the same page layout", differing only in page sequencing — so they
share the decoded [`GalileoINAVData`](@ref GNSSDecoder.GalileoINAVData)
container. The two decoders differ only in which signal-health facet of word
type 5 [`is_sat_healthy`](@ref) reports (E1-B/C or E5b), mirroring GPS
CNAV.

```@docs
GalileoINAVData
GNSSDecoder.GalileoReducedCED
GNSSDecoder.GalileoAlmanac
SignalHealth
DataValidityStatus
```

### Galileo E5a

```@docs
GalileoFNAVData
```

### Galileo E6-B (C/NAV — High Accuracy Service)

Galileo E6-B carries the C/NAV message, which is the Signal-in-Space channel of
the Galileo High Accuracy Service: PPP-grade orbit, clock, code-bias and
phase-bias corrections *to another signal's* broadcast ephemeris, rather than an
ephemeris of its own. A HAS message is spread across up to 32 encoded pages
distributed among satellites and reassembled with a Reed-Solomon erasure decode,
so
[`GalileoE6BData`](@ref) exposes both the latest complete
[`GalileoHASMessage`](@ref) and an accumulated view of the latest of each
correction block.

```@docs
GalileoE6BData
GalileoHASMessage
GalileoHASMask
GalileoHASSatelliteMask
GalileoHASCorrectionBlock
GalileoHASOrbitCorrection
GalileoHASClockCorrection
GalileoHASCodeBias
GalileoHASPhaseBias
HASStatus
```

### GPS L1C-D

```@docs
GPSCNAV2Data
GNSSDecoder.GPSL1C_DReducedAlmanac
GNSSDecoder.GPSL1C_DMidiAlmanac
GNSSDecoder.GPSL1C_DDifferentialCorrection
```

### GPS CNAV (shared by L5I and L2C)

GPS L5I and GPS L2C carry the identical CNAV message, so they share the decoded
[`GPSCNAVData`](@ref) container (and its almanac/correction records). The two
decoders differ only in which signal-health bit [`is_sat_healthy`](@ref)
reports.

```@docs
GPSCNAVData
GNSSDecoder.GPSCNAVReducedAlmanac
GNSSDecoder.GPSCNAVMidiAlmanac
GNSSDecoder.GPSCNAVClockDifferentialCorrection
GNSSDecoder.GPSCNAVEphemerisDifferentialCorrection
GNSSDecoder.GPSCNAVIntegritySupportMessage
```


### BeiDou D1/D2 NAV (shared by B1I and B3I)

BeiDou B1I and B3I carry the identical legacy navigation message — D1 NAV on
MEO/IGSO satellites, D2 NAV on GEO satellites, selected by PRN — so they share
the decoded [`BeiDouDNAVData`](@ref) container, mirroring GPS CNAV.

```@docs
BeiDouDNAVData
GNSSDecoder.BeiDouDNAVAlmanac
```

### BeiDou B1C

```@docs
BeiDouCNAV1Data
GNSSDecoder.BeiDouB1CBGTO
```

### BeiDou B2a

```@docs
BeiDouCNAV2Data
```

### BeiDou B2b

```@docs
BeiDouCNAV3Data
```

### BeiDou shared almanac records

The BDS-3 midi and reduced almanac blocks are bit-identical across B-CNAV1
(B1C), B-CNAV2 (B2a), and B-CNAV3 (B2b), so all three decoders produce the
same record types.

```@docs
GNSSDecoder.BeiDouMidiAlmanac
GNSSDecoder.BeiDouReducedAlmanac
```
