# API Reference

## Decoder State

```@docs
GNSSDecoderState
```

## Constructors

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

```@docs
decode
```

## State Management

```@docs
reset_decoder_state
```

## Health Status

```@docs
is_sat_healthy
```

## Positioning Readiness

Pair [`is_decoding_completed_for_positioning`](@ref) with [`is_sat_healthy`](@ref)
to gate use of a satellite in a fix: the first confirms the required navigation
data set has been decoded and validated, the second that the satellite is
broadcasting healthy. See the docstring for what it deliberately does *not*
gate on (ephemeris freshness, second-order corrections, the alert flag).

```@docs
is_decoding_completed_for_positioning
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
| `get_data_frequency` | navigation-message symbol rate | `50 Hz` |
| `get_time_system` / `get_time_system_id` / `get_time_system_name` | time scale the decoded week numbers and times of week are counted in | `GPST()` / `:GPST` / `"GPS Time"` |
| `get_system_start_time` / `get_tai_system_start_time` / `get_tai_offset` | that scale's epoch (as a UTC label, and on the continuous TAI scale) and offset from TAI — what turns a decoded WN/TOW pair into an absolute instant | `1980-01-06T00:00:00` / `1980-01-06T00:00:19` / `19 s` |

Dispatch is on the constants type, which keeps decoders that share a data
container distinct: GPS L5-I and L2C-M both decode into a `GPSCNAVData`, and
Galileo E1-B and E5b-I both into a `GalileoINAVData`, but each reports its own
band, ids and symbol rates.

```julia
using GNSSDecoder, GNSSSignals
get_data_frequency(GPSL5IDecoderState(1))     # 100 Hz
get_data_frequency(GPSL2CMDecoderState(1))    #  50 Hz
get_signal_name(GalileoE5aDecoderState(1))    # "Galileo E5a-I"
get_signal_name(GalileoE5bDecoderState(1))    # "Galileo E5b-I"
get_band_id(GalileoE6BDecoderState(1))        # :E6
get_band_id(GalileoE1BDecoderState(1))        # :L1 — bands are identified by RF
                                              #       frequency, not ICD label
get_system_start_time(GalileoE1BDecoderState(1))  # 1999-08-21T23:59:47
```

Accessors describing the spreading code (`get_code_length`,
`get_code_frequency`, `get_carrier_phase_offset`, …) are not forwarded — they
belong to acquisition and tracking rather than to a symbol-domain decoder — but
remain one step away via [`get_signal_type`](@ref):

```julia
get_code_length(get_signal_type(GPSL1CADecoderState(1)))  # 1023
```

```@docs
get_signal_type
```

## Satellite and Time Metadata

Three questions a positioning engine asks of every decoder, whatever signal it
is on. Each is answered from broadcast data by some signals and from
constellation structure by others, and each is reconstruction work the ICD
defines — so it is done here, once per signal, rather than in every consumer.

| Accessor | Answers | `nothing` when |
|:---|:---|:---|
| [`get_orbit_class`](@ref) | GEO / IGSO / MEO, as an [`OrbitClass`](@ref) | the signal cannot say — see the warning in its docstring |
| [`get_time_of_week`](@ref) | seconds of week at the epoch the message stamps | no time has been decoded and validated yet |
| [`get_time_offset`](@ref) | offset to another constellation's time scale, as a [`GNSSTimeOffset`](@ref) | this satellite is not broadcasting one for that target |

`get_time_of_week` folds away the four spellings the ICDs use — a plain `TOW`
or `SOW` on most signals, `ITOW·7200 + toi·18` on GPS L1C-D, and
`HOW·3600 + soh·18` on BeiDou B1C — and pairs with the state's symbol counter to
give the current time on every signal with one expression:

```julia
using GNSSDecoder, GNSSSignals
tow = get_time_of_week(state)
now = tow + state.num_bits_after_valid_syncro_sequence / get_data_frequency(state)
```

`get_time_offset` folds away the five shapes the same quantity is broadcast in —
Galileo's two-term GPS-only GGTO, the GPS CNAV and CNAV-2 GGTOs tagged by a
`GNSS_ID`/`GGTO_ID` code, BeiDou B2a/B2b's single tagged BGTO set, B1C's set
keyed per system, and D1/D2's three epochless pairs named per system — including
each ICD's "not available" sentinel:

```julia
offset = get_time_offset(state, GPST())     # `GPST()`, `GST()` or `BDT()`
# Δτ is the time since the offset's reference epoch: `t - t_0 + 604800(WN - WN_0)`,
# or simply `t` on BeiDou D1/D2, which broadcasts no epoch (`t_0 === nothing`).
t_gpst = t_bdt - (offset.A_0 + offset.A_1 * Δτ + offset.A_2 * Δτ^2)
```

It also folds in the part that is *not* broadcast. Two time scales differ by a
defined whole-second offset plus a steering residual, and only the residual is
on the air — the bias coefficient is 16 bits at 2⁻³⁵ s, a range of ±0.95 µs,
which cannot express a whole second. `A_0` carries both, so the subtraction
above is true as written. The defined part is zero between GPST and GST (both
TAI − 19 s) and −14 s from BDT to either, so a consumer that tested only Galileo
would never see it.

```@docs
OrbitClass
get_orbit_class
get_time_of_week
GNSSTimeOffset
get_time_offset
```

## Shared Utilities

Signal-independent building blocks used across the decoders (CRC-24Q, the
BCH(51,8) TOI codec, the block (de)interleaver, and the GF(2^8) Reed-Solomon
codec behind Galileo HAS).

These are public but unexported, so they are reached through the module: write
`GNSSDecoder.crc24q(…)`, or import what you need with
`using GNSSDecoder: crc24q`.

```@docs
GNSSDecoder.crc24q
GNSSDecoder.BCHToiSync
GNSSDecoder.sync_bch_toi
GNSSDecoder.pack_hard_codeword
GNSSDecoder.deinterleave!
GNSSDecoder.interleave!
GNSSDecoder.GaloisField256
GNSSDecoder.GALILEO_HAS_GF256
GNSSDecoder.rs_generator_polynomial
GNSSDecoder.rs_systematic_generator_matrix
GNSSDecoder.rs_erasure_decode
```

## Data Types

Every concrete per-signal data type subtypes the abstract supertype of its
constellation, which in turn subtypes `AbstractGNSSData`. The supertypes carry
the facts every signal of a constellation shares, stated once via subtype
dispatch. Three intermediate levels name families within a constellation:
`AbstractGalileoEphemerisData` separates the ephemeris-bearing Galileo messages
from E6-B's corrections-only data, and `AbstractGPSCNAVData` /
`AbstractBeiDouCNAVData` name the modern civil message families (CNAV/CNAV-2
and B-CNAV1/2/3), whose shared quasi-Keplerian ephemeris a consumer's
propagator dispatches on.

```@docs
GNSSDecoder.AbstractGPSData
GNSSDecoder.AbstractGPSCNAVData
GNSSDecoder.AbstractGalileoData
GNSSDecoder.AbstractGalileoEphemerisData
GNSSDecoder.AbstractBeiDouData
GNSSDecoder.AbstractBeiDouCNAVData
```

### GPS L1 C/A

```@docs
GNSSDecoder.GPSL1CAConstants
GPSL1CAData
GPSL1CAAlmanac
```

### Galileo I/NAV (shared by E1-B and E5b)

Galileo E1-B and E5b-I carry the identical I/NAV message — the OS SIS ICD states
they use "the same page layout", differing only in page sequencing — so they
share the decoded [`GalileoINAVData`](@ref GNSSDecoder.GalileoINAVData) container
and one constants struct,
[`GalileoINAVConstants`](@ref GNSSDecoder.GalileoINAVConstants). The per-signal
constants are type aliases that fix its signal tag; they differ only in which
signal-health facet of word type 5 [`is_sat_healthy`](@ref) reports (E1-B/C or
E5b), mirroring GPS CNAV.

```@docs
GNSSDecoder.GalileoINAVConstants
GNSSDecoder.GalileoE1BConstants
GNSSDecoder.GalileoE5bConstants
GNSSDecoder.GalileoINAVData
GalileoReducedCED
GalileoAlmanac
SignalHealth
DataValidityStatus
```

### Galileo E5a

```@docs
GNSSDecoder.GalileoE5aConstants
GNSSDecoder.GalileoE5aData
```

### Galileo E6-B (C/NAV — High Accuracy Service)

Galileo E6-B carries the C/NAV message, which is the Signal-in-Space channel of
the Galileo High Accuracy Service: PPP-grade orbit, clock, code-bias and
phase-bias corrections *to another signal's* broadcast ephemeris, rather than an
ephemeris of its own. A HAS message is spread across up to 32 encoded pages
distributed among satellites and reassembled with a Reed-Solomon erasure decode
(see [Shared Utilities](#Shared-Utilities)), so
[`GalileoE6BData`](@ref) exposes both the latest complete
[`GalileoHASMessage`](@ref) and an accumulated view of the latest of each
correction block.

```@docs
GNSSDecoder.GalileoE6BConstants
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
GNSSDecoder.GPSL1C_DConstants
GPSL1C_DData
GPSL1C_DReducedAlmanac
GPSL1C_DMidiAlmanac
GPSL1C_DDifferentialCorrection
```

### GPS CNAV (shared by L5I and L2C)

GPS L5I and GPS L2C carry the identical CNAV message, so they share the decoded
[`GPSCNAVData`](@ref) container (and its almanac/correction records) and one
constants struct, [`GPSCNAVConstants`](@ref GNSSDecoder.GPSCNAVConstants). The
per-signal constants are type aliases that fix its signal tag; they differ only
in which signal-health bit [`is_sat_healthy`](@ref) reports.

```@docs
GNSSDecoder.GPSCNAVConstants
GNSSDecoder.GPSL5IConstants
GNSSDecoder.GPSL2CMConstants
GPSCNAVData
GPSCNAVReducedAlmanac
GPSCNAVMidiAlmanac
GPSCNAVClockDifferentialCorrection
GPSCNAVEphemerisDifferentialCorrection
GPSCNAVIntegritySupportMessage
```


### BeiDou D1/D2 NAV (shared by B1I and B3I)

BeiDou B1I and B3I carry the identical legacy navigation message — D1 NAV on
MEO/IGSO satellites, D2 NAV on GEO satellites, selected by PRN — so they share
the decoded [`BeiDouDNAVData`](@ref) container and one constants struct,
[`BeiDouDNAVConstants`](@ref GNSSDecoder.BeiDouDNAVConstants). The per-signal
constants are type aliases that fix its signal tag, mirroring GPS CNAV.

```@docs
GNSSDecoder.BeiDouDNAVConstants
GNSSDecoder.BeiDouB1IConstants
GNSSDecoder.BeiDouB3IConstants
BeiDouDNAVData
BeiDouDNAVAlmanac
```

### BeiDou B1C

```@docs
GNSSDecoder.BeiDouB1CConstants
BeiDouB1CData
BeiDouB1CBGTO
```

### BeiDou B2a

```@docs
GNSSDecoder.BeiDouB2aConstants
BeiDouB2aData
```

### BeiDou B2b

```@docs
GNSSDecoder.BeiDouB2bConstants
BeiDouB2bData
```

### BeiDou shared almanac records

The BDS-3 midi and reduced almanac blocks are bit-identical across B-CNAV1
(B1C), B-CNAV2 (B2a), and B-CNAV3 (B2b), so all three decoders produce the
same record types.

```@docs
BeiDouMidiAlmanac
BeiDouReducedAlmanac
```
