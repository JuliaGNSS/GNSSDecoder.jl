# API Reference

## Decoder State

```@docs
GNSSDecoderState
```

## Constructors

```@docs
GPSL1CADecoderState
GalileoE1BDecoderState
GalileoE5aDecoderState
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
| `get_system_start_time` / `get_tai_offset` | that scale's epoch and offset from TAI — what turns a decoded WN/TOW pair into an absolute instant | `1980-01-06T00:00:00` / `19 s` |

Dispatch is on the constants type, which keeps decoders that share a data
container distinct: GPS L5-I and L2C-M both decode into a `GPSCNAVData` but
report their own band, ids and symbol rates.

```julia
using GNSSDecoder, GNSSSignals
get_data_frequency(GPSL5IDecoderState(1))     # 100 Hz
get_data_frequency(GPSL2CMDecoderState(1))    #  50 Hz
get_signal_name(GalileoE5aDecoderState(1))    # "Galileo E5a-I"
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

## Shared Utilities

Signal-independent building blocks used across the decoders (CRC-24Q, the
BCH(51,8) TOI codec, and the block (de)interleaver).

```@docs
crc24q
BCHToiSync
sync_bch_toi
soft_to_hard_codeword
pack_hard_codeword
deinterleave!
interleave!
```

## Data Types

Every concrete per-signal data type subtypes the abstract supertype of its
constellation, which in turn subtypes `AbstractGNSSData`. The supertypes carry
the facts every signal of a constellation shares (e.g. the Galileo
ephemeris/clock completeness checks), stated once via subtype dispatch.

```@docs
GNSSDecoder.AbstractGPSData
GNSSDecoder.AbstractGalileoData
GNSSDecoder.AbstractBeiDouData
```

### GPS L1 C/A

```@docs
GNSSDecoder.GPSL1CAConstants
GNSSDecoder.GPSL1CAData
```

### Galileo E1B

```@docs
GNSSDecoder.GalileoE1BConstants
GNSSDecoder.GalileoE1BData
GNSSDecoder.GalileoAlmanac
GNSSDecoder.SignalHealth
GNSSDecoder.DataValidityStatus
```

### Galileo E5a

```@docs
GNSSDecoder.GalileoE5aConstants
GNSSDecoder.GalileoE5aData
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
