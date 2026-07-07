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

`GNSSSignals.get_data_frequency` is extended for [`GNSSDecoderState`](@ref): it
returns the navigation-message symbol rate of the signal the decoder demodulates
(e.g. `50 Hz` for GPS L1 C/A, `100 Hz` for GPS L5-I, `50 Hz` for GPS L2C-M). It
forwards to the corresponding signal's rate in GNSSSignals, so the value stays
single-sourced. Dispatch is on the constants type, which keeps decoders that
share a data container distinct — GPS L5-I and L2C-M both decode into
`GPSCNAVData` but report their own rates.

```julia
using GNSSDecoder, GNSSSignals
get_data_frequency(GPSL5IDecoderState(1))   # 100 Hz
get_data_frequency(GPSL2CMDecoderState(1))  #  50 Hz
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
