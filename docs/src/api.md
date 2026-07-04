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
