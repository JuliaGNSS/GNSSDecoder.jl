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

### GPS L5I

```@docs
GNSSDecoder.GPSL5IConstants
GPSL5IData
GPSL5IReducedAlmanac
GPSL5IMidiAlmanac
GPSL5IClockDifferentialCorrection
GPSL5IEphemerisDifferentialCorrection
GPSL5IIntegritySupportMessage
```
