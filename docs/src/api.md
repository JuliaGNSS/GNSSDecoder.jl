# API Reference

## Decoder State

```@docs
GNSSDecoderState
```

## Constructors

```@docs
GPSL1DecoderState
GalileoE1BDecoderState
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

## Data Types

### GPS L1

```@docs
GNSSDecoder.GPSL1Constants
GNSSDecoder.GPSL1Data
```

### Galileo E1B

```@docs
GNSSDecoder.GalileoE1BConstants
GNSSDecoder.GalileoE1BData
GNSSDecoder.SignalHealth
GNSSDecoder.DataValidityStatus
```
