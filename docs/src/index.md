# GNSSDecoder.jl

A Julia package for decoding GNSS (Global Navigation Satellite System) navigation messages.

## Supported Systems

- **GPS L1 C/A**: Decodes the 50 bps LNAV data stream from GPS L1 civil signals
- **Galileo E1B**: Decodes the 250 bps I/NAV data stream from Galileo E1B Open Service signals

## Installation

```julia
using Pkg
Pkg.add("GNSSDecoder")
```

## Quick Start

### GPS L1 Decoding

Initialize a decoder and process bits from your tracking loop:

```jldoctest gps_example
julia> using GNSSDecoder

julia> state = GPSL1DecoderState(1);  # Initialize decoder for PRN 1

julia> state.prn  # Access PRN
1

julia> typeof(state)
GNSSDecoderState{GNSSDecoder.GPSL1Data, GNSSDecoder.GPSL1Constants, GNSSDecoder.GPSL1Cache, GNSSDecoder.UInt320}
```

Process incoming bits and check the decoder state:

```jldoctest gps_example
julia> state = decode(state, UInt8(0b10001011), 8);  # Decode 8 bits

julia> state.num_bits_buffered  # Bits are now buffered
8
```

In a real application, you would decode bits from a tracking loop:

```julia
for i in 1:iterations
    # Track signal (e.g., with Tracking.jl)
    track_res = track(signal, track_state, state.prn, sampling_freq)
    track_state = get_state(track_res)

    # Decode navigation message
    state = decode(state, get_bits(track_res), get_num_bits(track_res))
end

# After decoding completes, access the data
if !isnothing(state.data.TOW)
    println("Time of Week: $(state.data.TOW)")
end
```

### Galileo E1B Decoding

```jldoctest galileo_example
julia> using GNSSDecoder

julia> state = GalileoE1BDecoderState(1);  # Initialize decoder for PRN 1

julia> state.prn
1

julia> typeof(state)
GNSSDecoderState{GNSSDecoder.GalileoE1BData, GNSSDecoder.GalileoE1BConstants, GNSSDecoder.GalileoE1BCache, GNSSDecoder.UInt288}

julia> state = decode(state, UInt16(0b0101100000), 10);  # Decode 10 bits (preamble)

julia> state.num_bits_buffered
10
```

## State Management

### Resetting After Signal Loss

If signal tracking is lost and reacquired, use [`reset_decoder_state`](@ref) to clear
buffers while preserving previously decoded ephemeris:

```jldoctest reset_example
julia> using GNSSDecoder

julia> state = GPSL1DecoderState(1);

julia> state = decode(state, UInt8(0xff), 8);  # Some decoding

julia> state.num_bits_buffered
8

julia> state = reset_decoder_state(state);  # Reset after signal loss

julia> state.num_bits_buffered  # Buffers are cleared
0

julia> state.prn  # PRN is preserved
1
```

### Checking Satellite Health

```jldoctest health_example
julia> using GNSSDecoder

julia> state = GPSL1DecoderState(1);

julia> is_sat_healthy(state)  # Health not yet decoded
false

julia> state = GalileoE1BDecoderState(1);

julia> is_sat_healthy(state)  # Health not yet decoded
false
```

## Data Fields

### GPS L1 Data

After successful decoding, `state.data` contains:

| Field | Description |
|-------|-------------|
| `TOW` | Time of Week (seconds) |
| `trans_week` | Transmission week number |
| `svhealth` | Satellite health status |
| `t_0e`, `t_0c` | Reference times for ephemeris and clock |
| `e` | Eccentricity |
| `sqrt_A` | Square root of semi-major axis |
| `M_0` | Mean anomaly at reference time |
| `Ω_0`, `ω` | Longitude of ascending node, argument of perigee |
| `i_0`, `i_dot` | Inclination and rate |
| `Δn`, `Ω_dot` | Mean motion difference, rate of right ascension |
| `C_rs`, `C_rc`, `C_us`, `C_uc`, `C_is`, `C_ic` | Harmonic correction terms |
| `a_f0`, `a_f1`, `a_f2` | Clock correction coefficients |
| `T_GD` | Group delay differential |

### Galileo E1B Data

Similar ephemeris and clock parameters are available for Galileo, plus:

| Field | Description |
|-------|-------------|
| `WN` | Week number |
| `signal_health_e1b` | E1B signal health status |
| `data_validity_status_e1b` | Data validity status |
| `broadcast_group_delay_e1_e5a` | E1-E5a group delay |
| `broadcast_group_delay_e1_e5b` | E1-E5b group delay |
