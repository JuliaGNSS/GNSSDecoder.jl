# GNSSDecoder.jl

A Julia package for decoding GNSS (Global Navigation Satellite System) navigation messages.

## Supported Systems

- **GPS L1 C/A**: Decodes the 50 bps LNAV data stream from GPS L1 civil signals
- **GPS L1C-D**: Decodes the 100 sps CNAV-2 data stream from the modernized GPS L1C signal's data component
- **Galileo E1B**: Decodes the 250 bps I/NAV data stream from Galileo E1B Open Service signals

## Installation

```julia
using Pkg
Pkg.add("GNSSDecoder")
```

## Quick Start

### GPS L1 C/A Decoding

Initialize a decoder and process soft symbols from your tracking loop. The
decoder consumes `Float32` soft symbols where the sign carries the bit
decision (positive ⇒ bit 0, negative ⇒ bit 1) and the magnitude carries
confidence (AFF3CT LLR convention):

```jldoctest gps_example
julia> using GNSSDecoder

julia> state = GPSL1CADecoderState(1);  # Initialize decoder for PRN 1

julia> state.prn  # Access PRN
1

julia> typeof(state)
GNSSDecoderState{GNSSDecoder.GPSL1CAData, GNSSDecoder.GPSL1CAConstants, GNSSDecoder.GPSL1CACache}
```

Process incoming soft symbols and check the decoder state:

```jldoctest gps_example
julia> state = decode(state, Float32[+1, -1, -1, -1, +1, -1, +1, +1], 8);  # Decode 8 soft symbols

julia> GNSSDecoder.num_bits_buffered(state)  # Symbols are now buffered
8
```

In a real application, you would decode soft symbols from a tracking loop. With
`Tracking.jl` v2, take them from `get_soft_bits`, which returns the
polarity-corrected, amplitude-weighted soft bits for the tracked satellite:

```julia
for i in 1:iterations
    # Track signal (e.g., with Tracking.jl)
    track_state = track!(measurement, track_state)

    # Soft symbols for this satellite (Float32; sign = bit, magnitude = confidence)
    soft_symbols = get_soft_bits(track_state, state.prn)

    # Decode navigation message
    state = decode(state, soft_symbols, length(soft_symbols))
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
GNSSDecoderState{GNSSDecoder.GalileoE1BData, GNSSDecoder.GalileoE1BConstants, GNSSDecoder.GalileoE1BCache}

julia> state = decode(state, Float32[+1, -1, +1, +1, -1, -1, -1, -1, -1, +1], 10);  # Decode 10 soft symbols

julia> GNSSDecoder.num_bits_buffered(state)
10
```

### GPS L1C-D Decoding

The GPS L1C-D (CNAV-2) decoder synchronises on the BCH-encoded TOI counter
(no fixed preamble), LDPC-decodes subframes 2 and 3, and validates each with
CRC-24Q. Construction loads the LDPC parity-check matrices shipped with the
package:

```jldoctest l1cd_example
julia> using GNSSDecoder

julia> state = GPSL1C_DDecoderState(1);  # Initialize decoder for PRN 1

julia> state.prn
1

julia> typeof(state)
GNSSDecoderState{GPSL1C_DData, GNSSDecoder.GPSL1C_DConstants, GNSSDecoder.GPSL1C_DCache}

julia> state = decode(state, Float32[+1, -1, +1, +1, -1, -1, -1, -1, -1, +1], 10);  # Decode 10 soft symbols

julia> GNSSDecoder.num_bits_buffered(state)
10
```

## State Management

### Resetting After Signal Loss

If signal tracking is lost and reacquired, use [`reset_decoder_state`](@ref) to clear
buffers while preserving previously decoded ephemeris:

```jldoctest reset_example
julia> using GNSSDecoder

julia> state = GPSL1CADecoderState(1);

julia> state = decode(state, Float32[+1, +1, +1, +1, +1, +1, +1, +1], 8);  # Some decoding

julia> GNSSDecoder.num_bits_buffered(state)
8

julia> state = reset_decoder_state(state);  # Reset after signal loss

julia> GNSSDecoder.num_bits_buffered(state)  # Buffers are cleared
0

julia> state.prn  # PRN is preserved
1
```

### Checking Satellite Health

```jldoctest health_example
julia> using GNSSDecoder

julia> state = GPSL1CADecoderState(1);

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

### GPS L1C-D Data

CNAV-2 clock-and-ephemeris data plus the subframe-3 page payloads — see
[`GPSL1C_DData`](@ref) for the full field list:

| Field | Description |
|-------|-------------|
| `toi`, `ITOW`, `WN` | Time of interval, interval time of week, week number |
| `t_0e`, `ΔA`, `e`, `M_0`, `ω`, `Ω_0`, `i_0`, … | Clock and ephemeris (CED) parameters |
| `α0..α3`, `β0..β3` | Klobuchar ionospheric coefficients (subframe-3 page 1) |
| `A0_UTC`, `Δt_LS`, … | UTC parameters (page 1) |
| `A0_GGTO`, `t_GGTO`, … | GPS/GNSS time offset and EOP (page 2) |
| `reduced_almanacs`, `midi_almanacs` | Per-SV almanac dictionaries (pages 3/4) |
| `differential_corrections` | Per-SV differential corrections (page 5) |
| `text_message` | Broadcast text (page 6) |
