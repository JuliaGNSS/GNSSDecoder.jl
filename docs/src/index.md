# GNSSDecoder.jl

A Julia package for decoding GNSS (Global Navigation Satellite System) navigation messages.

## Supported Systems

- **GPS L1 C/A**: Decodes the 50 bps LNAV data stream from GPS L1 civil signals
- **GPS L1C-D**: Decodes the 100 sps CNAV-2 data stream from the modernized GPS L1C signal's data component
- **GPS L2C**: Decodes the 50 sps CNAV data stream from the GPS L2 CM (civil-moderate) signal component
- **GPS L5I**: Decodes the 100 sps CNAV data stream from the GPS L5 in-phase signal component
- **Galileo E1B**: Decodes the 250 bps I/NAV data stream from Galileo E1B Open Service signals
- **Galileo E5b**: Decodes the same 250 bps I/NAV data stream from the Galileo E5b in-phase (data) component
- **Galileo E5a**: Decodes the 50 sps F/NAV data stream from the Galileo E5a in-phase (data) component
- **Galileo E6B**: Decodes the 1000 sps C/NAV data stream from the Galileo E6-B component — the High Accuracy Service (HAS) orbit, clock and bias corrections
- **BeiDou B1I / B3I**: Decode the legacy D1 (50 bps, MEO/IGSO) and D2 (500 bps, GEO) NAV messages
- **BeiDou B1C**: Decodes the 100 sps B-CNAV1 data stream from the B1C data component
- **BeiDou B2a**: Decodes the 200 sps B-CNAV2 data stream from the B2a data component
- **BeiDou B2b**: Decodes the 1000 sps B-CNAV3 data stream from the B2b_I signal

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
GNSSDecoderState{GPSL1CAData, GNSSDecoder.GPSL1CAConstants, GNSSDecoder.GPSL1CACache}
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
GNSSDecoderState{GalileoINAVData, GNSSDecoder.GalileoINAVConstants{:GalileoE1B}, GNSSDecoder.GalileoINAVCache}

julia> state = decode(state, Float32[+1, -1, +1, +1, -1, -1, -1, -1, -1, +1], 10);  # Decode 10 soft symbols

julia> GNSSDecoder.num_bits_buffered(state)
10
```

### Galileo E5b Decoding

E5b-I carries the identical I/NAV message as E1-B (OS SIS ICD §4.3.1: "the same
page layout … only page sequencing is different"), so the decoder is the same
one, differing only in the signal it reports and in which health facet of word
type 5 [`is_sat_healthy`](@ref) reads:

```jldoctest galileo_e5b_example
julia> using GNSSDecoder, GNSSSignals

julia> state = GalileoE5bDecoderState(1);

julia> get_signal_name(state)
"Galileo E5b-I"

julia> get_band_id(state)
:E5b
```

### Galileo E6B (High Accuracy Service) Decoding

E6-B's C/NAV message carries the Galileo High Accuracy Service: orbit, clock,
code-bias and phase-bias corrections *to another signal's* broadcast ephemeris.
Each satellite broadcasts one 1000-symbol page per second, and a HAS message is
assembled from up to 32 such pages via a Reed-Solomon erasure decode — so a
single satellite may take up to 32 seconds per message, while a receiver pooling
pages from several E6-B satellites completes them much faster.

```jldoctest galileo_e6b_example
julia> using GNSSDecoder

julia> state = GalileoE6BDecoderState(1);

julia> state.constants.syncro_sequence_length   # 1000 symbols = 1 s per page
1000

julia> is_decoding_completed_for_positioning(state)   # corrections, not an ephemeris
false
```

Once a message completes, the corrections are in `state.data`:

```julia
# The most recently completed HAS message, exactly as the ICD defines it
message = state.data.message

# ...or the accumulated latest of each content block, which is usually what a
# correction consumer wants (HAS splits masks/orbits from clocks across messages)
for correction in state.data.orbit_corrections.corrections
    correction.GNSS_ID == 2 || continue          # 0 = GPS, 2 = Galileo
    isnothing(correction.δ_radial) && continue   # "data not available" sentinel
    @show correction.SVID, correction.IOD_ref, correction.δ_radial
end
```

`δ_radial`, `δ_in_track` and `δ_cross_track` are in the satellite-centred NTW
frame and must be rotated into ECEF before use (HAS SIS ICD §7.2); each block
carries its own `validity_interval` in seconds from the message's `TOH`, and the
`mask_id` / `IOD_set_id` that tie it to a satellite set and a broadcast
ephemeris issue.

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

### GPS L5I Decoding

The GPS L5I (CNAV) decoder consumes the 100 sps FEC-encoded channel symbols.
The rate-1/2 K=7 convolutional FEC runs continuously across message
boundaries, so each sync attempt Viterbi-decodes the buffered 616-symbol
window, looks for the 8-bit preamble at both ends of the decoded window, and
validates the 300-bit message with CRC-24Q:

```jldoctest l5i_example
julia> using GNSSDecoder

julia> state = GPSL5IDecoderState(1);  # Initialize decoder for PRN 1

julia> state.prn
1

julia> typeof(state)
GNSSDecoderState{GPSCNAVData, GNSSDecoder.GPSCNAVConstants{:GPSL5I}, GNSSDecoder.GPSCNAVCache}

julia> state = decode(state, Float32[+1, -1, +1, +1, -1, -1, -1, -1, -1, +1], 10);  # Decode 10 soft symbols

julia> GNSSDecoder.num_bits_buffered(state)
10
```

### GPS L2C Decoding

GPS L2C broadcasts the *same* CNAV message as GPS L5I (IS-GPS-200N §30), on the
L2 CM component at 50 sps. Decoding therefore reuses the shared GPS CNAV core,
and decoded fields land in the same [`GPSCNAVData`](@ref) container; only the
health check differs (it reports the L2 signal-health bit):

```jldoctest l2c_example
julia> using GNSSDecoder

julia> state = GPSL2CMDecoderState(1);  # Initialize decoder for PRN 1

julia> state.prn
1

julia> typeof(state)
GNSSDecoderState{GPSCNAVData, GNSSDecoder.GPSCNAVConstants{:GPSL2CM}, GNSSDecoder.GPSCNAVCache}

julia> state = decode(state, Float32[+1, -1, +1, +1, -1, -1, -1, -1, -1, +1], 10);  # Decode 10 soft symbols

julia> GNSSDecoder.num_bits_buffered(state)
10
```

### BeiDou Decoding

All five BeiDou open-service signals are decoded through the same API. B1I and
B3I carry the identical legacy message — D1 NAV (50 bps) on MEO/IGSO
satellites and D2 NAV (500 bps) on GEO satellites, selected automatically by
PRN — so they share the [`BeiDouDNAVData`](@ref) container the way GPS L5I and
L2C share `GPSCNAVData`:

```jldoctest b1i_example
julia> using GNSSDecoder

julia> state = BeiDouB1IDecoderState(20);  # PRN 20: MEO/IGSO, D1 NAV

julia> typeof(state)
GNSSDecoderState{BeiDouDNAVData, GNSSDecoder.BeiDouDNAVConstants{:BeiDouB1I}, GNSSDecoder.BeiDouDNAVCache}

julia> state = decode(state, Float32[+1, -1, +1, +1, -1, -1, -1, -1, -1, +1], 10);  # Decode 10 soft symbols

julia> GNSSDecoder.num_bits_buffered(state)
10
```

The modernized BDS-3 signals B1C (B-CNAV1), B2a (B-CNAV2), and B2b (B-CNAV3)
are coded with non-binary LDPC codes over GF(2⁶). They are decoded here through
the exact binary image of the ICD's parity-check matrix (`data/bcnv*.alist`, see
`scripts/generate_beidou_alist.jl`) with binary belief propagation, which makes
the code definition exact but the decoder weaker than a non-binary one: expect
of the order of a dB less sensitivity than an FFT-QSPA decoder would give,
visible as a raised frame-erasure rate at low C/N₀ rather than as bad data (a
failed decode is dropped by the CRC gate). The belief-propagation stage is also
scale-sensitive, so feed confidence-weighted soft symbols on a roughly LLR-like
scale (`≈ 2·r/σ²`) for best sensitivity — see [`decode`](@ref):

```jldoctest b2a_example
julia> using GNSSDecoder

julia> state = BeiDouB2aDecoderState(30);  # Initialize decoder for PRN 30

julia> typeof(state)
GNSSDecoderState{BeiDouB2aData, GNSSDecoder.BeiDouB2aConstants, GNSSDecoder.BeiDouB2aCache}

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
| `WN` | Transmission week number (modulo 1024) |
| `sv_health` | Raw 6-bit satellite health word (0 = healthy) |
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

### Galileo I/NAV Data (E1B and E5b)

Similar ephemeris and clock parameters are available for Galileo, plus:

| Field | Description |
|-------|-------------|
| `WN` | Week number |
| `E1B_SHS` / `E5b_SHS` | E1-B/C and E5b signal health status |
| `E1B_DVS` / `E5b_DVS` | Data validity status per component |
| `BGD_E1_E5a` | E1-E5a group delay |
| `BGD_E1_E5b` | E1-E5b group delay |
| `almanacs` | Per-SV almanac dictionary (word types 7-10) |
| `reduced_ced` | Reduced clock and ephemeris data (word type 16, E1-B only) |

### Galileo E6B Data (High Accuracy Service)

C/NAV carries corrections rather than an ephemeris — see
[`GalileoE6BData`](@ref) for the full field list:

| Field | Description |
|-------|-------------|
| `HAS_status` | HAS service status (test / operational / do not use) |
| `message` | The most recently completed [`GalileoHASMessage`](@ref), header included (`TOH`, `mask_id`, `IOD_set_id`) |
| `masks` | Every received [`GalileoHASMask`](@ref), keyed by Mask ID |
| `orbit_corrections` | Latest orbit corrections (radial / in-track / cross-track) |
| `clock_corrections`, `clock_subset_corrections` | Latest clock corrections |
| `code_biases`, `phase_biases` | Latest code and phase biases per satellite/signal cell |

### GPS L1C-D Data

CNAV-2 clock-and-ephemeris data plus the subframe-3 page payloads — see
[`GPSL1C_DData`](@ref) for the full field list:

| Field | Description |
|-------|-------------|
| `toi`, `ITOW`, `WN` | Time of interval, interval time of week, week number |
| `t_0e`, `ΔA`, `e`, `M_0`, `ω`, `Ω_0`, `i_0`, … | Clock and ephemeris (CED) parameters |
| `α_0..α_3`, `β_0..β_3` | Klobuchar ionospheric coefficients (subframe-3 page 1) |
| `A_0UTC`, `Δt_LS`, … | UTC parameters (page 1) |
| `A_0GGTO`, `t_GGTO`, … | GPS/GNSS time offset and EOP (page 2) |
| `reduced_almanacs`, `midi_almanacs` | Per-SV almanac dictionaries (pages 3/4) |
| `differential_corrections` | Per-SV differential corrections (page 5) |
| `text_message` | Broadcast text (page 6) |
