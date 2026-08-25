# Galileo E1-B (I/NAV) signal layer — Galileo OS SIS ICD, Issue 2.2 §4.3.
#
# The E1-B component of the E1 signal carries the I/NAV message at 250 sps in
# 250-symbol (1 s) page parts. Framing, FEC, sync, the CRC gate, and every
# word-type parser are the shared I/NAV core in `galileo/inav.jl` — the very
# same message is broadcast on E5b-I, which the ICD states uses "the same page
# layout" with only the page sequencing differing (§4.3.1). This file holds only
# what is specific to E1-B: the constants alias, the decoder-state constructors
# (including the BOC(1,1) approximation), the signal mapping, and the health
# check, which reports the E1-B/C facet of word type 5.

"""
    GalileoE1BConstants

Galileo E1-B specialization of [`GalileoINAVConstants`](@ref GNSSDecoder.GalileoINAVConstants)
(`GalileoINAVConstants{:GalileoE1B}`). Same field values as the E5b constants —
the I/NAV message is identical on both components; the distinct tag only selects
the reported signal and the E1-B/C health facet in [`is_sat_healthy`](@ref).
Reference: Galileo OS SIS ICD, Issue 2.2, Table 38 / Table 68.
"""
const GalileoE1BConstants = GalileoINAVConstants{:GalileoE1B}

"""
$(TYPEDSIGNATURES)

Create a decoder state for Galileo E1-B I/NAV navigation messages.

Initializes a [`GNSSDecoderState`](@ref) configured for decoding Galileo E1-B
(Open Service) navigation messages. The decoder extracts ephemeris, clock
correction, ionospheric parameters, and health data from the 250 sps I/NAV
symbol stream using Viterbi decoding, and publishes them in a
[`GalileoINAVData`](@ref GNSSDecoder.GalileoINAVData) (the shared I/NAV
container, also used by the E5b decoder).

# Arguments

  - `prn::Int`: Pseudo-Random Noise code identifier (1-36 for Galileo satellites)

# Returns

  - `GNSSDecoderState{GalileoINAVData}`: Initialized decoder state for Galileo E1-B

# Example

```julia
state = GalileoE1BDecoderState(1)  # Create decoder for PRN 1
state = decode(state, soft_symbols, num_symbols)
if is_sat_healthy(state)
    # Use state.data for positioning
end
```

# See Also

  - [`GNSSDecoderState`](@ref): The underlying state structure
  - [`GalileoE5bDecoderState`](@ref): The Galileo E5b decoder sharing this I/NAV core
  - [`decode`](@ref): Decode soft symbols using this state
  - [`reset_decoder_state`](@ref): Reset after signal loss
  - [`is_sat_healthy`](@ref): Check satellite health status
"""
function GalileoE1BDecoderState(prn)
    GNSSDecoderState(
        prn,
        GalileoINAVData(),
        GalileoINAVData(),
        GalileoE1BConstants(),
        GalileoINAVCache(),
        nothing,
        false,
    )
end

function GNSSDecoderState(system::GalileoE1B, prn)
    GalileoE1BDecoderState(prn)
end

# GalileoE1B_BOC11 is the BOC(1,1) approximation of Galileo E1B — a lower
# sampling-rate replica many software receivers substitute for the full
# CBOC(6,1,1/11) spec (see GNSSSignals.jl `GalileoE1B_BOC11`). It carries the
# *identical* I/NAV navigation message: same 4092-chip primary code, same
# 250 bps data rate, same page/word structure. The modulation difference is a
# tracking/acquisition concern only, so decoding reuses the E1B decoder
# unchanged.
function GNSSDecoderState(system::GalileoE1B_BOC11, prn)
    GalileoE1BDecoderState(prn)
end

# I/NAV rides on both E1-B and E5b-I but decodes into a single
# `GalileoINAVData`, so the signal is keyed on the constants type
# (`GalileoINAVConstants{:GalileoE1B}` vs `{:GalileoE5bI}`) — the only thing that
# tells the two decoders apart. Both the full CBOC E1B and its BOC(1,1)
# approximation report `GalileoE1B`: the approximation is a
# tracking/acquisition concern and the I/NAV stream is identical. Signal
# metadata is forwarded through this mapping (see `src/gps/l1ca.jl`).
get_signal_type(::GalileoE1BConstants) = GalileoE1B

"""
$(TYPEDSIGNATURES)

Check if the Galileo satellite is healthy and usable for positioning on E1-B.

Examines both the signal health status (`E1B_SHS`) and data validity
status (`E1B_DVS`) from I/NAV word type 5. A satellite is
considered healthy only if both conditions are met:

  - Signal health is `signal_ok`
  - Data validity is `navigation_data_valid`

This is the only decode-level difference from the Galileo E5b decoder, which
reports the E5b facet of the same word type 5.

!!! warning

    This function requires that word type 5 has been successfully decoded.
    Check that `state.data.E1B_SHS` is not `nothing` before relying
    on this result.

# Arguments

  - `state::GNSSDecoderState{<:GalileoINAVData,<:GalileoE1BConstants}`: Galileo E1-B decoder state with decoded data

# Returns

  - `Bool`: `true` if satellite health and data validity indicate normal operation

# Example

```julia
state = GalileoE1BDecoderState(1)
state = decode(state, soft_symbols, num_symbols)
if is_sat_healthy(state)
    # Safe to use for positioning
end
```

# See Also

  - [`GalileoE1BDecoderState`](@ref): Create decoder state
  - [`decode`](@ref): Decode navigation data
"""
function is_sat_healthy(state::GNSSDecoderState{<:GalileoINAVData,<:GalileoE1BConstants})
    state.data.E1B_SHS == signal_ok && state.data.E1B_DVS == navigation_data_valid
end
