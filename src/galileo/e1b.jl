# Galileo E1-B (I/NAV) signal layer — Galileo OS SIS ICD, Issue 2.2 §4.3.
#
# E1-B carries the I/NAV message at 250 sps in 250-symbol (1 s) page parts. Framing,
# FEC, sync, CRC gate and word-type parsers are the shared I/NAV core in
# `galileo/inav.jl` (E5b-I broadcasts the same message, §4.3.1). This file holds the
# E1-B constants alias, decoder-state constructors, signal mapping and health check
# (E1-B/C facet of word type 5).

"""
    GalileoE1BConstants

Galileo E1-B specialization of [`GalileoINAVConstants`](@ref GNSSDecoder.GalileoINAVConstants)
(`GalileoINAVConstants{:GalileoE1B}`). Same field values as the E5b constants; the
tag only selects the reported signal and the E1-B/C health facet in
[`is_sat_healthy`](@ref).
Reference: Galileo OS SIS ICD, Issue 2.2, Table 38 / Table 68.
"""
const GalileoE1BConstants = GalileoINAVConstants{:GalileoE1B}

"""
$(TYPEDSIGNATURES)

Create a decoder state for Galileo E1-B I/NAV navigation messages.

Decodes ephemeris, clock correction, ionospheric and health data from the 250 sps
I/NAV symbol stream into a [`GalileoINAVData`](@ref GNSSDecoder.GalileoINAVData)
(shared with the E5b decoder).

# Arguments

  - `prn::Int`: Pseudo-Random Noise code identifier (1-36 for Galileo satellites)

# Returns

  - `GNSSDecoderState{GalileoINAVData}`: Initialized decoder state for Galileo E1-B

# Example

```julia
state = GalileoE1BDecoderState(1)  # Create decoder for PRN 1
state = decode!(state, soft_bits, num_bits)
if is_sat_healthy(state)
    # Use state.data for positioning
end
```

# See Also

  - [`GNSSDecoderState`](@ref): The underlying state structure
  - [`GalileoE5bDecoderState`](@ref): The Galileo E5b decoder sharing this I/NAV core
  - [`decode!`](@ref): Decode soft symbols using this state
  - [`reset_decoder_state!`](@ref): Reset after signal loss
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

# The default BOC(1,1) and full CBOC replicas carry identical I/NAV symbols.
function GNSSDecoderState(::GalileoE1B, prn)
    GalileoE1BDecoderState(prn)
end

# Both E1-B and E5b-I decode into `GalileoINAVData`, so the signal is keyed on the
# constants type (`{:GalileoE1B}` vs `{:GalileoE5bI}`). Both the full CBOC E1B and its
# BOC(1,1) approximation report `GalileoE1B`; the decoder does not retain which
# tracking replica produced the stream (see `src/gnss.jl`).
get_signal_type(::GalileoE1BConstants) = GalileoE1B

"""
$(TYPEDSIGNATURES)

Check if the Galileo satellite is healthy and usable for positioning on E1-B.

Healthy only if both the signal health status (`E1B_SHS`) and data validity status
(`E1B_DVS`) from I/NAV word type 5 are normal:

  - Signal health is `signal_ok`
  - Data validity is `navigation_data_valid`

This is the only decode-level difference from the E5b decoder, which reports the E5b
facet of the same word type.

!!! warning

    A `false` result can also mean word type 5 has not been decoded yet. Check that
    `state.data.E1B_SHS` is not `nothing` to tell the two apart.

# Arguments

  - `state::GNSSDecoderState{<:GalileoINAVData,<:GalileoE1BConstants}`: Galileo E1-B decoder state with decoded data

# Returns

  - `Bool`: `true` if satellite health and data validity indicate normal operation

# Example

```julia
state = GalileoE1BDecoderState(1)
state = decode!(state, soft_bits, num_bits)
if is_sat_healthy(state)
    # Safe to use for positioning
end
```

# See Also

  - [`GalileoE1BDecoderState`](@ref): Create decoder state
  - [`decode!`](@ref): Decode navigation data
"""
function is_sat_healthy(state::GNSSDecoderState{<:GalileoINAVData,<:GalileoE1BConstants})
    state.data.E1B_SHS == signal_ok && state.data.E1B_DVS == navigation_data_valid
end
