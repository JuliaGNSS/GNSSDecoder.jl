# Galileo E5b (I/NAV) signal layer — Galileo OS SIS ICD, Issue 2.2 §4.3.
#
# I/NAV is broadcast on E5b-I (data component, 1207.14 MHz, Table 2) as well as E1-B;
# E5b-Q is a dataless pilot (§2.3.1.1, Table 5). The two components "use the same page
# layout ... Only page sequencing is different" (§4.3.1): same 250 sps, page parts,
# 10-bit page sync `0101100000`, K=7 NSC FEC over the same 30×8 interleaver, even/odd
# pairing into 128-bit words and CRC-24Q over the 196 protected bits.
#
# Page-part layouts differ only in fields this decoder never reads (Table 38): E1-B's
# odd part carries OSNMA(40) + SAR(22) + Spare(2) and ends in the 8-bit SSP, E5b-I has
# "Reserved 1" (64) and "Reserved 2" (8). The CRC-protected prefix (82 bits) and CRC
# (odd-part bits 83-106) are identical, and the trailing 8 bits are not protected. The
# shared I/NAV core in `galileo/inav.jl` is therefore reused whole, as `gps/l2c.jl`
# reuses the GPS CNAV core.
#
# Sequencing consequence: word types 16 (Reduced CED) and 17-20 (FEC2 RS CED) are E1-B
# only (Table 40), so `state.data.reduced_ced` stays empty on E5b. Word types 1-6 and
# 7-10 (almanac chain, offset half a constellation from E1-B's, Table 41) all appear.
#
# The only decode-level difference from E1-B is the facet of word type 5 that
# `is_sat_healthy` reports.

"""
    GalileoE5bConstants

Galileo E5b specialization of [`GalileoINAVConstants`](@ref GNSSDecoder.GalileoINAVConstants)
(`GalileoINAVConstants{:GalileoE5bI}`). Same field values as the E1-B constants; the
tag only selects the reported signal (`GalileoE5bI`; `GalileoE5bQ` is the dataless
pilot) and the E5b health facet in [`is_sat_healthy`](@ref).
Reference: Galileo OS SIS ICD, Issue 2.2, Table 38 / Table 68.
"""
const GalileoE5bConstants = GalileoINAVConstants{:GalileoE5bI}

"""
$(TYPEDSIGNATURES)

Create a decoder state for Galileo E5b I/NAV navigation messages.

Decodes the I/NAV message from the 250 sps soft symbols of E5b-I with the shared I/NAV
core (identical to E1-B's: page sync at both ends of the 260-symbol window, 30×8
deinterleave and Viterbi to a 114-bit page part, even + odd stitching, CRC-24Q gate).
Fields land in a [`GalileoINAVData`](@ref GNSSDecoder.GalileoINAVData).

# Arguments

  - `prn::Int`: Pseudo-Random Noise code identifier (1-36 for Galileo satellites)

# Returns

  - `GNSSDecoderState{GalileoINAVData}`: Initialized decoder state for Galileo E5b

# Example

```julia
state = GalileoE5bDecoderState(1)  # Create decoder for PRN 1
state = decode!(state, soft_bits, num_bits)
if is_sat_healthy(state)
    # Use state.data for positioning
end
```

# See Also

  - [`GNSSDecoderState`](@ref): The underlying state structure
  - [`GalileoE1BDecoderState`](@ref): The Galileo E1-B decoder sharing this I/NAV core
  - [`decode!`](@ref): Decode soft symbols using this state
  - [`reset_decoder_state!`](@ref): Reset after signal loss
  - [`is_sat_healthy`](@ref): Check satellite health status
"""
function GalileoE5bDecoderState(prn)
    GNSSDecoderState(
        prn,
        GalileoINAVData(),
        GalileoINAVData(),
        GalileoE5bConstants(),
        GalileoINAVCache(),
        nothing,
        false,
    )
end

# Only `GalileoE5bI` maps to a decoder; `GalileoE5bQ` is the dataless pilot (OS SIS ICD
# Table 5).
function GNSSDecoderState(system::GalileoE5bI, prn)
    GalileoE5bDecoderState(prn)
end

# Keyed on the constants type to stay distinct from E1-B (see `e1b.jl`); metadata is
# forwarded through it, so each decoder reports its own band and ids.
get_signal_type(::GalileoE5bConstants) = GalileoE5bI

"""
$(TYPEDSIGNATURES)

Check if the Galileo satellite is healthy and usable for positioning on E5b.

Healthy only if both the E5b signal health status (`E5b_SHS`) and data validity status
(`E5b_DVS`) from I/NAV word type 5 are normal (OS SIS ICD Tables 81 and 84):

  - Signal health is `signal_ok`
  - Data validity is `navigation_data_valid`

This is the only decode-level difference from the E1-B decoder, which reports the
E1-B/C facet of the same word type; word type 5 carries both facets.

!!! warning

    A `false` result can also mean word type 5 has not been decoded yet. Check that
    `state.data.E5b_SHS` is not `nothing` to tell the two apart.

# Arguments

  - `state::GNSSDecoderState{<:GalileoINAVData,<:GalileoE5bConstants}`: Galileo E5b decoder state with decoded data

# Returns

  - `Bool`: `true` if satellite health and data validity indicate normal operation

# Example

```julia
state = GalileoE5bDecoderState(1)
state = decode!(state, soft_bits, num_bits)
if is_sat_healthy(state)
    # Safe to use for positioning
end
```

# See Also

  - [`GalileoE5bDecoderState`](@ref): Create decoder state
  - [`decode!`](@ref): Decode navigation data
"""
function is_sat_healthy(state::GNSSDecoderState{<:GalileoINAVData,<:GalileoE5bConstants})
    state.data.E5b_SHS == signal_ok && state.data.E5b_DVS == navigation_data_valid
end
