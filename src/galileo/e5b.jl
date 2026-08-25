# Galileo E5b (I/NAV) signal layer — Galileo OS SIS ICD, Issue 2.2 §4.3.
#
# I/NAV is broadcast on the E5b-I (in-phase, data) component of the E5b
# sideband at 1207.14 MHz (Table 2) as well as on E1-B; E5b-Q is a dataless
# pilot (§2.3.1.1, Table 5 — "no data ('pilot component')"). The ICD is explicit that the two I/NAV components "use
# the same page layout since the service provided on these frequencies is a
# dual frequency service, using frequency diversity. Only page sequencing is
# different" (§4.3.1) — same 250 sps, same 250-symbol page parts, same 10-bit
# page sync `0101100000`, same K=7 rate-1/2 NSC FEC over the same 30×8 block
# interleaver, same even/odd pairing into 128-bit words, same CRC-24Q over the
# 196 protected bits of a page pair.
#
# The page-part layouts differ only in fields this decoder never reads
# (Table 38): where E1-B's odd part carries OSNMA(40) + SAR(22) + Spare(2),
# E5b-I carries one 64-bit "Reserved 1" field, and where E1-B ends in the 8-bit
# SSP, E5b-I has "Reserved 2". Both are 64 + 8 bits, both leave the CRC-protected
# prefix at 82 bits and the CRC at odd-part bits 83-106, and neither trailing
# 8-bit field is CRC-protected. So the symbol-domain decoder is unchanged and
# this file reuses the whole shared I/NAV core in `galileo/inav.jl`, exactly as
# `gps/l2c.jl` reuses the GPS CNAV core.
#
# The one sequencing consequence: word types 16 (Reduced CED) and 17-20 (FEC2
# RS CED) are E1-B only (Table 40), so `state.data.reduced_ced` stays empty on
# E5b. Word types 1-6 (ephemeris, clock, iono, health, GST-UTC) and 7-10 (the
# almanac chain, offset half a constellation from E1-B's, Table 41) all appear.
#
# The only decode-level difference from E1-B is which facet of word type 5
# `is_sat_healthy` reports — E5b signal health and data validity instead of
# E1-B/C's.

"""
    GalileoE5bConstants

Galileo E5b specialization of [`GalileoINAVConstants`](@ref GNSSDecoder.GalileoINAVConstants)
(`GalileoINAVConstants{:GalileoE5bI}`). Same field values as the E1-B constants —
the I/NAV message is identical on both components; the distinct tag only selects
the reported signal (`GalileoE5bI`, the data-bearing component; `GalileoE5bQ` is
the dataless pilot) and the E5b health facet in [`is_sat_healthy`](@ref).
Reference: Galileo OS SIS ICD, Issue 2.2, Table 38 / Table 68.
"""
const GalileoE5bConstants = GalileoINAVConstants{:GalileoE5bI}

"""
$(TYPEDSIGNATURES)

Create a decoder state for Galileo E5b I/NAV navigation messages.

Initializes a [`GNSSDecoderState`](@ref) configured for decoding the Galileo
I/NAV message from the FEC-encoded 250 sps soft symbols of the E5b-I component.
The I/NAV message is identical to E1-B's, so decoding reuses the shared I/NAV
core: the 10-bit page-sync pattern is matched at both ends of the 260-symbol
window, the 240 encoded symbols are 30×8 deinterleaved and Viterbi-decoded to a
114-bit page part, consecutive even + odd parts are stitched into a 128-bit
nominal word, and the pair is gated on CRC-24Q before its word type is parsed.
Decoded fields land in a [`GalileoINAVData`](@ref GNSSDecoder.GalileoINAVData)
(the shared I/NAV container).

# Arguments

  - `prn::Int`: Pseudo-Random Noise code identifier (1-36 for Galileo satellites)

# Returns

  - `GNSSDecoderState{GalileoINAVData}`: Initialized decoder state for Galileo E5b

# Example

```julia
state = GalileoE5bDecoderState(1)  # Create decoder for PRN 1
state = decode(state, soft_symbols, num_symbols)
if is_sat_healthy(state)
    # Use state.data for positioning
end
```

# See Also

  - [`GNSSDecoderState`](@ref): The underlying state structure
  - [`GalileoE1BDecoderState`](@ref): The Galileo E1-B decoder sharing this I/NAV core
  - [`decode`](@ref): Decode soft symbols using this state
  - [`reset_decoder_state`](@ref): Reset after signal loss
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

# Dispatch from a GNSSSignals system type, mirroring `GNSSDecoderState(::GalileoE1B, …)`.
# I/NAV rides on the E5b-I (data) component — `GalileoE5bI` — while `GalileoE5bQ`
# is the dataless pilot (OS SIS ICD Table 5), so only `GalileoE5bI` maps to a
# decoder.
function GNSSDecoderState(system::GalileoE5bI, prn)
    GalileoE5bDecoderState(prn)
end

# The signal this decoder demodulates — the E5b-I data component, keyed on the
# constants type so it stays distinct from E1-B despite sharing every line of
# decoding and one `GalileoINAVData` container (see `src/galileo/e1b.jl`). Signal
# metadata is forwarded through it, so the two decoders report their own band and
# ids (E5b vs E1) while agreeing on the 250 sps symbol rate.
get_signal_type(::GalileoE5bConstants) = GalileoE5bI

"""
$(TYPEDSIGNATURES)

Check if the Galileo satellite is healthy and usable for positioning on E5b.

Examines both the E5b signal health status (`E5b_SHS`) and E5b data
validity status (`E5b_DVS`) from I/NAV word type 5
(OS SIS ICD Tables 81 and 84). A satellite is considered healthy only if both
conditions are met:

  - Signal health is `signal_ok`
  - Data validity is `navigation_data_valid`

This is the only decode-level difference from the Galileo E1-B decoder, which
reports the E1-B/C facet of the same word type 5. Word type 5 carries both
facets, so a decoder on either component always has both available.

!!! warning

    This function requires that word type 5 has been successfully decoded.
    Check that `state.data.E5b_SHS` is not `nothing` before relying
    on this result.

# Arguments

  - `state::GNSSDecoderState{<:GalileoINAVData,<:GalileoE5bConstants}`: Galileo E5b decoder state with decoded data

# Returns

  - `Bool`: `true` if satellite health and data validity indicate normal operation

# Example

```julia
state = GalileoE5bDecoderState(1)
state = decode(state, soft_symbols, num_symbols)
if is_sat_healthy(state)
    # Safe to use for positioning
end
```

# See Also

  - [`GalileoE5bDecoderState`](@ref): Create decoder state
  - [`decode`](@ref): Decode navigation data
"""
function is_sat_healthy(state::GNSSDecoderState{<:GalileoINAVData,<:GalileoE5bConstants})
    state.data.E5b_SHS == signal_ok && state.data.E5b_DVS == navigation_data_valid
end
