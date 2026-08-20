# BeiDou B1I signal layer — BDS-SIS-ICD-B1I-3.0.
#
# B1I broadcasts the legacy navigation message on the 1561.098 MHz B1
# carrier: format D1 (50 bps, under the NH20 secondary code) on MEO/IGSO
# satellites and format D2 (500 bps, no secondary code) on GEO satellites
# (§5.1.1). The message core — BCH(15,11,1) word coding, subframe layouts,
# and all parameter scalings — is shared with B3I through `beidou/dnav.jl`,
# exactly as the two ICDs mirror each other's §5.
#
# The only substantive per-signal difference is which group delay applies:
# the broadcast clock parameters are referenced to the B3I frequency, so the
# single-frequency B1I user must correct the satellite clock offset with
# TGD1 ((Δt_sv)_B1I = Δt_sv − TGD1, BDS-SIS-ICD-B1I-3.0 §5.2.4.10), while a
# B3I user applies the clock parameters directly (BDS-SIS-ICD-B3I-1.0
# §5.2.4.10 broadcasts TGD1/TGD2 for the *other* signals' users). Both TGDs
# are decoded either way; applying them is the positioning consumer's job.

"""
    BeiDouB1IConstants

BeiDou B1I specialization of [`BeiDouDNAVConstants`](@ref GNSSDecoder.BeiDouDNAVConstants)
(`BeiDouDNAVConstants{:BeiDouB1I}`). Same field values as the B3I constants —
the legacy D1/D2 message is identical on both signals; the distinct tag only
selects the signal identity reported by [`get_signal_type`](@ref).
Reference: BDS-SIS-ICD-B1I-3.0 §5.
"""
const BeiDouB1IConstants = BeiDouDNAVConstants{:BeiDouB1I}

"""
$(TYPEDSIGNATURES)

Create a decoder state for BeiDou B1I legacy navigation messages (D1/D2 NAV).

Initializes a [`GNSSDecoderState`](@ref) configured for decoding the BeiDou
legacy navigation message from B1I soft symbols. For MEO/IGSO satellites
(PRN 6-58) that is the D1 message: 50 bps data bits after wipe-off of the
NH20 secondary code. For GEO satellites (PRN 1-5 and 59-63) it is the D2
message at 500 bps (no secondary code); the decoder selects the format from
the PRN. Each 300-bit subframe is synchronized via the 11-bit preamble
`11100010010`, BCH(15,11,1)-decoded per word, and parsed into a
[`BeiDouDNAVData`](@ref) (the container shared with B3I).

# Arguments

  - `prn::Int`: Pseudo-Random Noise code identifier (1-63 for BeiDou satellites)

# Returns

  - `GNSSDecoderState{BeiDouDNAVData}`: Initialized decoder state for BeiDou B1I

# Example

```julia
state = BeiDouB1IDecoderState(20)         # PRN 20 (MEO/IGSO ⇒ D1)
state = decode(state, soft_symbols, num_symbols)
if is_sat_healthy(state)
    # Use state.data for positioning
end
```

# See Also

  - [`GNSSDecoderState`](@ref): The underlying state structure
  - [`BeiDouB3IDecoderState`](@ref): The B3I decoder sharing this legacy NAV core
  - [`decode`](@ref): Decode soft symbols using this state
  - [`reset_decoder_state`](@ref): Reset after signal loss
  - [`is_sat_healthy`](@ref): Check satellite health status
"""
function BeiDouB1IDecoderState(prn)
    1 <= prn <= 63 || throw(ArgumentError("BeiDou PRN must be in 1..63"))
    GNSSDecoderState(
        prn,
        BeiDouDNAVData(),
        BeiDouDNAVData(),
        BeiDouB1IConstants(),
        BeiDouDNAVCache(),
        nothing,
        false,
    )
end

# Dispatch from the GNSSSignals system type, mirroring `GNSSDecoderState(::GPSL1CA, …)`.
function GNSSDecoderState(system::BeiDouB1I, prn)
    BeiDouB1IDecoderState(prn)
end

# The signal this decoder demodulates. Keyed on the constants type, which is
# 1:1 with the signal and tells apart the two decoders sharing
# `BeiDouDNAVData` (see `src/beidou/b3i.jl`).
get_signal_type(::BeiDouB1IConstants) = BeiDouB1I
