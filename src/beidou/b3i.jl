# BeiDou B3I signal layer — BDS-SIS-ICD-B3I-1.0.
#
# B3I broadcasts the same legacy navigation message as B1I on the
# 1268.52 MHz B3 carrier: format D1 (50 bps, NH20 secondary code) on
# MEO/IGSO satellites, format D2 (500 bps) on GEO satellites. The B3I ICD's
# §5 mirrors BDS-SIS-ICD-B1I-3.0 §5 word for word — same BCH(15,11,1) word
# coding, same subframe/page layouts, same parameter scalings — so this file
# reuses the entire shared legacy NAV core in `beidou/dnav.jl` and the
# `BeiDouDNAVData` container.
#
# Group-delay semantics (BDS-SIS-ICD-B3I-1.0 §5.2.4.9/.10): the broadcast
# clock correction parameters are referenced to the B3I signal itself, so the
# single-frequency B3I user applies them directly and TGD1/TGD2 (equipment
# group delay differentials of B1I/B2I relative to B3I) are *not* applied on
# B3I. Both TGDs are decoded regardless; applying the right one is the
# positioning consumer's job.

"""
    BeiDouB3IConstants

BeiDou B3I specialization of [`BeiDouDNAVConstants`](@ref GNSSDecoder.BeiDouDNAVConstants)
(`BeiDouDNAVConstants{:BeiDouB3I}`). Same field values as the B1I constants —
the legacy D1/D2 message is identical on both signals; the distinct tag only
selects the signal identity reported by [`get_signal_type`](@ref).
Reference: BDS-SIS-ICD-B3I-1.0 §5.
"""
const BeiDouB3IConstants = BeiDouDNAVConstants{:BeiDouB3I}

"""
$(TYPEDSIGNATURES)

Create a decoder state for BeiDou B3I legacy navigation messages (D1/D2 NAV).

Initializes a [`GNSSDecoderState`](@ref) configured for decoding the BeiDou
legacy navigation message from B3I soft symbols. The message is identical to
B1I's (see [`BeiDouB1IDecoderState`](@ref)), so decoding reuses the shared
legacy NAV core: D1 (50 bps, after NH20 wipe-off) for MEO/IGSO PRNs, D2
(500 bps) for GEO PRNs, 11-bit preamble sync, per-word BCH(15,11,1)
correction, and parsing into the shared [`BeiDouDNAVData`](@ref) container.

# Arguments

  - `prn::Int`: Pseudo-Random Noise code identifier (1-63 for BeiDou satellites)

# Returns

  - `GNSSDecoderState{BeiDouDNAVData}`: Initialized decoder state for BeiDou B3I

# Example

```julia
state = BeiDouB3IDecoderState(30)         # PRN 30 (MEO/IGSO ⇒ D1)
state = decode(state, soft_symbols, num_symbols)
if is_sat_healthy(state)
    # Use state.data for positioning
end
```

# See Also

  - [`GNSSDecoderState`](@ref): The underlying state structure
  - [`BeiDouB1IDecoderState`](@ref): The B1I decoder sharing this legacy NAV core
  - [`decode`](@ref): Decode soft symbols using this state
  - [`reset_decoder_state`](@ref): Reset after signal loss
  - [`is_sat_healthy`](@ref): Check satellite health status
"""
function BeiDouB3IDecoderState(prn)
    1 <= prn <= 63 || throw(ArgumentError("BeiDou PRN must be in 1..63"))
    GNSSDecoderState(
        prn,
        BeiDouDNAVData(),
        BeiDouDNAVData(),
        BeiDouB3IConstants(),
        BeiDouDNAVCache(),
        nothing,
        false,
    )
end

# Dispatch from the GNSSSignals system type.
function GNSSDecoderState(system::BeiDouB3I, prn)
    BeiDouB3IDecoderState(prn)
end

# The signal this decoder demodulates — keyed on the constants type so it
# stays distinct from B1I despite sharing `BeiDouDNAVData`.
get_signal_type(::BeiDouB3IConstants) = BeiDouB3I
