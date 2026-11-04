# BeiDou B1I signal layer — BDS-SIS-ICD-B1I-3.0.
#
# B1I broadcasts the legacy navigation message on the 1561.098 MHz B1 carrier: D1
# (50 bps, NH20 secondary code) on MEO/IGSO satellites, D2 (500 bps, no secondary code)
# on GEO satellites (§5.1.1). Word coding, subframe layouts, and scalings are shared
# with B3I through `beidou/dnav.jl`.
#
# The only per-signal difference is the group delay: the clock parameters are referenced
# to B3I, so a single-frequency B1I user corrects with TGD1
# ((Δt_sv)_B1I = Δt_sv − TGD1, BDS-SIS-ICD-B1I-3.0 §5.2.4.10), while a B3I user applies
# them directly (BDS-SIS-ICD-B3I-1.0 §5.2.4.10). Both TGDs are decoded; applying them is
# the positioning consumer's job.

"""
    BeiDouB1IConstants

BeiDou B1I specialization of [`BeiDouDNAVConstants`](@ref GNSSDecoder.BeiDouDNAVConstants)
(`BeiDouDNAVConstants{:BeiDouB1I}`). Same field values as the B3I constants; the tag
only selects the signal reported by [`get_signal_type`](@ref).
Reference: BDS-SIS-ICD-B1I-3.0 §5.
"""
const BeiDouB1IConstants = BeiDouDNAVConstants{:BeiDouB1I}

"""
$(TYPEDSIGNATURES)

Create a decoder state for BeiDou B1I legacy navigation messages (D1/D2 NAV).

The format is selected from the PRN: D1 for MEO/IGSO (PRN 6-58; 50 bps after NH20
wipe-off), D2 for GEO (PRN 1-5 and 59-63; 500 bps, no secondary code). Each 300-bit
subframe is synchronized via the 11-bit preamble `11100010010`, BCH(15,11,1)-decoded
per word, and parsed into a [`BeiDouDNAVData`](@ref) (shared with B3I).

# Arguments

  - `prn::Int`: Pseudo-Random Noise code identifier (1-63 for BeiDou satellites)

# Returns

  - `GNSSDecoderState{BeiDouDNAVData}`: Initialized decoder state for BeiDou B1I

# Example

```julia
state = BeiDouB1IDecoderState(20)         # PRN 20 (MEO/IGSO ⇒ D1)
state = decode!(state, soft_bits, num_bits)
if is_sat_healthy(state)
    # Use state.data for positioning
end
```

# See Also

  - [`GNSSDecoderState`](@ref): The underlying state structure
  - [`BeiDouB3IDecoderState`](@ref): The B3I decoder sharing this legacy NAV core
  - [`decode!`](@ref): Decode soft symbols using this state
  - [`reset_decoder_state!`](@ref): Reset after signal loss
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

# The signal this decoder demodulates, keyed on the constants type, which distinguishes
# the two decoders sharing `BeiDouDNAVData`.
get_signal_type(::BeiDouB1IConstants) = BeiDouB1I
