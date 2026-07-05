# GPS L2C (CNAV) signal layer — IS-GPS-200N §30.
#
# GPS L2C broadcasts the civil navigation (CNAV) message on the L2 CM code
# (the time-multiplexed L2 CL code is a dataless pilot, IS-GPS-200N §3.2.1.4 /
# §3.3.3.1). The CNAV message on L2C is *bit-for-bit identical* to the one on
# GPS L5I: same 8-bit preamble `10001011`, same rate-1/2 K=7 continuous
# convolutional FEC (G1 = 0o171, G2 = 0o133), same CRC-24Q, same
# per-message-type bit layouts, same π and `TOW × 6` scaling (IS-GPS-200N §30
# ≡ IS-GPS-705J §20.3.3). The signal-layer differences — 25 bps → 50 sps and
# 12-second messages (vs L5's 50 bps → 100 sps, 6-second messages) — are
# purely in the time domain and do not affect this symbol-domain decoder: the
# message is always 600 symbols and the sync window 616. The TOW count is in
# 6-second units on both signals (IS-GPS-200N §30.3.3: "multiplied by 6 …
# start of the next 12-second message"), so the shared `decode_syncro_sequence`
# (`TOW × 6`) applies unchanged.
#
# Consequently this file reuses the entire shared GPS CNAV core in
# `gps/cnav.jl` (FEC, sync, message parsing, and the `GPSCNAVData` container).
# The only L2C-specific behaviour is `is_sat_healthy`, which reports the L2
# signal-health bit (IS-GPS-200N §30.3.3.1.1.2) instead of the L5 bit. This
# mirrors how gnss-sdr and pocketsdr drive L2C and L5 through one CNAV decoder.

"""
    GPSL2CMConstants

GPS L2C specialization of [`GPSCNAVConstants`](@ref GNSSDecoder.GPSCNAVConstants) (`GPSCNAVConstants{:GPSL2CM}`).
Same field values as the GPS L5I constants — the CNAV message is identical on
both signals; the distinct tag only selects the L2 health bit in
[`is_sat_healthy`](@ref). The data-bearing L2C component is the L2 CM code
(`GPSL2CM`); the L2 CL code is a dataless pilot. Reference: IS-GPS-200N
§30.3.2 / §3.3.3.1.
"""
const GPSL2CMConstants = GPSCNAVConstants{:GPSL2CM}

"""
$(TYPEDSIGNATURES)

Create a decoder state for GPS L2C CNAV navigation messages.

Initializes a [`GNSSDecoderState`](@ref) configured for decoding GPS L2C civil
navigation (CNAV) messages from the FEC-encoded 50 sps soft symbols of the L2
CM component. The CNAV message is identical to GPS L5I's, so decoding reuses
the shared GPS CNAV core: each sync attempt Viterbi-decodes the buffered
616-symbol window, locates the 8-bit preamble (`0b10001011`) at both ends of
the decoded bit window, validates the 300-bit message with CRC-24Q, and
dispatches it to per-type parsers (message types 10-15, 30-37, and 40,
IS-GPS-200N §30.3.3). Decoded fields land in a [`GPSCNAVData`](@ref) (the shared
CNAV container).

# Arguments

  - `prn::Int`: Pseudo-Random Noise code identifier (1-63 for GPS satellites)

# Returns

  - `GNSSDecoderState{GPSCNAVData}`: Initialized decoder state for GPS L2C

# Example

```julia
state = GPSL2CMDecoderState(1)           # PRN 1
state = decode(state, soft_symbols, num_symbols)
if is_sat_healthy(state)
    # Use state.data for positioning
end
```

# See Also

  - [`GNSSDecoderState`](@ref): The underlying state structure
  - [`GPSL5IDecoderState`](@ref): The GPS L5I decoder sharing this CNAV core
  - [`decode`](@ref): Decode soft symbols using this state
  - [`reset_decoder_state`](@ref): Reset after signal loss
  - [`is_sat_healthy`](@ref): Check satellite health status
"""
function GPSL2CMDecoderState(prn)
    GNSSDecoderState(
        prn,
        GPSCNAVData(),
        GPSCNAVData(),
        GPSL2CMConstants(),
        GPSCNAVCache(),
        nothing,
        false,
    )
end

# Dispatch from a GNSSSignals system type, mirroring `GNSSDecoderState(::GPSL5I, …)`.
# CNAV rides on the L2 CM (data) component — `GPSL2CM` — while `GPSL2CL` is the
# dataless pilot (IS-GPS-200N §3.3.3.1), so only `GPSL2CM` maps to a decoder.
function GNSSDecoderState(system::GPSL2CM, prn)
    GPSL2CMDecoderState(prn)
end

# L2C-M's CNAV symbol rate; keyed on the constants type so it stays distinct from
# L5-I despite sharing `GPSCNAVData` (see `src/gps/l5.jl`).
GNSSSignals.get_data_frequency(::GNSSDecoderState{<:Any,GPSL2CMConstants}) =
    get_data_frequency(GPSL2CM)

"""
$(TYPEDSIGNATURES)

Check if the GPS L2C satellite is healthy and usable for positioning.

Examines the L2 signal health bit decoded from the most recent message
type 10 (IS-GPS-200N §30.3.3.1.1.2): a satellite is healthy iff the L2 health
bit is 0 (some or all codes and data on the L2 carrier are OK). This is the
only decode-level difference from GPS L5I, which reports the L5 health bit.

!!! warning

    Requires message type 10 to have been decoded and the positioning set to
    have been validated; returns `false` until then.

# Arguments

  - `state::GNSSDecoderState{<:GPSCNAVData,<:GPSL2CMConstants}`: GPS L2C decoder state.

# Returns

  - `Bool`: `true` iff the L2 signal-health bit indicates OK.
"""
function is_sat_healthy(state::GNSSDecoderState{<:GPSCNAVData,<:GPSL2CMConstants})
    state.data.l2_health === false
end
