# GPS L2C (CNAV) signal layer — IS-GPS-200N §30.
#
# CNAV is broadcast on the L2 CM code (the time-multiplexed L2 CL code is a dataless
# pilot, IS-GPS-200N §3.2.1.4 / §3.3.3.1). The message is bit-for-bit identical to
# GPS L5I: same preamble `10001011`, rate-1/2 K=7 FEC (G1 = 0o171, G2 = 0o133),
# CRC-24Q, message-type layouts and `TOW × 6` scaling (IS-GPS-200N §30 ≡
# IS-GPS-705J §20.3.3). L2C differs only in the time domain (25 bps → 50 sps,
# 12-second messages vs L5's 50 bps → 100 sps, 6-second), which does not affect this
# symbol-domain decoder: messages are 600 symbols, the sync window 616. TOW counts
# 6-second units on both (IS-GPS-200N §30.3.3), so `decode_sync_unit` applies as is.
#
# This file therefore reuses the shared CNAV core in `gps/cnav.jl`. The only
# L2C-specific behaviour is `is_sat_healthy`, which reports the L2 signal-health bit
# (IS-GPS-200N §30.3.3.1.1.2) instead of the L5 bit.

"""
    GPSL2CMConstants

GPS L2C specialization of [`GPSCNAVConstants`](@ref GNSSDecoder.GPSCNAVConstants)
(`GPSCNAVConstants{:GPSL2CM}`). Same field values as the L5I constants; the tag only selects
the L2 health bit in [`is_sat_healthy`](@ref). The data component is the L2 CM code
(`GPSL2CM`); the L2 CL code is a dataless pilot. Reference: IS-GPS-200N §30.3.2 / §3.3.3.1.
"""
const GPSL2CMConstants = GPSCNAVConstants{:GPSL2CM}

"""
$(TYPEDSIGNATURES)

Create a decoder state for GPS L2C CNAV navigation messages.

Initializes a [`GNSSDecoderState`](@ref) for decoding GPS L2C CNAV messages from
the FEC-encoded 50 sps soft symbols of the L2 CM component. The message is
identical to GPS L5I's, so the shared CNAV core is reused: each sync attempt
Viterbi-decodes the buffered 616-symbol window, locates the 8-bit preamble
(`0b10001011`) at both ends, validates the 300-bit message with CRC-24Q, and
dispatches to per-type parsers (message types 10-15, 30-37, 40, IS-GPS-200N
§30.3.3). Fields land in a [`GPSCNAVData`](@ref).

# Arguments

  - `prn::Int`: Pseudo-Random Noise code identifier (1-63 for GPS satellites)

# Returns

  - `GNSSDecoderState{GPSCNAVData}`: Initialized decoder state for GPS L2C

# Example

```julia
state = GPSL2CMDecoderState(1)           # PRN 1
state = decode!(state, soft_bits, num_bits)
if is_sat_healthy(state)
    # Use state.data for positioning
end
```

# See Also

  - [`GNSSDecoderState`](@ref): The underlying state structure
  - [`GPSL5IDecoderState`](@ref): The GPS L5I decoder sharing this CNAV core
  - [`decode!`](@ref): Decode soft symbols using this state
  - [`reset_decoder_state!`](@ref): Reset after signal loss
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

# Dispatch from a GNSSSignals system type. Only `GPSL2CM` maps to a decoder;
# `GPSL2CL` is the dataless pilot (IS-GPS-200N §3.3.3.1).
function GNSSDecoderState(system::GPSL2CM, prn)
    GPSL2CMDecoderState(prn)
end

# Keyed on the constants type so L2C-M stays distinct from L5-I despite sharing
# `GPSCNAVData` (see `src/gps/l5.jl`, `src/gps/l1ca.jl`).
get_signal_type(::GPSL2CMConstants) = GPSL2CM

"""
$(TYPEDSIGNATURES)

Check if the GPS L2C satellite is healthy and usable for positioning.

Examines the L2 health bit of the most recent message type 10 (IS-GPS-200N
§30.3.3.1.1.2): healthy iff it is 0 (some or all codes and data on L2 are OK).
GPS L5I reports the L5 bit instead.

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
