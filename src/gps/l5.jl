# GPS L5I (CNAV) signal layer — IS-GPS-705J §20.3.
#
# The GPS L5 in-phase component carries the CNAV message at 50 bps, FEC-encoded
# to 100 sps over 6-second 300-bit messages (IS-GPS-705J §3.3.3.1.1, §20.3).
# The message format, FEC, CRC, sync, and per-message-type parsing are the
# shared GPS CNAV core in `gps/cnav.jl` (the very same message is broadcast on
# GPS L2C, IS-GPS-200N §30). This file holds only what is specific to L5I: the
# constants type, the decoder-state constructor, and the health check (which
# reports the L5 signal-health bit).

"""
    GPSL5IConstants

GPS L5I specialization of [`GPSCNAVConstants`](@ref GNSSDecoder.GPSCNAVConstants) (`GPSCNAVConstants{:GPSL5I}`).
Same field values as the GPS L2C constants — the distinct tag only selects the
L5 health bit in [`is_sat_healthy`](@ref). Reference: IS-GPS-705J §20.3.3 /
§20.3.4.3.
"""
const GPSL5IConstants = GPSCNAVConstants{:GPSL5I}

"""
$(TYPEDSIGNATURES)

Create a decoder state for GPS L5I CNAV navigation messages.

Initializes a [`GNSSDecoderState`](@ref) configured for decoding GPS L5I
civil navigation (CNAV) messages from FEC-encoded 100 sps soft symbols. Each
sync attempt Viterbi-decodes the buffered 616-symbol window, locates the
8-bit preamble (`0b10001011`) at both ends of the decoded bit window,
validates the 300-bit message with CRC-24Q, and dispatches it to per-type
parsers (message types 10-15, 30-37, and 40, IS-GPS-705J §20.3.3).

# Arguments

  - `prn::Int`: Pseudo-Random Noise code identifier (1-63 for GPS satellites)

# Returns

  - `GNSSDecoderState{GPSCNAVData}`: Initialized decoder state for GPS L5I

# Example

```julia
state = GPSL5IDecoderState(1)            # PRN 1
state = decode(state, soft_symbols, num_symbols)
if is_sat_healthy(state)
    # Use state.data for positioning
end
```

# See Also

  - [`GNSSDecoderState`](@ref): The underlying state structure
  - [`GPSL2CMDecoderState`](@ref): The GPS L2C decoder sharing this CNAV core
  - [`decode`](@ref): Decode soft symbols using this state
  - [`reset_decoder_state`](@ref): Reset after signal loss
  - [`is_sat_healthy`](@ref): Check satellite health status
"""
function GPSL5IDecoderState(prn)
    GNSSDecoderState(
        prn,
        GPSCNAVData(),
        GPSCNAVData(),
        GPSL5IConstants(),
        GPSCNAVCache(),
        nothing,
        false,
    )
end

function GNSSDecoderState(system::GPSL5I, prn)
    GPSL5IDecoderState(prn)
end

"""
$(TYPEDSIGNATURES)

Check if the GPS L5 satellite is healthy and usable for positioning.

Examines the L5 signal health bit decoded from the most recent message
type 10 (IS-GPS-705J §20.3.3.1.1.2): a satellite is healthy iff the health
bit is 0 (all navigation data on the L5 signal are OK).

!!! warning

    Requires message type 10 to have been decoded and the positioning set to
    have been validated; returns `false` until then.

# Arguments

  - `state::GNSSDecoderState{<:GPSCNAVData,<:GPSL5IConstants}`: GPS L5I decoder state.

# Returns

  - `Bool`: `true` iff the L5 signal-health bit indicates OK.
"""
function is_sat_healthy(state::GNSSDecoderState{<:GPSCNAVData,<:GPSL5IConstants})
    state.data.l5_health === false
end
