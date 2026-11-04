# Definitions shared across GPS signals — the counterpart of
# `galileo/galileo.jl` and `beidou/beidou.jl`. Included before the per-signal
# GPS decoders.

# WGS-84 constants (IS-GPS-200N §20.3.3.4.3; identical in IS-GPS-705J and
# IS-GPS-800J). Each GPS `*Constants` struct exposes them as field defaults.
const GPS_μ = 3.986005e14        # Earth gravitational parameter (m³/s²)
const GPS_F = -4.442807633e-10   # relativistic correction constant (s/√m)

# The one assigned GNSS Type ID in the GGTO messages (CNAV message type 35,
# IS-GPS-705J §20.3.3.8.1 / IS-GPS-200N §30.3.3.8.1; CNAV-2 subframe 3 page 2,
# IS-GPS-800J §3.5.4.2.1): 000 = no data available, 001 = Galileo,
# 010 = GLONASS, 011-111 reserved and to be read as presently unusable. GLONASS
# is assigned but unreachable here: `GNSSSignals` defines no GLONASS
# `TimeSystem` to ask for.
const GPS_GGTO_ID_GALILEO = 1

# `WN_GGTO` is 13 bits, the same width as the CNAV / CNAV-2 week number, so
# resolving it changes nothing except across a rollover. See
# `resolve_reference_week`.
const GPS_GGTO_WN_MODULUS = 8192

"""
    gps_ggto_offset(state::GNSSDecoderState, target::TimeSystem) -> Union{Nothing,GNSSTimeOffset}

Normalise the GGTO parameter set that CNAV (message type 35) and CNAV-2
(subframe 3 page 2) broadcast into a [`GNSSTimeOffset`](@ref), or `nothing`
when the set is absent, marked unavailable, or refers to a system other than
`target` — the GPS counterpart of `beidou_bgto_offset`.

Both messages carry one set at a time, tagged with the system it refers to, so
which target is answerable is the satellite's choice at the moment of asking.
The tag is the `GNSS_ID` field, which both messages' data types name alike.
"""
function gps_ggto_offset(state::GNSSDecoderState, target::TimeSystem)
    data = state.data
    # GNSS Type ID 001 is the only assigned reachable target: 000 is "no data
    # available" and the reserved codes are to be read as unusable rather than
    # guessed at.
    (data.GNSS_ID == GPS_GGTO_ID_GALILEO && target === GST()) || return nothing
    # Decoded together from one message, so all present or all `nothing`. Each is
    # checked anyway, to narrow it for the call below (see `galileo_ggto_offset`).
    (; A_0GGTO, A_1GGTO, A_2GGTO, t_GGTO, WN_GGTO) = data
    (
        isnothing(A_0GGTO) ||
        isnothing(A_1GGTO) ||
        isnothing(A_2GGTO) ||
        isnothing(t_GGTO) ||
        isnothing(WN_GGTO)
    ) && return nothing
    broadcast_time_offset(
        state,
        target,
        A_0GGTO,
        A_1GGTO,
        A_2GGTO,
        t_GGTO,
        WN_GGTO,
        data.WN,
        GPS_GGTO_WN_MODULUS,
    )
end
