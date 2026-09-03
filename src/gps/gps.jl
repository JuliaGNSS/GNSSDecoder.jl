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
    gps_ggto_offset(state::GNSSDecoderState, target::TimeSystem, ggto_id) -> Union{Nothing,GNSSTimeOffset}

Normalise the GGTO parameter set that CNAV (message type 35) and CNAV-2
(subframe 3 page 2) broadcast into a [`GNSSTimeOffset`](@ref), or `nothing`
when the set is absent, marked unavailable, or refers to a system other than
`target` — the GPS counterpart of `beidou_bgto_offset`.

Both messages carry one set at a time, tagged with the system it refers to, so
which target is answerable is the satellite's choice at the moment of asking.
`ggto_id` is that broadcast tag, passed in because the two messages name the
same field differently (`GNSS_ID` on CNAV, `GGTO_ID` on CNAV-2).
"""
function gps_ggto_offset(state::GNSSDecoderState, target::TimeSystem, ggto_id)
    data = state.data
    # GNSS Type ID 001 is the only assigned reachable target: 000 is "no data
    # available" and the reserved codes are to be read as unusable rather than
    # guessed at.
    (ggto_id == GPS_GGTO_ID_GALILEO && target === GST()) || return nothing
    isnothing(data.A_0GGTO) && return nothing
    broadcast_time_offset(
        state,
        target,
        data.A_0GGTO,
        data.A_1GGTO,
        data.A_2GGTO;
        t_0 = data.t_GGTO,
        WN_0 = data.WN_GGTO,
        WN = data.WN,
        WN_0_modulus = GPS_GGTO_WN_MODULUS,
    )
end
