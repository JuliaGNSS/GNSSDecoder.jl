# Definitions shared across GPS signals (counterpart of `galileo/galileo.jl` and
# `beidou/beidou.jl`). Included before the per-signal GPS decoders.

# WGS-84 constants (IS-GPS-200N §20.3.3.4.3; identical in IS-GPS-705J and IS-GPS-800J).
const GPS_μ = 3.986005e14        # Earth gravitational parameter (m³/s²)
const GPS_F = -4.442807633e-10   # relativistic correction constant (s/√m)

# GGTO Type ID (symbol `GGTO ID` since IRN-003; CNAV message type 35, IS-GPS-705J
# §20.3.3.8.1 / IS-GPS-200N §30.3.3.8.1; CNAV-2 subframe 3 page 2, IS-GPS-800J
# §3.5.4.2.1): 000 = no data, 001 = Galileo, 010 = GLONASS, 011-111 reserved
# (unusable). GLONASS is unreachable here: `GNSSSignals` defines no GLONASS `TimeSystem`.
const GPS_GGTO_ID_GALILEO = 1

# `WN_GGTO` is 13 bits, like the CNAV / CNAV-2 week number, so resolving it only
# matters across a rollover. See `resolve_reference_week`.
const GPS_GGTO_WN_MODULUS = 8192

"""
    gps_ggto_offset(state::GNSSDecoderState, target::TimeSystem) -> Union{Nothing,GNSSTimeOffset}

Normalise the GGTO parameter set of CNAV (message type 35) and CNAV-2
(subframe 3 page 2) into a [`GNSSTimeOffset`](@ref), or `nothing` when the set is
absent, marked unavailable, or refers to a system other than `target`. GPS
counterpart of `beidou_bgto_offset`.

Only one set is carried at a time, tagged by `GNSS_ID`, so which target is
answerable is the satellite's choice.
"""
function gps_ggto_offset(state::GNSSDecoderState, target::TimeSystem)
    data = state.data
    # Only ID 001 is reachable; 000 is "no data", reserved codes are unusable.
    (data.GNSS_ID == GPS_GGTO_ID_GALILEO && target === GST()) || return nothing
    # Decoded together, so all present or all `nothing`; each is checked anyway to
    # narrow its type for the call below (see `galileo_ggto_offset`).
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
