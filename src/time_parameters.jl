# Broadcast time parameters (UTC polynomial, leap seconds, inter-system offsets) with
# every truncated reference week lifted onto the decoder's own week count, so consumers
# need not know the ICD field widths. The conversions themselves (UTC from system time,
# between systems) are the consumer's.

"""
$(TYPEDEF)

The broadcast offset of a system's time from UTC: the polynomial `A_0 + A_1·Δt + A_2·Δt²`
(s) at `Δt` seconds from its reference epoch `(WN_0, t_0)`, on top of the leap seconds of
[`LeapSecondParameters`](@ref). [`get_utc_parameters`](@ref) returns it. The fields are
named as in [`InterSystemTimeOffsetParameters`](@ref), the same polynomial toward another
system's time. The decoders' data containers spell the UTC epoch and the coefficients as
their ICD does: GPS `(WN_t, t_ot)` on LNAV and `(WN_ot, t_ot)` on CNAV and CNAV-2 with
`A_0`, `A_1`, `A_2`; Galileo `(WN_0t, t_0t)` with `A_0`, `A_1`; BeiDou B-CNAV
`(WN_ot, t_ot)` with `A_0UTC`, `A_1UTC`, `A_2UTC`.

`WN_0` is a full week number on the decoder's own count, not the truncated field the
message carries. `A_2` is broadcast by GPS CNAV and CNAV-2 and BeiDou B-CNAV, and is `0.0`
elsewhere. BeiDou D1/D2 broadcasts no reference epoch: `t_0` and `WN_0` are `nothing` and
the polynomial's argument is the BDT time of week itself (BDS-SIS-ICD-B1I-3.0 §5.2.4.18).

# Fields

$(TYPEDFIELDS)
"""
struct UTCParameters
    """
    Bias (s)
    """
    A_0::Float64
    """
    Drift (s/s)
    """
    A_1::Float64
    """
    Drift rate (s/s²); `0.0` where not broadcast
    """
    A_2::Float64
    """
    Reference time of week (s); `nothing` on BeiDou D1/D2
    """
    t_0::Union{Nothing,Int}
    """
    Reference week, on the decoder's own week count; `nothing` on BeiDou D1/D2
    """
    WN_0::Union{Nothing,Int}
end

"""
$(TYPEDEF)

The broadcast leap seconds of a system's time against UTC: the count now and after the
announced event, and the week and day at whose end the event takes effect.
[`get_leap_second_parameters`](@ref) returns it.

The counts are the system's own, accumulated since its epoch: 18 s on GPS and Galileo (since
2017) and 4 s on BeiDou, 14 fewer. `event_week` is a full week number on the decoder's own
count, and `event_day` counts the days of that week from 0 (Sunday) to 6 on every system, so
the event takes effect at `(event_day + 1) · 86400` s of `event_week`. The ICDs write the
fields `Δt_LS`, `Δt_LSF`, `WN_LSF` and `DN`, and count `DN` from 1 on GPS and Galileo and
from 0 on BeiDou.

# Fields

$(TYPEDFIELDS)
"""
struct LeapSecondParameters
    """
    Leap seconds before the event (s)
    """
    leap_seconds::Int
    """
    Leap seconds after the event (s)
    """
    leap_seconds_after_event::Int
    """
    Week the event falls in, on the decoder's own week count
    """
    event_week::Int
    """
    Day of that week at whose end the event takes effect, 0 (Sunday) to 6
    """
    event_day::Int
end

"""
    get_utc_parameters(state::GNSSDecoderState; approximate_year::Integer)
        -> Union{Nothing,UTCParameters}

The broadcast UTC parameters of the decoder's validated data, with the reference week
resolved, or `nothing` until they and the decoder's week are decoded.

Every message here truncates the reference week (the ICDs' `WN_t`, `WN_ot` or `WN_0t`):
GPS LNAV and Galileo to 8 bits, GPS CNAV, CNAV-2 and BeiDou B-CNAV to 13. It is lifted
to the week nearest the decoder's own [`get_broadcast_week`](@ref), so a reference epoch
a few days into the next week, while the decoded week still lags, lands in the next week
rather than a modulus of weeks back. `approximate_year` is passed on to
[`get_broadcast_week`](@ref), where it resolves the 1024-week rollover on GPS L1 C/A; it
is ignored elsewhere.

BeiDou D1/D2 reports its UTC set without a reference epoch (see [`UTCParameters`](@ref));
Galileo E6-B broadcasts none and answers `nothing`.
"""
function get_utc_parameters(state::GNSSDecoderState; approximate_year::Integer)
    data = state.data
    week = get_broadcast_week(data; approximate_year)
    isnothing(week) && return nothing
    utc_parameters(data, week)
end
get_utc_parameters(::GNSSDecoderState{<:GalileoE6BData}; approximate_year::Integer) =
    nothing

"""
    get_leap_second_parameters(state::GNSSDecoderState; approximate_year::Integer)
        -> Union{Nothing,LeapSecondParameters}

The broadcast leap seconds of the decoder's validated data, with the event's week resolved
and its day counted from 0, or `nothing` until they and the decoder's week are decoded.

The event's week, the ICDs' `WN_LSF`, is truncated (8 bits on GPS LNAV, Galileo and BeiDou
D1/D2; 13 on GPS CNAV, CNAV-2 and BeiDou B-CNAV) and lifted like the UTC reference week (see
[`get_utc_parameters`](@ref)). Galileo E6-B broadcasts none and answers `nothing`.
"""
function get_leap_second_parameters(state::GNSSDecoderState; approximate_year::Integer)
    data = state.data
    week = get_broadcast_week(data; approximate_year)
    isnothing(week) && return nothing
    leap_second_parameters(data, week)
end
get_leap_second_parameters(
    ::GNSSDecoderState{<:GalileoE6BData};
    approximate_year::Integer,
) = nothing

"""
    get_inter_system_time_offset_parameters(state::GNSSDecoderState, target::TimeSystem)
        -> Union{Nothing,InterSystemTimeOffsetParameters}

The broadcast offset from the satellite's own time system to `target` (`GPST()`, `GST()` or
`BDT()`), with its reference week resolved to the week nearest the decoder's own. `nothing`
in the cases [`get_broadcast_inter_system_time_offset`](@ref) lists, which evaluates the
same set.

| Message                 | Broadcast                         | Targets           |
|:----------------------- |:--------------------------------- |:----------------- |
| Galileo I/NAV, F/NAV    | GGTO                              | `GPST()`          |
| GPS CNAV, CNAV-2        | GGTO, for the system of `GGTO_ID` | `GST()`           |
| BeiDou B2a, B2b         | BGTO, for the system of `GNSS_ID` | `GPST()`, `GST()` |
| BeiDou B1C              | BGTO, one set kept per system     | `GPST()`, `GST()` |
| BeiDou D1/D2 (B1I, B3I) | without a reference epoch         | `GPST()`, `GST()` |

GPS L1 C/A and Galileo E6-B broadcast none. No `approximate_year` is needed: the decoder's
week on every message that broadcasts an offset needs no rollover resolution.
"""
get_inter_system_time_offset_parameters(state::GNSSDecoderState, target::TimeSystem) =
    inter_system_time_offset_parameters(state, target)

# Each family defines the moduli of its own fields beside them and reads the set under its
# ICD's names: GPS LNAV `A_0`, `A_1`, `t_ot`, `WN_t`; GPS CNAV and CNAV-2 `A_0`, `A_1`,
# `A_2`, `t_ot`, `WN_ot`; Galileo `A_0`, `A_1`, `t_0t`, `WN_0t`; BeiDou B-CNAV `A_0UTC`,
# `A_1UTC`, `A_2UTC`, `t_ot`, `WN_ot`.
utc_parameters(data::GPSLNAVData, week::Int) = utc_parameters(
    data.A_0,
    data.A_1,
    0.0,
    data.t_ot,
    data.WN_t,
    week,
    LNAV_UTC_WN_MODULUS,
)
utc_parameters(data::GPSCNAVData, week::Int) = utc_parameters(
    data.A_0,
    data.A_1,
    data.A_2,
    data.t_ot,
    data.WN_ot,
    week,
    CNAV_UTC_WN_MODULUS,
)
utc_parameters(data::GPSCNAV2Data, week::Int) = utc_parameters(
    data.A_0,
    data.A_1,
    data.A_2,
    data.t_ot,
    data.WN_ot,
    week,
    CNAV2_UTC_WN_MODULUS,
)
utc_parameters(data::AbstractGalileoEphemerisData, week::Int) = utc_parameters(
    data.A_0,
    data.A_1,
    0.0,
    data.t_0t,
    data.WN_0t,
    week,
    GALILEO_UTC_WN_MODULUS,
)
utc_parameters(data::AbstractBeiDouCNAVData, week::Int) = utc_parameters(
    data.A_0UTC,
    data.A_1UTC,
    data.A_2UTC,
    data.t_ot,
    data.WN_ot,
    week,
    BEIDOU_CNAV_UTC_WN_MODULUS,
)
# D1/D2 broadcasts no reference epoch, so there is no week to resolve.
function utc_parameters(data::BeiDouDNAVData, week::Int)
    (isnothing(data.A_0UTC) || isnothing(data.A_1UTC)) && return nothing
    UTCParameters(something(data.A_0UTC), something(data.A_1UTC), 0.0, nothing, nothing)
end

# GPS and Galileo count `DN` from 1 (IS-GPS-200 §20.3.3.5.2.4, Galileo OS SIS ICD §5.1.7),
# BeiDou from 0 (BDS-SIS-ICD-B1I-3.0 §5.2.4.18, B1C-1.0 Table 7-20); Sunday either way.
leap_second_parameters(data::GPSLNAVData, week::Int) =
    leap_second_parameters(data, week, LNAV_LSF_WN_MODULUS, 1)
leap_second_parameters(data::GPSCNAVData, week::Int) =
    leap_second_parameters(data, week, CNAV_LSF_WN_MODULUS, 1)
leap_second_parameters(data::GPSCNAV2Data, week::Int) =
    leap_second_parameters(data, week, CNAV2_LSF_WN_MODULUS, 1)
leap_second_parameters(data::AbstractGalileoEphemerisData, week::Int) =
    leap_second_parameters(data, week, GALILEO_LSF_WN_MODULUS, 1)
leap_second_parameters(data::AbstractBeiDouCNAVData, week::Int) =
    leap_second_parameters(data, week, BEIDOU_CNAV_LSF_WN_MODULUS, 0)
leap_second_parameters(data::BeiDouDNAVData, week::Int) =
    leap_second_parameters(data, week, DNAV_LSF_WN_MODULUS, 0)

# A set with a reference epoch; `A_2` is `0.0` where none is broadcast. The fields arrive
# as positional arguments, each checked in turn, so every `Union{Nothing,T}` narrows
# without boxing (see `has_klobuchar_coefficients`).
function utc_parameters(A_0, A_1, A_2, t_0, WN_0, week::Int, WN_0_modulus::Int)
    (
        isnothing(A_0) ||
        isnothing(A_1) ||
        isnothing(A_2) ||
        isnothing(t_0) ||
        isnothing(WN_0)
    ) && return nothing
    UTCParameters(
        something(A_0),
        something(A_1),
        something(A_2),
        Int(something(t_0)),
        resolve_reference_week(something(WN_0), week, WN_0_modulus),
    )
end

function leap_second_parameters(data, week::Int, WN_LSF_modulus::Int, first_DN::Int)
    (
        isnothing(data.Δt_LS) ||
        isnothing(data.Δt_LSF) ||
        isnothing(data.WN_LSF) ||
        isnothing(data.DN)
    ) && return nothing
    LeapSecondParameters(
        Int(something(data.Δt_LS)),
        Int(something(data.Δt_LSF)),
        resolve_reference_week(something(data.WN_LSF), week, WN_LSF_modulus),
        Int(something(data.DN)) - first_DN,
    )
end
