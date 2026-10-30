# The broadcast time parameters a navigation file records alongside the ephemeris —
# the UTC polynomial, the leap seconds and Galileo's GGTO — with every truncated
# reference week lifted onto the decoder's own week count, so that no consumer has
# to know how wide the ICD makes each field. The conversions themselves (UTC from
# system time, one system's time from another's) are the consumer's.

"""
$(TYPEDEF)

The broadcast offset of a system's time from UTC, as broadcast: the polynomial
`A_0 + A_1·Δt + A_2·Δt²` (s) at `Δt` seconds from its reference epoch `(WN_0t, t_0t)`, on top of the leap seconds of [`LeapSecondParameters`](@ref).
[`get_utc_parameters`](@ref) returns it.

`WN_0t` is a full week number on the decoder's own count (GPS, Galileo or BeiDou
weeks), not the truncated field the message carries. `A_2` is broadcast by GPS
CNAV and CNAV-2 and the BeiDou B-CNAV messages, and is `0.0` on the others.
BeiDou D1/D2 broadcasts no reference epoch: there `t_0t` and `WN_0t` are
`nothing`, and the polynomial's argument is the BDT time of week itself
(BDS-SIS-ICD-B1I-3.0 §5.2.4.18).

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
    t_0t::Union{Nothing,Int}
    """
    Reference week, on the decoder's own week count; `nothing` on BeiDou D1/D2
    """
    WN_0t::Union{Nothing,Int}
end

"""
$(TYPEDEF)

The broadcast leap seconds of a system's time against UTC: the count now and
after the announced event, and the week and day the event falls at the end of.
[`get_leap_second_parameters`](@ref) returns it.

The counts are the system's own, the leap seconds since its epoch: GPS and
Galileo 18 s since 2017, BeiDou 4 s, 14 fewer. `WN_LSF` is a full week number on
the decoder's own count, not the truncated field the message carries. `DN` is as
broadcast: GPS and Galileo count the day of week from 1 (Sunday) to 7, BeiDou
from 0 to 6.

# Fields

$(TYPEDFIELDS)
"""
struct LeapSecondParameters
    """
    Leap seconds before the event (s)
    """
    Δt_LS::Int
    """
    Leap seconds after the event (s)
    """
    Δt_LSF::Int
    """
    Week the event falls in, on the decoder's own week count
    """
    WN_LSF::Int
    """
    Day of that week at whose end the event takes effect, as broadcast
    """
    DN::Int
end

"""
$(TYPEDEF)

Galileo's broadcast offset of Galileo System Time from GPS time, as broadcast:
the polynomial `A_0G + A_1G·Δt` (s) at `Δt` seconds from its reference epoch
`(WN_0G, t_0G)` (Galileo OS SIS ICD §5.1.8). [`get_ggto_parameters`](@ref)
returns it.

`WN_0G` is a full Galileo week number, not the 6-bit field the message carries.
The epoch is the one [`get_broadcast_time_offset_epoch`](@ref)`(state, GPST())`
dates on TAI, and [`get_broadcast_time_offset`](@ref) evaluates the same
polynomial.

# Fields

$(TYPEDFIELDS)
"""
struct GGTOParameters
    """
    Bias (s)
    """
    A_0G::Float64
    """
    Drift (s/s)
    """
    A_1G::Float64
    """
    Reference time of week (s)
    """
    t_0G::Int
    """
    Reference week, a full Galileo week number
    """
    WN_0G::Int
end

"""
    get_utc_parameters(state::GNSSDecoderState; approximate_year::Integer)
        -> Union{Nothing,UTCParameters}

The broadcast UTC parameters of the decoder's validated data, with the reference
week `WN_0t` resolved, or `nothing` until they and the decoder's week are
decoded.

Every message here truncates `WN_0t`: GPS LNAV and Galileo to 8 bits, GPS CNAV,
CNAV-2 and BeiDou B-CNAV to 13. It is lifted to the week nearest the decoder's own
[`get_broadcast_week`](@ref), on the same count, so
a reference epoch a few days into the next week, while the decoded week still
lags behind, lands in the next week rather than a modulus of weeks back.
`approximate_year` is that of `get_broadcast_week`, which on GPS L1 C/A resolves the
1024-week rollover and is ignored elsewhere.

BeiDou D1/D2 reports its UTC set without a reference epoch (see
[`UTCParameters`](@ref)); Galileo E6-B, which broadcasts none, answers `nothing`.
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

The broadcast leap seconds of the decoder's validated data, with the event's week
`WN_LSF` resolved, or `nothing` until they and the decoder's week are decoded.

`WN_LSF` is truncated like the UTC reference week — 8 bits on GPS LNAV, Galileo
and BeiDou D1/D2, 13 on GPS CNAV, CNAV-2 and BeiDou B-CNAV — and is lifted the
same way, to the week nearest the decoder's own on the same count (see
[`get_utc_parameters`](@ref)). Galileo E6-B, which broadcasts none, answers
`nothing`.
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
    get_ggto_parameters(state::GNSSDecoderState) -> Union{Nothing,GGTOParameters}

Galileo's broadcast GGTO from an I/NAV or F/NAV decoder's validated data, with the
6-bit reference week `WN_0G` resolved to the Galileo week nearest the decoder's own,
or `nothing` until it and the decoder's week are decoded, while the satellite
marks it unavailable, and on every other decoder.

It is the same set [`get_broadcast_time_offset`](@ref) evaluates toward
`GPST()`, resolved the same way. GPS CNAV and CNAV-2 broadcast a GPS-to-GNSS
offset for the system their `GNSS_ID` names instead, which
`get_broadcast_time_offset` evaluates too. No `approximate_year` is needed:
every Galileo week number is a full one.
"""
function get_ggto_parameters(state::GNSSDecoderState{<:AbstractGalileoEphemerisData})
    offset = broadcast_time_offset(state, GPST())
    isnothing(offset) && return nothing
    (; A_0, A_1, t_0, WN_0) = offset
    (isnothing(t_0) || isnothing(WN_0)) && return nothing
    GGTOParameters(A_0, A_1, t_0, WN_0)
end
get_ggto_parameters(
    ::GNSSDecoderState{<:Union{GalileoE6BData,AbstractGPSData,AbstractBeiDouData}},
) = nothing

# Each family names the moduli of its own fields, beside them in its file; the sets
# are built alike from there.
utc_parameters(data::GPSLNAVData, week::Int) =
    utc_parameters(data, week, 0.0, LNAV_UTC_WN_MODULUS)
utc_parameters(data::GPSCNAVData, week::Int) =
    utc_parameters(data, week, data.A_2UTC, CNAV_UTC_WN_MODULUS)
utc_parameters(data::GPSCNAV2Data, week::Int) =
    utc_parameters(data, week, data.A_2UTC, CNAV2_UTC_WN_MODULUS)
utc_parameters(data::AbstractGalileoEphemerisData, week::Int) =
    utc_parameters(data, week, 0.0, GALILEO_UTC_WN_MODULUS)
utc_parameters(data::AbstractBeiDouCNAVData, week::Int) =
    utc_parameters(data, week, data.A_2UTC, BEIDOU_CNAV_UTC_WN_MODULUS)
# D1/D2 broadcasts no reference epoch, so there is no week to resolve.
function utc_parameters(data::BeiDouDNAVData, week::Int)
    (isnothing(data.A_0UTC) || isnothing(data.A_1UTC)) && return nothing
    UTCParameters(something(data.A_0UTC), something(data.A_1UTC), 0.0, nothing, nothing)
end

leap_second_parameters(data::GPSLNAVData, week::Int) =
    leap_second_parameters(data, week, LNAV_LSF_WN_MODULUS)
leap_second_parameters(data::GPSCNAVData, week::Int) =
    leap_second_parameters(data, week, CNAV_LSF_WN_MODULUS)
leap_second_parameters(data::GPSCNAV2Data, week::Int) =
    leap_second_parameters(data, week, CNAV2_LSF_WN_MODULUS)
leap_second_parameters(data::AbstractGalileoEphemerisData, week::Int) =
    leap_second_parameters(data, week, GALILEO_LSF_WN_MODULUS)
leap_second_parameters(data::AbstractBeiDouCNAVData, week::Int) =
    leap_second_parameters(data, week, BEIDOU_CNAV_LSF_WN_MODULUS)
leap_second_parameters(data::BeiDouDNAVData, week::Int) =
    leap_second_parameters(data, week, DNAV_LSF_WN_MODULUS)

# A set with a reference epoch, `A_2` being `0.0` where the message broadcasts none.
# Checked field by field, as `has_klobuchar_coefficients` says why.
function utc_parameters(data, week::Int, A_2, WN_0t_modulus::Int)
    (
        isnothing(data.A_0UTC) ||
        isnothing(data.A_1UTC) ||
        isnothing(A_2) ||
        isnothing(data.t_0t) ||
        isnothing(data.WN_0t)
    ) && return nothing
    UTCParameters(
        something(data.A_0UTC),
        something(data.A_1UTC),
        something(A_2),
        Int(something(data.t_0t)),
        resolve_reference_week(something(data.WN_0t), week, WN_0t_modulus),
    )
end

function leap_second_parameters(data, week::Int, WN_LSF_modulus::Int)
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
        Int(something(data.DN)),
    )
end
