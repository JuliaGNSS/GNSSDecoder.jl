"""
    BroadcastEphemerisValidity

Time references and broadcast fit metadata for one validated ephemeris and
clock set. Both references are absolute `GNSSSignals.TAITime` instants on the
GNSS time scale's TAI timeline. `fit_interval_seconds` is the broadcast curve
fit duration when the message carries enough information to determine it;
`nothing` means the duration is not known here. This metadata does not decide
whether a receiver should accept an ephemeris at a particular epoch.
"""
struct BroadcastEphemerisValidity
    ephemeris_reference::GNSSSignals.TAITime
    clock_reference::GNSSSignals.TAITime
    fit_interval_seconds::Union{Nothing,Int}
end

"""
    get_broadcast_ephemeris_validity(state; approximate_year) -> Union{Nothing,BroadcastEphemerisValidity}

Return the absolute ephemeris and clock reference epochs, plus the broadcast
curve-fit duration when this message makes it knowable. The references are
placed in the week nearest the broadcast stamp (`get_broadcast_time_of_week`),
so a reference just across a week boundary is represented by the adjacent
week's absolute instant. `approximate_year` is used to resolve truncated week
numbers such as GPS LNAV's 10-bit week.

Returns `nothing` until the positioning data set and the week/time stamp are
available, and for messages without their own ephemeris (such as Galileo E6-B).
This reports broadcast metadata only; consumers choose any freshness policy.
"""
function get_broadcast_ephemeris_validity(
    state::GNSSDecoderState;
    approximate_year::Integer,
)
    isnothing(get_broadcast_ephemeris(state)) && return nothing

    week = get_broadcast_week(state; approximate_year)
    stamp_tow = get_broadcast_time_of_week(state)
    (isnothing(week) || isnothing(stamp_tow)) && return nothing

    data = state.data
    week_start = get_week_start_epoch(state, week)
    stamp = week_start + stamp_tow
    toe = _nearest_broadcast_reference(week_start + something(data.t_0e), stamp)
    toc = _nearest_broadcast_reference(week_start + something(data.t_0c), stamp)
    BroadcastEphemerisValidity(toe, toc, _broadcast_fit_interval_seconds(data))
end

@inline _nearest_broadcast_reference(reference, stamp) =
    stamp + GNSSSignals.fold_week_crossover(reference - stamp)

# GPS LNAV is the only currently decoded format here whose curve-fit duration
# is selected by a broadcast fit flag and IODC table. IS-GPS-200N §20.3.3.4.3.1
# and Table 20-XII give 4 hours for flag 0; flag 1 selects 6, 8, 14, or 26
# hours from the IODC ranges. Other formats' nominal transmission cadences do
# not, by themselves, state a curve-fit duration, so they remain unknown here.
@inline function _broadcast_fit_interval_seconds(data::GPSLNAVData)
    flag = data.fit_interval
    iodc = data.IODC
    (isnothing(flag) || isnothing(iodc)) && return nothing
    !flag && return 4 * 3600

    0 <= iodc <= 239 && return 6 * 3600
    256 <= iodc <= 495 && return 6 * 3600
    512 <= iodc <= 751 && return 6 * 3600
    768 <= iodc <= 1007 && return 6 * 3600
    240 <= iodc <= 247 && return 8 * 3600
    (248 <= iodc <= 255 || iodc == 496) && return 14 * 3600
    (497 <= iodc <= 503 || 1021 <= iodc <= 1023) && return 26 * 3600
    nothing
end

@inline _broadcast_fit_interval_seconds(::AbstractGNSSData) = nothing
