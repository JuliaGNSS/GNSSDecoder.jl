"""
    get_decoded_time_of_week(state::GNSSDecoderState) -> Union{Nothing,Float64}

Time at the last channel symbol consumed by [`decode!`](@ref), in seconds
relative to the week of [`get_broadcast_time_of_week`](@ref). This advances
the validated broadcast stamp by the decoder's channel-symbol count and rate.
It returns `nothing` when either the stamp or its time anchor is unavailable,
including on Galileo E6-B.

The result is deliberately not reduced modulo one week: it can exceed 604800
when symbols arrive past a week boundary before a new stamp is validated. Keep
it paired with the broadcast week; wrapping it alone would lose the week carry.
Like the broadcast week, it may inherit a temporarily lagging week field.

This is the uncorrected satellite time. A ranging consumer adds the residual
code/carrier phase and applies satellite clock and group-delay corrections.
Vector tracking may adjust the state's time-reference counter to align it with
its code phase. The counter does not report the age of the last preamble.
"""
function get_decoded_time_of_week(state::GNSSDecoderState)
    stamp = get_broadcast_time_of_week(state)
    count = state.num_bits_since_time_reference
    (isnothing(stamp) || isnothing(count)) && return nothing
    Float64(stamp) + count / ustrip(Hz, get_bit_rate(state))
end
