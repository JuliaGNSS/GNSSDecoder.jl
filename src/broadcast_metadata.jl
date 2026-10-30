# ---- Satellite and time metadata --------------------------------------------
#
# The questions a positioning engine asks of every decoder, whatever signal it
# is on: which orbit is this satellite in, what week and time of week does its
# decoded message stamp, and how does its system time relate to another
# constellation's.
# Each is answered from broadcast data by some signals and from constellation
# structure by others. Deriving a broadcast quantity from a decoded one is
# decoding work; it belongs here, next to the ICD knowledge, and is stated once
# per signal rather than once per consumer.

"""
    OrbitClass

Orbit class of a navigation satellite: [`get_orbit_class`](@ref) reports it.

The three classes the ICDs of the constellations decoded here distinguish. The
values are deliberately *not* exported, like the health enums' values:
`medium_earth_orbit` is exactly the name a receiver or orbit-propagation
package alongside this one may define, and consumers compare against them
rarely enough that `GNSSDecoder.geostationary_orbit` costs nothing.

# Values

  - `geostationary_orbit`: GEO
  - `inclined_geosynchronous_orbit`: IGSO
  - `medium_earth_orbit`: MEO

The integer values are this package's own; nothing here relies on them matching
any ICD's encoding, and BeiDou's 2-bit `sat_type` (1 GEO / 2 IGSO / 3 MEO) is
mapped across explicitly.
"""
@enum OrbitClass begin
    geostationary_orbit
    inclined_geosynchronous_orbit
    medium_earth_orbit
end

"""
    get_orbit_class(state::GNSSDecoderState) -> Union{Nothing,OrbitClass}

Orbit class of the satellite this decoder is tracking, or `nothing` when the
signal does not make it knowable.

Answered from whichever source the signal actually has:

  - **GPS and Galileo** — [`medium_earth_orbit`](@ref OrbitClass) for every
    satellite. Both constellations are MEO-only by design, so this needs no
    decoded data and is available before the first subframe.
  - **BeiDou B1C, B2a, B2b** — from the broadcast `sat_type` field of the
    satellite's own Ephemeris I block (1 GEO / 2 IGSO / 3 MEO, reserved 0;
    Table 7-6 of each ICD). `nothing` until that block has been decoded and
    promoted, and `nothing` for the reserved code — the same screen
    `is_known_sat_type` applies to positioning.
  - **BeiDou B1I and B3I** — [`geostationary_orbit`](@ref OrbitClass) for the
    GEO PRNs, `nothing` otherwise. D1/D2 NAV broadcasts no orbit-type field, and
    the PRN partition the ICD fixes (BDS-SIS-ICD-B1I-3.0 Table 4-1) separates
    GEO from non-GEO only, so IGSO and MEO cannot be told apart on these
    signals.

!!! warning "`nothing` is not `medium_earth_orbit`"

    Because of that last case, a `nothing` means "this signal cannot say", not
    "not GEO" and not "MEO". Testing `=== GNSSDecoder.geostationary_orbit` is
    exact on every signal — a satellite that is GEO always reports so. Testing
    `=== GNSSDecoder.medium_earth_orbit` is *not*: a B1I MEO satellite reports
    `nothing`. Treat `nothing` as unknown and branch on the GEO case.

# Example

```julia
# Which reference does this BeiDou satellite's ephemeris use?
if get_orbit_class(state) === GNSSDecoder.geostationary_orbit
    # GEO/IGSO semi-major-axis reference
end
```

# See Also

  - [`OrbitClass`](@ref): the returned enumeration
"""
function get_orbit_class end

# GPS and Galileo fly MEO-only constellations, so the answer is a constellation
# fact rather than a decoded one and is stated on the constellation supertype.
# There is deliberately no `AbstractGNSSData` fallback: a constellation with
# GEO or IGSO members (QZSS, NavIC) must state its own answer rather than
# silently inherit "MEO".
get_orbit_class(::GNSSDecoderState{<:AbstractGPSData}) = medium_earth_orbit
get_orbit_class(::GNSSDecoderState{<:AbstractGalileoData}) = medium_earth_orbit

"""
$(TYPEDSIGNATURES)

Seconds of week that the decoder's validated navigation data stamps, or
`nothing` before a time has been decoded.

It is the time of week the message broadcasts for the epoch it stamps, not the
current one, which the section below reads from it.

Every signal here broadcasts one, but four of them do not broadcast it as
seconds: GPS L1C-D counts two-hour intervals plus 18-second TOI steps, BeiDou
B1C counts hours plus 18-second SOH steps, and the rest carry a plain `TOW` or
`SOW`. The reconstruction is the ICD's, so it is done here rather than in every
consumer:

| Signal                                 | Broadcast fields | Seconds of week                                                                               |
|:-------------------------------------- |:---------------- |:--------------------------------------------------------------------------------------------- |
| GPS L1 C/A, GPS CNAV (L2C/L5)          | `TOW`            | `TOW`                                                                                         |
| GPS L1C-D                              | `ITOW`, `toi`    | `ITOW·7200 + toi·18`, `toi = 0` carrying into the next interval (IS-GPS-800 §3.5.2, §3.5.3.2) |
| Galileo I/NAV (E1-B, E5b-I), E5a F/NAV | `TOW`            | `TOW`                                                                                         |
| BeiDou B1I, B3I (D1/D2)                | `SOW`            | `SOW`                                                                                         |
| BeiDou B1C                             | `HOW`, `soh`     | `HOW·3600 + soh·18` (BDS-SIS-ICD-B1C §7.3)                                                    |
| BeiDou B2a, B2b                        | `SOW`            | `SOW`                                                                                         |

Galileo E6-B answers `nothing` always: C/NAV carries HAS corrections stamped
with a time *of hour* per correction block, not a time of week, and those stay
on the blocks that own them (see [`GalileoHASCorrectionBlock`](@ref)).

# The epoch it denotes, and reading the current time

Every signal's value is the time of the epoch that
[`num_bits_since_time_reference`](@ref GNSSDecoderState) is anchored to —
which is why the anchor is set per signal at promotion time, to the frame the
broadcast field actually stamps (the *next* subframe for GPS LNAV, CNAV and
L1C-D; the *current* frame for the BeiDou signals). So the current time of week
is available through [`get_decoded_time_of_week`](@ref):

```julia
using Unitful: Hz, ustrip
time_of_week = get_decoded_time_of_week(state)
```

Because it reads [`data`](@ref GNSSDecoderState) rather than `raw_data`, this
returns a time only once the message set has passed CRC/parity and the
cross-subframe consistency check — the same gate
[`has_positioning_data`](@ref) reports. Call it on
`state.raw_data` directly to see an unvalidated time.

# See Also

  - [`get_broadcast_time_offset`](@ref): relating that time to another constellation's
"""
get_broadcast_time_of_week(state::GNSSDecoderState) = get_broadcast_time_of_week(state.data)

"""
    get_broadcast_week(state::GNSSDecoderState; approximate_year::Integer) -> Union{Nothing,Int}
    get_broadcast_week(data::AbstractGNSSData; approximate_year::Integer) -> Union{Nothing,Int}

Absolute week number of the decoder's validated navigation data, counted from
the origin of its own time system (`get_system_start_time`), or `nothing`
before a week number has been decoded.

It is the week the message broadcasts, not the current one: across a week
crossover it lags behind (see below). GNSSSignals' `get_week_start_epoch`
answers for a decoder state, so the broadcast week and time of week date an
instant on TAI:

```julia
week = get_broadcast_week(state; approximate_year = 2026)
time_of_week = get_broadcast_time_of_week(state)
if !isnothing(week) && !isnothing(time_of_week)
    epoch = get_week_start_epoch(state, week) + time_of_week  # a GNSSSignals.TAITime
end
```

Every message here broadcasts its week number wide enough to be absolute for
decades (Galileo I/NAV and F/NAV 12 bits; GPS CNAV, CNAV-2, BeiDou D1/D2 and
B-CNAV 13 bits), except GPS L1 C/A: LNAV sends only the week modulo 1024
(IS-GPS-200 §20.3.3.3.1.1), which rolls over every 19.6 years (1999-08-22,
2019-04-07, 2038-11-21, …). That rollover is resolved here, against
`approximate_year`: the 1024-week cycle that puts the broadcast week nearest
the middle of that year. Any year within ±9 years of the observation selects
the right cycle: the current year for a live signal, the year a recording was
made in for a replay.

`approximate_year` is required on every signal, though only GPS L1 C/A reads
it, so that a call reads the same whatever decoder it is handed. There is no
default: the decoder has no notion of the current time, and the year an
observation was made in is the caller's to know.

Galileo E6-B answers `nothing` always: C/NAV broadcasts no week number (see
[`get_broadcast_time_of_week`](@ref)).

The week is the one the decoded `WN` field states, which the message that
carries it updates on its own cadence. Across a week crossover it can lag the
time of week by a subframe, as the signals' `reset_decoder_state!` notes say.

Because it reads [`data`](@ref GNSSDecoderState) rather than `raw_data`, the
state form answers only once the message set has been validated; call it on
`state.raw_data` to see an unvalidated week.

# See Also

  - [`get_broadcast_time_of_week`](@ref):
    the time within that week
"""
get_broadcast_week(state::GNSSDecoderState; approximate_year::Integer) =
    get_broadcast_week(state.data; approximate_year)

# Every signal states its own `get_broadcast_week(data; approximate_year)`, beside its
# `get_broadcast_time_of_week`, and there is no fallback reading `data.WN`: a future signal
# with a truncated week number must resolve it rather than silently inherit the
# truncated value.

"""
    get_broadcast_time_offset(state::GNSSDecoderState, target::TimeSystem, time_of_week) -> Union{Nothing,Float64}

The broadcast offset from this satellite's own time scale to `target` (s), at
`time_of_week` (s) on the satellite's own scale, or `nothing` when this decoder
has no such offset.

`target` is a `GNSSSignals.TimeSystem` — `GPST()`, `GST()` or `BDT()`. The
satellite's own scale is `get_time_system(state)`; asking for that returns
`nothing`, because the offset from a scale to itself is not a broadcast
quantity.

The offset `Δt` is **the broadcasting system's time minus the target system's**
— `t_own - t_target`, the direction every ICD here defines the quantity in
(Galileo OS SIS ICD §5.1.8; IS-GPS-800/200 GGTO; BDS-SIS-ICD-B1C §7.13.2 Eq.
7-30). To convert a seconds-of-week reading, *subtract* it, and add the leap
seconds between the two scales, GNSSSignals'
`get_time_scale_offset(get_time_system(state), target)`, which the note below
explains: `t_target = t_own - Δt + leap_seconds`.

It is the broadcast polynomial

```julia
Δτ = time_of_week - t_0 + 604800 * (WN - WN_0)
Δt = A_0 + A_1 * Δτ + A_2 * Δτ^2
```

evaluated against the decoder's own week number `WN`. The reference week `WN_0`
is sent truncated, and not always to the width of the week it is subtracted
from — Galileo's `WN_0G` is 6 bits against a 12-bit `WN` — so it is lifted into
the decoder's own week numbering first; pairing the two here, the caller needs
only the time of week. The week number is the decoded `WN`, so right after a
week crossover, until the message carrying it is decoded again, it can still be
the previous week's. Galileo's GGTO has no quadratic term. BeiDou D1/D2's
two-term offsets carry no reference epoch: their argument is the time of week
itself (BDS-SIS-ICD-B1I-3.0 §5.2.4.19, `Δt_GPS = A_0GPS + A_1GPS · t_E`).

!!! note "The steering residual alone"

    Two systems differ by the leap seconds between them plus a broadcast
    steering residual, and only the second of those is on the air. The bias
    coefficient these messages carry is 16 bits at 2⁻³⁵ s — a range of ±0.95 µs —
    or on BeiDou D1/D2 14 bits at 0.1 ns, ±0.82 µs, so it cannot express a
    whole-second offset even in principle, and does not try to: real BeiDou
    satellites broadcast tens of nanoseconds against a BDT-to-GPST offset of 14
    seconds. IS-GPS-800J
    §3.5.4.2.1.1 states the split outright: the GGTO parameters provide the
    offset "modulo one second", and "users must also apply any integer seconds
    difference between the systems using definitions of each system time scale".

    `Δt` is evaluated from the coefficients as broadcast, as RINEX navigation
    files and other receivers report them, so it is the steering residual alone.
    The leap seconds follow from the scales' definitions:
    every one of these scales is a fixed offset from TAI, set to TAI − UTC when
    it was defined and never adjusted since (19 s for GPST, which GST adopted,
    and 33 s for BDT), so two of them differ by the leap seconds UTC gained in
    between. That is zero between GPST and GST and +14 s from BDT toward either
    (a BDT reading is 14 s *behind*), so a conversion that is only ever tried on
    Galileo never notices them missing.

!!! warning "Seconds of week, not week numbers"

    This converts a time *of week*; it says nothing about which week. The scales
    do not share a week origin — BDT week 0 begins 14 s into GPS week 1356 — so
    a full `(WN, TOW)` conversion goes through TAI, which carries the leap
    seconds, the epoch offset and the week over in one step:

    ```julia
    epoch = get_week_start_epoch(BDT(), WN_bdt) + t_bdt - Δt  # a GNSSSignals.TAITime
    WN_gps, t_gps = get_week(GPST(), epoch), get_time_of_week(GPST(), epoch)
    ```

    GNSSSignals takes the whole weeks off in integer seconds, so the nanoseconds
    of `Δt` survive the round trip.

`nothing` covers, and does not distinguish between, the three ways a signal can
fail to have an offset:

  - it broadcasts no inter-system offset at all (GPS L1 C/A, Galileo E6-B),
  - it does, but has not decoded one yet,
  - it has, but for a different system than `target` — the GPS, B2a and B2b
    messages carry a single offset set tagged with the system it refers to, so
    which one is on the air is the satellite's choice, while Galileo's is
    GPS-only by definition and BeiDou B1C keys a set per system.

The ICDs' "not available" sentinels are screened here, so a returned offset is
always a real one: a GPS/BeiDou `GNSS_ID`/`GGTO_ID` of 0 means "no offset in
this message" (IS-GPS-800 Table 6.2-7, BDS-SIS-ICD-B1C §7.13.1), and Galileo
marks an absent GGTO by setting all four fields to all ones (§5.1.8, handled at
decode by `galileo_ggto`).

The coefficients as broadcast stay on the decoded data (`data.A_0G`,
`data.A_0BGTO`, `data.A_0GGTO`, their reference weeks `data.WN_0G`,
`data.WN_0BGTO`, `data.WN_GGTO`, …), where a RINEX writer reads them. GLONASS
offsets are decoded and published there too, but cannot be asked for here:
`GNSSSignals` defines no GLONASS `TimeSystem`, this package decodes no GLONASS
signal, and inventing a target type for a constellation neither package handles
would be a name with nothing behind it.

# Example

```julia
using Unitful: s, ustrip
# Steer a BeiDou observable onto GPS time.
Δt = get_broadcast_time_offset(state, GPST(), t_bdt)
if !isnothing(Δt)
    # The leap seconds between the scales are not broadcast: +14 s from BDT to GPST.
    t_gps = t_bdt - Δt + ustrip(s, get_time_scale_offset(BDT(), GPST()))
end
```

# See Also

  - [`get_broadcast_time_of_week`](@ref):
    the time of week the offset is evaluated at
  - [`get_broadcast_time_offset_epoch`](@ref): the reference epoch it is evaluated
    from, to judge its age by
"""
function get_broadcast_time_offset(
    state::GNSSDecoderState,
    target::TimeSystem,
    time_of_week::Real,
)
    offset = broadcast_time_offset(state, target)
    isnothing(offset) && return nothing
    (; A_0, A_1, A_2, t_0, WN_0) = offset
    # A resolved `WN_0` implies a decoded `WN`: the record is only built with one.
    Δτ =
        isnothing(t_0) || isnothing(WN_0) ? Float64(time_of_week) :
        time_of_week - t_0 + SECONDS_PER_WEEK * (something(state.data.WN) - WN_0)
    return A_0 + A_1 * Δτ + A_2 * Δτ^2
end

"""
    get_broadcast_time_offset_epoch(state::GNSSDecoderState, target::TimeSystem) -> Union{Nothing,TAITime}

The reference epoch of the broadcast offset to `target` that
[`get_broadcast_time_offset`](@ref) evaluates: its reference week and time of week,
dated on TAI as GNSSSignals' `get_week_start_epoch` dates the satellite's own scale,
or `nothing` when there is no such offset or it carries no reference epoch.

The broadcast polynomial is extrapolated from this epoch, and a decoder keeps an
offset set until the satellite broadcasts a newer one — BeiDou B1C keeps one per
system, which a satellite may stop refreshing. How far from it an offset is still
good enough is the caller's to decide, not the decoder's:

```julia
epoch = get_week_start_epoch(BDT(), WN_bdt) + t_bdt
reference = get_broadcast_time_offset_epoch(state, GPST())
age = isnothing(reference) ? nothing : epoch - reference  # s
```

BeiDou D1/D2 answers `nothing` here even where it has an offset: its two-term
offsets carry no reference epoch, and their argument is the time of week itself.
"""
function get_broadcast_time_offset_epoch(state::GNSSDecoderState, target::TimeSystem)
    offset = broadcast_time_offset(state, target)
    isnothing(offset) && return nothing
    (; t_0, WN_0) = offset
    (isnothing(t_0) || isnothing(WN_0)) && return nothing
    return GNSSSignals.get_week_start_epoch(offset.source, WN_0) + t_0
end

"""
    BroadcastTimeOffsetParameters

Broadcast steering-residual polynomial from `source` to `target`, as returned
by [`get_broadcast_time_offset_parameters`](@ref). The offset is
`source_time - target_time`; the integer seconds between the time scales are
not included.

`A_0`, `A_1`, `A_2` are bias (s), drift (s/s), and quadratic coefficient
(s/s²). `t_0` is the reference seconds of week and `WN_0` its week on the
source system's count, resolved against the decoder's week. Together they date
an absolute reference with `get_week_start_epoch(parameters.source, WN_0) + t_0`.
Both are `nothing` on BeiDou D1/D2, whose polynomial argument is the source
seconds of week itself. Terms absent from the message are zero.

This is a snapshot of the coefficients; it contains no decoder storage.
Original ICD spellings remain on the decoded message containers. Galileo's
[`GGTOParameters`](@ref) is also available for consumers of that format.
"""
struct BroadcastTimeOffsetParameters{S<:TimeSystem,T<:TimeSystem}
    source::S
    target::T
    A_0::Float64
    A_1::Float64
    A_2::Float64
    t_0::Union{Nothing,Int}
    WN_0::Union{Nothing,Int}
end

"""
    get_broadcast_time_offset_parameters(state::GNSSDecoderState, target::TimeSystem)
        -> Union{Nothing,BroadcastTimeOffsetParameters}

Read a broadcast inter-system time polynomial with explicit source and target
scales, normalized coefficient names and a resolved reference week. Returns
`nothing` for missing, unavailable, or differently targeted offsets, using the
same checks as [`get_broadcast_time_offset`](@ref).

Use this to store or write coefficients; use `get_broadcast_time_offset` to
evaluate the polynomial at a time of week, and
[`get_broadcast_time_offset_epoch`](@ref) to date its reference.
"""
get_broadcast_time_offset_parameters(state::GNSSDecoderState, target::TimeSystem) =
    broadcast_time_offset(state, target)

# The broadcast offset from the satellite's own scale to `target` as a
# `BroadcastTimeOffsetParameters`, or `nothing`, in the cases `get_broadcast_time_offset` lists.
# Every signal states its own method, in its own file. There is no blanket `nothing`
# fallback for the same reason `get_orbit_class` has none: a future signal that does
# broadcast an offset would silently inherit "no offset available" and its
# correction would go quietly unapplied.
function broadcast_time_offset end

"""
    assemble_broadcast_time_offset(source, target, A_0, A_1, A_2, t_0, WN_0, WN,
                                   WN_0_modulus) -> Union{Nothing,BroadcastTimeOffsetParameters}

Assemble a `BroadcastTimeOffsetParameters` from one message's broadcast
coefficients, resolving the truncated reference week `WN_0` against the
decoder's own week `WN` (see `resolve_reference_week`). Every
`get_broadcast_time_offset` method routes through this rather than calling the
constructor, so that resolution has one place to be right instead of six.

The arguments are positional on purpose. `t_0`, `WN_0` and `WN` are
`Union{Nothing,Int}` fields at every call site, and as keyword arguments they
would be collected into a `NamedTuple` of that non-concrete type, which is
heap-allocated on every call — and `get_broadcast_time_offset` is read once per
satellite per navigation epoch, on a path that must stay allocation-free.
"""
function assemble_broadcast_time_offset(
    source::TimeSystem,
    target::TimeSystem,
    A_0::Real,
    A_1::Real,
    A_2::Real,
    t_0::Union{Nothing,Integer},
    WN_0::Union{Nothing,Integer},
    WN::Union{Nothing,Integer},
    WN_0_modulus::Int,
)
    reference_week = resolve_reference_week(WN_0, WN, WN_0_modulus)
    # A reference epoch that cannot be placed in the decoder's own week
    # numbering cannot be evaluated against, and `Δτ` needs `WN` anyway, so
    # there is no usable offset yet rather than one with an unresolved epoch.
    isnothing(reference_week) && !isnothing(WN_0) && return nothing
    BroadcastTimeOffsetParameters(
        source,
        target,
        A_0,
        A_1,
        A_2,
        isnothing(t_0) ? nothing : Int(t_0),
        reference_week,
    )
end

"""
    resolve_reference_week(WN_0, WN, modulus) -> Union{Nothing,Int}

Lift a truncated broadcast reference week `WN_0` into the same numbering as the
decoder's own week `WN`, by picking the representative nearest to `WN`.

Every one of these messages sends the reference week truncated, and not all of
them truncate it to the width of the week number it will be subtracted from:

| field             | width   | against a `WN` of |
|:----------------- |:------- |:----------------- |
| Galileo `WN_0G`   | 6 bits  | 12 bits           |
| GPS `WN_GGTO`     | 13 bits | 13 bits           |
| BeiDou `WN_0BGTO` | 13 bits | 13 bits           |

So `WN - WN_0` as broadcast is meaningful on GPS and BeiDou and meaningless on
Galileo, where it subtracts a 0..63 count from a 0..4095 one. Left unresolved it
does not fail visibly: it inflates `Δτ` by up to 4032 weeks, which a full-scale
`A_1G` turns into 0.7 ms — 200 km of range — on a quantity whose true value is
tens of nanoseconds.

Resolving against the *nearest* representative rather than the most recent one
matters at a week boundary: a reference epoch a few hours ahead of the current
week is legal, and rounding it down would place it `modulus - 1` weeks in the
past instead. On the two 13-bit fields the result is `WN_0` unchanged except
across a week-number rollover, where it is what keeps the difference small
instead of hugely negative.

`nothing` when either week is unknown — the reference cannot be placed without
the current week, and nothing can be evaluated without it either.
"""
function resolve_reference_week(
    WN_0::Union{Nothing,Integer},
    WN::Union{Nothing,Integer},
    modulus::Int,
)
    (isnothing(WN_0) || isnothing(WN)) && return nothing
    half = modulus ÷ 2
    return Int(WN) - (mod(Int(WN) - Int(WN_0) + half, modulus) - half)
end

"""
    get_group_delay(state::GNSSDecoderState, signal::AbstractGNSSSignal) -> Union{Nothing,Float64}

Group delay (s) of the ranging `signal` against the reference the broadcast
satellite clock is stated for, from the decoder's validated data: subtract it
from the clock polynomial's offset `Δt_SV` to get the clock offset on `signal`,
`Δt_SV(signal) = Δt_SV - get_group_delay(state, signal)`. `nothing` until
[`has_positioning_data`](@ref), and until every term of the
signal's correction has been decoded (see below).

!!! note "Combining signals"

    This is each signal's own correction, the ICDs' single-frequency one. A
    combination of ranges, such as the ionosphere-free one, takes the same
    combination of their group delays: with `γ = (f₁/f₂)²`, the ionosphere-free
    pair's is `(g₂ - γ·g₁) / (1 - γ)` for the group delays `g₁`, `g₂` this returns,
    which is how the ICDs state their dual-frequency corrections
    (BDS-SIS-ICD-B2a-1.0 Eqs. 7-18 to 7-21), and which vanishes for the pair a
    clock is stated for (GPS L1/L2 P(Y), Galileo E1/E5b on I/NAV and E1/E5a on
    F/NAV).

`signal` is the GNSSSignals signal instance the range was measured on, as the
tracking loop holds it. It may be a dataless pilot whose navigation data this decoder takes from
the band's data component (GPS L5-Q and L1C-P, Galileo E1-C, E5a-Q and E5b-Q,
BeiDou B1C and B2a pilots), and the two need not share a correction, so the
signal is asked for rather than read off the decoder:

| Message                   | Ranging signal       | Group delay                                   |
|:------------------------- |:-------------------- |:--------------------------------------------- |
| GPS LNAV (L1 C/A)         | L1 C/A               | `T_GD`                                        |
| GPS CNAV (L5, L2C)        | L5-I, L5-Q, L2 CM/CL | `T_GD - ISC_L5I5`, `- ISC_L5Q5`, `- ISC_L2C`  |
| GPS CNAV-2 (L1C)          | L1C-D, L1C-P, L1 C/A | `T_GD - ISC_L1CD`, `- ISC_L1CP`, `- ISC_L1CA` |
| Galileo I/NAV (E1-B, E5b) | E1-B/C, E5b-I/Q      | `(f_E1 / f)² · BGD_E1_E5b`                    |
| Galileo F/NAV (E5a)       | E5a-I/Q              | `(f_E1 / f)² · BGD_E1_E5a`                    |
| BeiDou D1/D2 (B1I, B3I)   | B1I, B3I             | `T_GD1`, `0`                                  |
| BeiDou B-CNAV1 (B1C)      | B1C data, pilot      | `T_GD_B1Cp + ISC_B1Cd`, `T_GD_B1Cp`           |
| BeiDou B-CNAV2 (B2a)      | B2a data, pilot      | `T_GD_B2ap + ISC_B2ad`, `T_GD_B2ap`           |
| BeiDou B-CNAV3 (B2b)      | B2b-I                | `T_GD_B2bI`                                   |

GPS states its clock for the L1 P(Y)/L2 P(Y) ionosphere-free pair, with an
inter-signal correction per civil signal (IS-GPS-200 §20.3.3.3.3.2 and
§30.3.3.3.1.1, IS-GPS-705 §20.3.3.3.1.2, IS-GPS-800 §3.5.3.9.1). Galileo states
it for an ionosphere-free pair, E1/E5b for I/NAV and E1/E5a for F/NAV, and a
single-frequency user on the pair's second band scales the BGD by the squared
frequency ratio, 1.79 on E5a and 1.70 on E5b, exactly 1 on E1 (OS SIS ICD
§5.1.5). BeiDou states it for B3I, which so needs none, and publishes a group
delay per signal from there, with the data components adding an ISC to their
band's pilot (BDS-SIS-ICD-B1I-3.0 §5.2.4.10, BDS-SIS-ICD-B1C-1.0 §7.6 and the
B2a/B2b counterparts).

Most terms ride with the clock and are there once the gate passes, CNAV-2's
`T_GD`, `ISC_L1CD` and `ISC_L1CP` among them (subframe 2, with the ephemeris).
GPS CNAV's `T_GD` and ISCs (message type 30), CNAV-2's ISCs of the other civil
signals (`ISC_L1CA` and the L2/L5 ones, subframe 3 page 1) and B2a's `T_GD_B2ap`
and `ISC_B2ad` (message type 30) arrive on their own cadence, so the gate does
not wait for them (see [`has_positioning_data`](@ref)), and
this returns `nothing` until they have. What to do meanwhile is the caller's choice: a single-frequency
fix that takes the correction as zero, `something(get_group_delay(state, signal), 0.0)`,
is off by about a metre, and an ionosphere-free combination scales that by the
missing range's coefficient (1.3 to 2.6 for the GPS and BeiDou pairs), so it may
rather wait.

Each message has methods for exactly the ranging signals it broadcasts a
correction for, the rows of the table above, whichever of its signals the
decoder is on: an I/NAV decoder on E1-B answers for E5b ranges too, and a
CNAV-2 decoder for L1 C/A. Any other pairing, such as an E5a range (corrected by
F/NAV's BGD) on an I/NAV decoder, or any range on Galileo E6-B, which carries no
group delay, throws a `MethodError` rather than returning a plausible number.
"""
function get_group_delay end
