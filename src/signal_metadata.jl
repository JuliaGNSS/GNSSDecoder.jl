# ---- Signal metadata --------------------------------------------------------
#
# A decoder state is 1:1 with the signal whose navigation message it
# demodulates, so everything GNSSSignals already knows about that signal — what
# it is called, which constellation and band it belongs to, how fast its
# navigation message runs, which time scale its week numbers and times of week
# are counted in — is knowable from the state alone. Restating any of it here
# would be a second source of truth to keep in sync, so instead each signal file
# states the *one* fact that is genuinely ours (the constants → signal mapping,
# `get_signal_type`) and every accessor below forwards through it.

"""
    get_signal_type(state::GNSSDecoderState) -> Type{<:AbstractGNSSSignal}

Get the GNSSSignals signal type whose navigation message a decoder demodulates,
e.g. `GPSL1CA` for a [`GPSL1CADecoderState`](@ref). This is the single mapping
from this package's decoders back into GNSSSignals; every signal accessor listed
under "Signal metadata" in the docs is forwarded through it, and it is the entry
point for anything not forwarded (`get_signal_type(state)` accepts every
GNSSSignals accessor that takes a signal type).

Returns the *type*, not an instance: constructing a signal builds its spreading
code matrix and SIMD lookup table, which a decoder — working purely in the
symbol domain — never needs. Being a type, it also folds to a compile-time
constant, so the forwarded accessors cost nothing at run time.

Each signal states it for its own decoder, so decoders that share a data
container still answer their own signal: GPS L5-I and L2C-M both decode into a
`GPSCNAVData`, and report `GPSL5I` and `GPSL2CM`.

The mapping names the data-bearing component of a signal pair, since that is
what carries the navigation message: `GPSL2CM` (not the `GPSL2CL` pilot) and
`GalileoE5aI` (not the `GalileoE5aQ` pilot). A decoder built from an
approximation of a signal reports the signal it approximates — the E1B BOC(1,1)
decoder state is a `GalileoE1B`, since the approximation is a
tracking/acquisition concern and the I/NAV stream it decodes is identical.

# Examples

```julia
using GNSSDecoder, GNSSSignals
get_signal_type(GPSL1CADecoderState(1))            # GPSL1CA
get_signal_type(GPSL2CMDecoderState(1))            # GPSL2CM
get_code_length(get_signal_type(GPSL1CADecoderState(1)))  # 1023
```

# See Also

  - [`GNSSDecoderState`](@ref): The state this is queried from
"""
@inline get_signal_type(state::GNSSDecoderState) = get_signal_type(state.constants)

# The GNSSSignals accessors a decoder state can answer for itself, forwarded
# through `get_signal_type`: identity (`get_signal_id` / `get_signal_name`,
# `get_constellation_id` / `get_constellation_name`, `get_band` / `get_band_id` /
# `get_band_name`),
# and the time scale the decoded week numbers and times of week are referenced
# to (`get_time_system` and friends: that scale's epoch as a UTC label,
# `get_system_start_time`, and on the continuous TAI scale,
# `get_tai_system_start_time`, and its offset from TAI, `get_tai_offset`; the
# Galileo epoch is 13 s shy of a UTC midnight *and* anchored to GPS Time rather
# than to UTC, so sourcing them rather than restating them matters).
#
# These extend GNSSSignals' own functions rather than introducing GNSSDecoder
# names, so a caller already doing `using GNSSSignals` needs no new import and
# can query a decoder state exactly as it queries a signal. Accessors that
# describe the *code* (`get_code_length`, `get_chip_rate`,
# `get_carrier_phase_offset`, …) are deliberately not forwarded: they belong to
# acquisition and tracking, not to a symbol-domain decoder, and remain one
# `get_signal_type(state)` away.
for accessor in (
    :get_signal_id,
    :get_signal_name,
    :get_constellation_id,
    :get_constellation_name,
    :get_band,
    :get_band_id,
    :get_band_name,
    :get_time_system,
    :get_time_system_id,
    :get_time_system_name,
    :get_system_start_time,
    :get_tai_system_start_time,
    :get_tai_offset,
)
    @eval @inline GNSSSignals.$accessor(state::GNSSDecoderState) =
        GNSSSignals.$accessor(get_signal_type(state))
end

# The symbol rate is the satellite's: a BeiDou B1I or B3I GEO satellite broadcasts its
# D2 message at 500 Hz, where the signal's rate is D1's 50 Hz, so the PRN goes along.
@inline GNSSSignals.get_bit_rate(state::GNSSDecoderState) =
    GNSSSignals.get_bit_rate(get_signal_type(state), state.prn)

# The TAI epoch a week of the decoder's time system begins at, forwarded like the
# accessors above, so that a decoded week and time of week date an instant:
# `get_week_start_epoch(state, week) + time_of_week`.
@inline GNSSSignals.get_week_start_epoch(state::GNSSDecoderState, week::Integer) =
    GNSSSignals.get_week_start_epoch(get_signal_type(state), week)
