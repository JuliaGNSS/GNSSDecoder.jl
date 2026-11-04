# Definitions shared across BeiDou signals (B1I/B3I D1-D2 NAV and the B-CNAV family on
# B1C, B2a, and B2b): the BDCS constants (all five decoders) and the midi and reduced
# almanac records (the three B-CNAV decoders). Per-signal framing, layout, and parsing
# live in `dnav.jl` + `b1i.jl`/`b3i.jl`, `b1c.jl`, `b2a.jl`, `b2b.jl`.
# Something belongs here only if no single signal file can own it: the GEO-orbit PRN
# partition selects the D1/D2 message format, so it lives in `dnav.jl`.

# BDCS (BeiDou Coordinate System) constants (BDS-SIS-ICD-B1I-3.0 Table 3-1 / §5.2.4.12;
# identical across the B1I, B3I, B1C, B2a, and B2b ICDs). Only the ones differing from
# the GPS WGS-84 values live here; π is `GNSS_PI` in `gnss.jl`. The Earth rotation rate
# differs from GPS/Galileo (7.2921151467e-5 rad/s), so the package-wide
# `EARTH_ROTATION_RATE` must *not* be used in BeiDou's ephemeris.
const BEIDOU_μ = 3.986004418e14                # geocentric gravitational constant (m³/s²)
const BEIDOU_EARTH_ROTATION_RATE = 7.2921150e-5 # Earth rotation rate of BDCS (rad/s)
# Relativistic correction constant F = −2√μ/c² (s/√m; B1C-1.0 §7.5.2 and the
# B2a/B2b counterparts): Galileo's value, because μ agrees.
const BEIDOU_F = -4.442807309e-10

# Whether a broadcast B-CNAV satellite orbit type names an orbit class: 1 = GEO,
# 2 = IGSO, 3 = MEO, 0 reserved (BDS-SIS-ICD-B1C-1.0 and B2a-1.0 Table 7-8, B2b-1.0
# Table 7-6).
#
# This gates positioning: `sat_type` selects the reference of the broadcast ephemeris,
# `A = A_ref + ΔA` with `A_ref = 27 906 100 m` (MEO) or `42 162 200 m` (IGSO/GEO), and
# `i = i_ref + δi` in the almanacs. The references differ by 14 256 100 m, so a
# satellite broadcasting the reserved code must not be used. B1I/B3I are unaffected:
# D1/D2 NAV broadcasts `sqrt_A` outright and has no orbit-type field.
#
# The almanac records keep the raw field; a reserved type in another satellite's almanac
# is no reason to discard the rest of the page.
is_known_sat_type(sat_type::Integer) = sat_type in 1:3
is_known_sat_type(::Nothing) = false

"""
    BeiDouMidiAlmanac

Midi almanac for one BeiDou satellite.

The BDS-3 midi almanac is one 156-bit block, identical across the three B-CNAV
messages, so one record serves the B1C decoder (B-CNAV1 subframe 3 page type 4,
BDS-SIS-ICD-B1C-1.0 Figure 6-21 / Table 7-13), the B2a decoder (B-CNAV2 message
type 40, BDS-SIS-ICD-B2a-1.0 Figure 6-20 / Table 7-13), and the B2b decoder
(B-CNAV3 message type 40, BDS-SIS-ICD-B2b-1.0 Figure 6-15 / Table 7-11), as
[`GalileoAlmanac`](@ref) does for I/NAV and F/NAV. Each record is keyed by `PRN_a`
in the per-signal data containers.

Reference value (ICD notes to the midi-almanac tables): `i = i_ref + δi` with
`i_ref = 0.30π rad` (54°) for MEO/IGSO and `0.0` for GEO, selected by `sat_type`.
Semi-circle fields are converted to radians on decode.

# Fields

  - `PRN_a::Int`: Almanac satellite PRN (1-63; a zero PRN marks an empty
    block and is never stored).
  - `sat_type::Int`: Satellite orbit type (2 bits: 1 = GEO, 2 = IGSO,
    3 = MEO, 0 reserved).
  - `WN_a::Int`: Almanac reference week number (BDT week).
  - `t_0a::Int`: Almanac reference time of week (seconds, LSB 2¹²).
  - `e::Float64`: Eccentricity (dimensionless).
  - `δi::Float64`: Inclination delta from the `sat_type` reference (rad).
  - `sqrt_A::Float64`: Square root of the semi-major axis (√m).
  - `Ω_0::Float64`: Longitude of ascending node at weekly epoch (rad).
  - `Ω_dot::Float64`: Rate of right ascension (rad/s).
  - `ω::Float64`: Argument of perigee (rad).
  - `M_0::Float64`: Mean anomaly at reference time (rad).
  - `a_f0::Float64`, `a_f1::Float64`: Satellite clock bias / drift (s, s/s).
    The ICDs spell these `af0`/`af1` (B2b-1.0 Table 7-11, B1C/B2a-1.0 Table 7-13),
    unlike the clock-correction block's `a_0`/`a_1`/`a_2`.
  - `health::Int`: Raw 8-bit satellite health word (BDS-SIS-ICD-B2b-1.0
    Table 7-12, the fullest definition: bit 8 (MSB) = satellite clock,
    bit 7 = B1C signal, bit 6 = B2a signal, bit 5 = B2b_I signal, bits 4-1
    reserved; 0 = healthy. The earlier B1C/B2a ICDs' Table 7-14 defines the
    same word without the B2b bit).
"""
Base.@kwdef struct BeiDouMidiAlmanac
    PRN_a::Int
    sat_type::Int
    WN_a::Int
    t_0a::Int
    e::Float64
    δi::Float64
    sqrt_A::Float64
    Ω_0::Float64
    Ω_dot::Float64
    ω::Float64
    M_0::Float64
    a_f0::Float64
    a_f1::Float64
    health::Int
end

"""
    BeiDouReducedAlmanac

Reduced almanac for one BeiDou satellite.

The BDS-3 reduced almanac is one 38-bit block, identical across the three B-CNAV
messages, so one record serves the B1C decoder (B-CNAV1 subframe 3 page type 2, four
blocks per page, BDS-SIS-ICD-B1C-1.0 Figure 6-18 / Table 7-16), the B2a decoder
(B-CNAV2 message types 31 and 33, BDS-SIS-ICD-B2a-1.0 Figure 6-17 / Table 7-16), and
the B2b decoder (B-CNAV3 message type 40, five blocks per message, BDS-SIS-ICD-B2b-1.0
Figure 6-12 / Table 7-14). The block carries no epoch; the almanac reference week/time
of the carrying page/message is copied into every record. Each record is keyed by
`PRN_a` in the per-signal data containers.

Reference values to apply (ICD notes to the reduced-almanac tables):
`A = A_ref + δA` with `A_ref = 27 906 100 m` (MEO) or `42 162 200 m`
(IGSO/GEO) selected by `sat_type`; `Φ₀ = M₀ + ω` relative to `e = 0` and
`δi = 0` with `i = 55°` (MEO/IGSO) or `i = 0°` (GEO). The user algorithm is
the midi almanac's with the missing parameters set to zero.

# Fields

  - `PRN_a::Int`: Almanac satellite PRN (1-63; a zero PRN marks an empty
    block and is never stored).
  - `sat_type::Int`: Satellite orbit type (2 bits: 1 = GEO, 2 = IGSO,
    3 = MEO, 0 reserved).
  - `WN_a::Int`: Almanac reference week number (from the carrying page/message).
  - `t_0a::Int`: Almanac reference time of week (seconds, from the carrying page/message).
  - `δA::Float64`: Semi-major-axis correction to the `sat_type` reference (m).
  - `Ω_0::Float64`: Longitude of ascending node at weekly epoch (rad).
  - `Φ_0::Float64`: Argument of latitude at reference time, `M₀ + ω` (rad).
  - `health::Int`: Raw 8-bit satellite health word (same layout as
    [`BeiDouMidiAlmanac`](@ref)'s).
"""
Base.@kwdef struct BeiDouReducedAlmanac
    PRN_a::Int
    sat_type::Int
    WN_a::Int
    t_0a::Int
    δA::Float64
    Ω_0::Float64
    Φ_0::Float64
    health::Int
end

# ---- Shared B-CNAV data blocks ----------------------------------------------
#
# B1C, B2a and B2b build their messages from the same eight data blocks, each defined
# once per ICD (layout, widths, signedness, scale factors) and placed at a
# message-specific bit offset. Parsing each block once keeps the decoders from
# drifting (they had diverged in scale-factor spelling and integer types).
#
# Each function takes the packed message word, its logical bit length, and the 1-based
# position of the block's first bit, and returns a `NamedTuple` keyed by the container
# field names, which the three containers share. Callers splat it into their container:
#
#     BeiDouCNAV3Data(raw; beidou_clock_block(word, word_length, 44)...)
#
# B1C figure/table numbers are cited; the B2a-1.0 and B2b-1.0 counterparts are named at
# each call site.

"""
Clock correction block, 69 bits: `t_0c`(11, LSB 300 s), `a_f0`(25, 2⁻³⁴ s),
`a_f1`(22, 2⁻⁵⁰ s/s), `a_f2`(11, 2⁻⁶⁶ s/s²), the ICD's `t_oc`, `a_0`, `a_1`, `a_2` —
BDS-SIS-ICD-B1C-1.0 Figure 6-13, Table 7-5.
"""
beidou_clock_block(word::Unsigned, word_length::Int, start::Int) = (;
    t_0c = Int64(get_bits(word, word_length, start, 11)) * 300,
    a_f0 = get_twos_complement_num(word, word_length, start + 11, 25) * 2.0^-34,
    a_f1 = get_twos_complement_num(word, word_length, start + 36, 22) * 2.0^-50,
    a_f2 = get_twos_complement_num(word, word_length, start + 58, 11) * 2.0^-66,
)

"""
SISAIoc block, 22 bits: `t_op`(11, LSB 300 s), `SISAI_ocb`(5), `SISAI_oc1`(3),
`SISAI_oc2`(3) — Figure 6-14. The three indices are raw broadcast values: the
ICDs defer the accuracy semantics to a future update (B1C §7.16).

`t_op` is scaled to seconds although no BDS table states its factor. The ICDs define
it as "the time of week for data prediction" (B2b §7.15 item 5); 11 bits in seconds
cannot span a week, while LSB 300 s gives 0..614100, the convention of `t_0c`/`t_0e`
(Table 7-5) and of GPS CNAV's same-named `t_op` (IS-GPS-705J). Leaving it raw would
give the same field name in two units (CONTEXT.md naming doctrine, rule 1).
"""
beidou_sisai_oc_block(word::Unsigned, word_length::Int, start::Int) = (;
    t_op = Int64(get_bits(word, word_length, start, 11)) * 300,
    SISAI_ocb = Int64(get_bits(word, word_length, start + 11, 5)),
    SISAI_oc1 = Int64(get_bits(word, word_length, start + 16, 3)),
    SISAI_oc2 = Int64(get_bits(word, word_length, start + 19, 3)),
)

"""
BDGIM ionospheric block, 74 bits: `α₁`(10) then `α₂`…`α₉`(8 each), all LSB 2⁻³
TECu — Figure 6-16, Table 7-10.

Signedness is not uniform: `α₁`, `α₃`, `α₄` and `α₅` are unsigned, `α₂` and
`α₆`…`α₉` are two's complement, and `α₅` carries the *negative* scale factor −2⁻³.
"""
beidou_bdgim_block(word::Unsigned, word_length::Int, start::Int) = (;
    α_bdgim_1 = Int64(get_bits(word, word_length, start, 10)) * 2.0^-3,
    α_bdgim_2 = get_twos_complement_num(word, word_length, start + 10, 8) * 2.0^-3,
    α_bdgim_3 = Int64(get_bits(word, word_length, start + 18, 8)) * 2.0^-3,
    α_bdgim_4 = Int64(get_bits(word, word_length, start + 26, 8)) * 2.0^-3,
    α_bdgim_5 = Int64(get_bits(word, word_length, start + 34, 8)) * -(2.0^-3),
    α_bdgim_6 = get_twos_complement_num(word, word_length, start + 42, 8) * 2.0^-3,
    α_bdgim_7 = get_twos_complement_num(word, word_length, start + 50, 8) * 2.0^-3,
    α_bdgim_8 = get_twos_complement_num(word, word_length, start + 58, 8) * 2.0^-3,
    α_bdgim_9 = get_twos_complement_num(word, word_length, start + 66, 8) * 2.0^-3,
)

"""
BDT-UTC offset block, 97 bits: `A_0UTC`(16, 2⁻³⁵), `A_1UTC`(13, 2⁻⁵¹),
`A_2UTC`(7, 2⁻⁶⁸), `Δt_LS`(8), `t_0t`(16, LSB 2⁴ s), `WN_0t`(13),
`WN_LSF`(13), `DN`(3), `Δt_LSF`(8) — Figure 6-17, Table 7-20.
"""
beidou_bdt_utc_block(word::Unsigned, word_length::Int, start::Int) = (;
    A_0UTC = get_twos_complement_num(word, word_length, start, 16) * 2.0^-35,
    A_1UTC = get_twos_complement_num(word, word_length, start + 16, 13) * 2.0^-51,
    A_2UTC = get_twos_complement_num(word, word_length, start + 29, 7) * 2.0^-68,
    Δt_LS = Int64(get_twos_complement_num(word, word_length, start + 36, 8)),
    t_0t = Int64(get_bits(word, word_length, start + 44, 16)) * 2^4,
    WN_0t = Int64(get_bits(word, word_length, start + 60, 13)),
    WN_LSF = Int64(get_bits(word, word_length, start + 73, 13)),
    DN = Int64(get_bits(word, word_length, start + 86, 3)),
    Δt_LSF = Int64(get_twos_complement_num(word, word_length, start + 89, 8)),
)

# B-CNAV `WN_0t` and `WN_LSF` are 13 bits, the width of `WN` (BDS-SIS-ICD-B1C-1.0
# Table 7-20): see `get_utc_parameters`.
const BEIDOU_CNAV_UTC_WN_MODULUS = 8192
const BEIDOU_CNAV_LSF_WN_MODULUS = 8192

"""
Earth orientation parameter block, 138 bits: `t_EOP`(16, LSB 2⁴ s), `PM_X`(21,
2⁻²⁰ arcsec), `PM_X_dot`(15, 2⁻²¹), `PM_Y`(21, 2⁻²⁰), `PM_Y_dot`(15, 2⁻²¹),
`ΔUT1`(31, 2⁻²⁴ s), `ΔUT1_dot`(19, 2⁻²⁵) — Figure 6-19, Table 7-18.
"""
beidou_eop_block(word::Unsigned, word_length::Int, start::Int) = (;
    t_EOP = Int64(get_bits(word, word_length, start, 16)) * 2^4,
    PM_X = get_twos_complement_num(word, word_length, start + 16, 21) * 2.0^-20,
    PM_X_dot = get_twos_complement_num(word, word_length, start + 37, 15) * 2.0^-21,
    PM_Y = get_twos_complement_num(word, word_length, start + 52, 21) * 2.0^-20,
    PM_Y_dot = get_twos_complement_num(word, word_length, start + 73, 15) * 2.0^-21,
    ΔUT1 = get_twos_complement_num(word, word_length, start + 88, 31) * 2.0^-24,
    ΔUT1_dot = get_twos_complement_num(word, word_length, start + 119, 19) * 2.0^-25,
)

"""
    beidou_orbit_class(sat_type) -> Union{Nothing,OrbitClass}

Map a broadcast B-CNAV `sat_type` onto an [`OrbitClass`](@ref): 1 GEO, 2 IGSO,
3 MEO, with the reserved 0 — and an undecoded `nothing` — reporting `nothing`.

The same screen as `is_known_sat_type`, in the form [`get_orbit_class`](@ref)
answers in. B1I/B3I do not reach this: D1/D2 NAV has no orbit-type field (see
`dnav.jl`).
"""
beidou_orbit_class(sat_type) =
    sat_type == 1 ? geostationary_orbit :
    sat_type == 2 ? inclined_geosynchronous_orbit :
    sat_type == 3 ? medium_earth_orbit : nothing

# Orbit class from the broadcast `sat_type` (Ephemeris I block); `nothing` until decoded
# and for the reserved code. Dispatched on the B-CNAV family type because only B-CNAV
# carries `sat_type`; B1I/B3I use the GEO PRN partition (see `beidou/dnav.jl`).
get_orbit_class(state::GNSSDecoderState{<:AbstractBeiDouCNAVData}) =
    beidou_orbit_class(state.data.sat_type)

# `WN_0BGTO` is 13 bits, the width of the BDT week number, so resolving it matters only
# across a week-number rollover. See `resolve_reference_week`.
const BEIDOU_BGTO_WN_MODULUS = 8192

"""
    beidou_bgto_target(GNSS_ID) -> Union{Nothing,TimeSystem}

Map a broadcast BGTO `GNSS_ID` onto the time system the offset refers to: 1 GPS,
2 Galileo, 3 GLONASS, 0 "not available" (§7.13.1 of the B1C, B2a and B2b ICDs).

GLONASS, code 0 and an undecoded `nothing` report `nothing`: `GNSSSignals` defines no
GLONASS `TimeSystem`, so a GLONASS offset is read off the data fields directly. See
[`get_broadcast_time_offset`](@ref).
"""
beidou_bgto_target(GNSS_ID) = GNSS_ID == 1 ? GPST() : GNSS_ID == 2 ? GST() : nothing

"""
    beidou_bgto_offset(state::GNSSDecoderState, target::TimeSystem) -> Union{Nothing,TimeOffsetParameters}

Normalise the flat BGTO parameter set that B2a (message type 33) and B2b (message
type 40) broadcast into a `TimeOffsetParameters`, or `nothing` when the set is absent,
marked unavailable, or refers to a system other than `target`.

Both signals carry one set at a time, tagged with its system, unlike B1C, which keys a
set per system (see [`BeiDouB1CBGTO`](@ref GNSSDecoder.BeiDouB1CBGTO)). The
`GNSS_ID == 0` "not available" sentinel is screened here because both decoders keep the
raw field; B1C screens it at decode, having no key to file an unavailable set under.
"""
function beidou_bgto_offset(state::GNSSDecoderState, target::TimeSystem)
    data = state.data
    beidou_bgto_target(data.GNSS_ID) === target || return nothing
    # Decoded together, so all present or all `nothing`; each is checked to narrow its
    # type for the call below (see `galileo_ggto_offset`).
    (; A_0BGTO, A_1BGTO, A_2BGTO, t_0BGTO, WN_0BGTO) = data
    (
        isnothing(A_0BGTO) ||
        isnothing(A_1BGTO) ||
        isnothing(A_2BGTO) ||
        isnothing(t_0BGTO) ||
        isnothing(WN_0BGTO)
    ) && return nothing
    assemble_time_offset_parameters(
        target,
        A_0BGTO,
        A_1BGTO,
        A_2BGTO,
        t_0BGTO,
        WN_0BGTO,
        data.WN,
        BEIDOU_BGTO_WN_MODULUS,
    )
end

"""
BDT-GNSS time offset block, 68 bits: `GNSS_ID`(3), `WN_0BGTO`(13),
`t_0BGTO`(16, LSB 2⁴ s), `A_0BGTO`(16, 2⁻³⁵), `A_1BGTO`(13, 2⁻⁵¹),
`A_2BGTO`(7, 2⁻⁶⁸) — Figure 6-20, Table 7-21. `GNSS_ID` 0 means the parameters
are unavailable (1 GPS, 2 Galileo, 3 GLONASS); the caller decides what to store.
"""
beidou_bgto_block(word::Unsigned, word_length::Int, start::Int) = (;
    GNSS_ID = Int64(get_bits(word, word_length, start, 3)),
    WN_0BGTO = Int64(get_bits(word, word_length, start + 3, 13)),
    t_0BGTO = Int64(get_bits(word, word_length, start + 16, 16)) * 2^4,
    A_0BGTO = get_twos_complement_num(word, word_length, start + 32, 16) * 2.0^-35,
    A_1BGTO = get_twos_complement_num(word, word_length, start + 48, 13) * 2.0^-51,
    A_2BGTO = get_twos_complement_num(word, word_length, start + 61, 7) * 2.0^-68,
)

"""
    beidou_midi_almanac(word, word_length, start) -> Union{Nothing,BeiDouMidiAlmanac}

One 156-bit midi-almanac block (Figure 6-21, Table 7-13), or `nothing` for an empty
slot (`PRN_a == 0`); following blocks may still be valid.
"""
function beidou_midi_almanac(word::Unsigned, word_length::Int, start::Int)
    PRN_a = Int(get_bits(word, word_length, start, 6))
    PRN_a == 0 && return nothing
    BeiDouMidiAlmanac(;
        PRN_a,
        sat_type = Int(get_bits(word, word_length, start + 6, 2)),
        WN_a = Int(get_bits(word, word_length, start + 8, 13)),
        t_0a = Int(get_bits(word, word_length, start + 21, 8)) * 2^12,
        e = Int(get_bits(word, word_length, start + 29, 11)) * 2.0^-16,
        δi = get_twos_complement_num(word, word_length, start + 40, 11) * 2.0^-14 * GNSS_PI,
        sqrt_A = Int(get_bits(word, word_length, start + 51, 17)) * 2.0^-4,
        Ω_0 = get_twos_complement_num(word, word_length, start + 68, 16) *
              2.0^-15 *
              GNSS_PI,
        Ω_dot = get_twos_complement_num(word, word_length, start + 84, 11) *
                2.0^-33 *
                GNSS_PI,
        ω = get_twos_complement_num(word, word_length, start + 95, 16) * 2.0^-15 * GNSS_PI,
        M_0 = get_twos_complement_num(word, word_length, start + 111, 16) *
              2.0^-15 *
              GNSS_PI,
        a_f0 = get_twos_complement_num(word, word_length, start + 127, 11) * 2.0^-20,
        a_f1 = get_twos_complement_num(word, word_length, start + 138, 10) * 2.0^-37,
        health = Int(get_bits(word, word_length, start + 148, 8)),
    )
end

"""
    beidou_reduced_almanac(word, word_length, start, WN_a, t_0a)
        -> Union{Nothing,BeiDouReducedAlmanac}

One 38-bit reduced-almanac block (Figure 6-18, Table 7-16), or `nothing` for the
`PRN_a == 0` empty-slot encoding. The block has no epoch, so the almanac reference
week and time of the carrying page are passed in and copied into the record.
"""
function beidou_reduced_almanac(
    word::Unsigned,
    word_length::Int,
    start::Int,
    WN_a::Integer,
    t_0a::Integer,
)
    PRN_a = Int(get_bits(word, word_length, start, 6))
    PRN_a == 0 && return nothing
    BeiDouReducedAlmanac(;
        PRN_a,
        sat_type = Int(get_bits(word, word_length, start + 6, 2)),
        WN_a = Int(WN_a),
        t_0a = Int(t_0a),
        δA = get_twos_complement_num(word, word_length, start + 8, 8) * 2.0^9,
        Ω_0 = get_twos_complement_num(word, word_length, start + 16, 7) * 2.0^-6 * GNSS_PI,
        Φ_0 = get_twos_complement_num(word, word_length, start + 23, 7) * 2.0^-6 * GNSS_PI,
        health = Int(get_bits(word, word_length, start + 30, 8)),
    )
end

"""
True when `IODE` equals the 8 LSBs of `IODC` (the B-CNAV "matched pair" rule),
`false` while `IODC` is `nothing`.
"""
# Function barrier: Julia 1.10 does not narrow a field read through `isnothing`, so
# `data.IODC & 0xff` would dispatch dynamically.
beidou_iode_matches_iodc(IODE, IODC::Integer) = IODE == IODC & 0xff
beidou_iode_matches_iodc(IODE, ::Nothing) = false
