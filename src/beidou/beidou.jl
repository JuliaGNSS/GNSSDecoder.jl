# Definitions shared across BeiDou signals (B1I/B3I D1-D2 NAV and the B-CNAV
# family on B1C, B2a, and B2b): the BDCS physical constants, shared by all five
# decoders, and the midi and reduced almanac records, shared by the three
# B-CNAV ones. Each individual signal's framing, message layout, and parser
# live in their own files (`dnav.jl` + `b1i.jl`/`b3i.jl`, `b1c.jl`, `b2a.jl`,
# `b2b.jl`). This mirrors how `src/gnss.jl` holds the definitions shared across
# *all* signals and `src/galileo/galileo.jl` those shared across Galileo.
#
# The test for belonging here is that no single signal file can own it. The
# GEO-orbit PRN partition looks like it qualifies and does not: it selects
# between the D1 and D2 *message formats*, so it is used by `dnav.jl` alone and
# lives there.

# BDCS (BeiDou Coordinate System) constants specific to BeiDou (BDS-SIS-ICD-
# B1I-3.0 Table 3-1 / §5.2.4.12; identical across the B1I, B3I, B1C, B2a, and
# B2b ICDs). Only the two that differ from the GPS WGS-84 values live here; π
# and the speed of light are shared package-wide (`GNSS_PI`, `SPEED_OF_LIGHT`
# in `gnss.jl`). Note the Earth rotation rate genuinely differs from the
# WGS-84/GTRF value used by GPS and Galileo (7.2921151467e-5 rad/s), so the
# package-wide `EARTH_ROTATION_RATE` must *not* be used for BeiDou.
const BEIDOU_μ = 3.986004418e14                # geocentric gravitational constant (m³/s²)
const BEIDOU_EARTH_ROTATION_RATE = 7.2921150e-5 # Earth rotation rate of BDCS (rad/s)

# Whether a broadcast B-CNAV satellite orbit type names an orbit class. The
# 2-bit field encodes 1 = GEO, 2 = IGSO, 3 = MEO and reserves 0 (Table 7-6 of
# the B1C, B2a, and B2b ICDs).
#
# This gates positioning, because `sat_type` is what selects the reference the
# broadcast ephemeris is expressed against: `A = A_ref + ΔA` with
# `A_ref = 27 906 100 m` for MEO and `42 162 200 m` for IGSO/GEO, and
# `i = i_ref + δi` in the almanacs. A satellite broadcasting the reserved code
# leaves its own semi-major axis unknowable — the two references differ by
# 14 256 100 m — so it must not be used rather than be positioned against a
# guessed orbit class. B1I and B3I are unaffected: D1/D2 NAV broadcasts `sqrt_A`
# outright and carries no orbit-type field.
#
# The almanac records keep the raw field either way; a reserved orbit type in
# somebody else's almanac is not a reason to discard the rest of the page.
is_known_sat_type(sat_type::Integer) = sat_type in 1:3
is_known_sat_type(::Nothing) = false

"""
    BeiDouMidiAlmanac

Midi almanac for one BeiDou satellite.

The BDS-3 midi almanac is one 156-bit block whose layout, scale factors, and
reference values are identical across the three B-CNAV messages, so the same
record is produced by the B1C decoder (B-CNAV1 subframe 3 page type 4,
BDS-SIS-ICD-B1C-1.0 Figure 6-21 / Table 7-13), the B2a decoder (B-CNAV2
message type 40, BDS-SIS-ICD-B2a-1.0 Figure 6-20 / Table 7-13), and the B2b
decoder (B-CNAV3 message type 40, BDS-SIS-ICD-B2b-1.0 Figure 6-15 /
Table 7-11) — mirroring how [`GalileoAlmanac`](@ref) is shared by the I/NAV
and F/NAV decoders. Each block is complete in itself and keyed by `PRN_a` in
the per-signal data containers.

Reference value to apply (ICD notes to the midi-almanac tables): inclination
is `i = i_ref + δi` with `i_ref = 0.30π rad` (54°) for MEO/IGSO satellites and
`i_ref = 0.0` for GEO satellites, selected by `sat_type`. Semi-circle fields
are converted to radians on decode.

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
    The midi-almanac terms are the one BeiDou clock pair the ICDs already spell
    `af0`/`af1` (B2b-1.0 Table 7-11, B1C/B2a-1.0 Table 7-13) — unlike the
    clock-correction block, which they write `a_0`/`a_1`/`a_2`.
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

The BDS-3 reduced almanac is one 38-bit block whose layout, scale factors,
and reference values are identical across the three B-CNAV messages, so the
same record is produced by the B1C decoder (B-CNAV1 subframe 3 page type 2,
four blocks per page, BDS-SIS-ICD-B1C-1.0 Figure 6-18 / Table 7-16), the B2a
decoder (B-CNAV2 message types 31 and 33, BDS-SIS-ICD-B2a-1.0 Figure 6-17 /
Table 7-16), and the B2b decoder (B-CNAV3 message type 40, five blocks per
message, BDS-SIS-ICD-B2b-1.0 Figure 6-12 / Table 7-14). The 38-bit block
itself carries no epoch; the almanac reference week/time broadcast alongside
it in the carrying page/message is copied into every record, so each record
is complete in itself and keyed by `PRN_a` in the per-signal data containers.

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
# B1C, B2a and B2b are three different messages built from the same eight data
# blocks. Each ICD defines a block once — layout, widths, signedness, scale
# factors — and then places it at whatever bit offset its own message needs, so
# the only thing that differs between the three decoders is that offset.
#
# Parsing each block once is what keeps the three from drifting. They already had:
# `δA` was `* 2.0^9` in two files and `* 512.0` in the third, `t_0t` and `t_EOP`
# were `* 2^4` in two and `* 16` in the third, and the integer conversions were
# variously `Int`, `Int64` and `Float64` — the same numbers by luck rather than by
# construction, in code verified against three separate ICD documents.
#
# Each function takes the packed message word, its logical bit length, and the
# 1-based position of the block's first bit, and returns a `NamedTuple` keyed by
# container field names — which the three containers agree on, that being what
# the package-wide naming rule bought. Callers splat it into their container:
#
#     BeiDouB2bData(raw; beidou_clock_block(word, word_length, 44)...)
#
# References are given per block; the B1C figure/table numbers are cited, with the
# B2a-1.0 and B2b-1.0 counterparts named at each call site.

"""
Clock correction block, 69 bits: `t_oc`(11, LSB 300 s), `a_0`(25, 2⁻³⁴ s),
`a_1`(22, 2⁻⁵⁰ s/s), `a_2`(11, 2⁻⁶⁶ s/s²) — BDS-SIS-ICD-B1C-1.0 Figure 6-13,
Table 7-5.
"""
beidou_clock_block(word::Unsigned, word_length::Int, start::Int) = (;
    t_0c = Int64(get_bits(word, word_length, start, 11)) * 300,
    a_f0 = get_twos_complement_num(word, word_length, start + 11, 25) * 2.0^-34,
    a_f1 = get_twos_complement_num(word, word_length, start + 36, 22) * 2.0^-50,
    a_f2 = get_twos_complement_num(word, word_length, start + 58, 11) * 2.0^-66,
)

"""
SISAIoc block, 22 bits: `t_op`(11), `SISAI_ocb`(5), `SISAI_oc1`(3),
`SISAI_oc2`(3) — Figure 6-14. Raw broadcast values: the ICDs defer the accuracy
semantics to a future update (B1C §7.16).
"""
beidou_sisai_oc_block(word::Unsigned, word_length::Int, start::Int) = (;
    t_op = Int64(get_bits(word, word_length, start, 11)),
    SISAI_ocb = Int64(get_bits(word, word_length, start + 11, 5)),
    SISAI_oc1 = Int64(get_bits(word, word_length, start + 16, 3)),
    SISAI_oc2 = Int64(get_bits(word, word_length, start + 19, 3)),
)

"""
BDGIM ionospheric block, 74 bits: `α₁`(10) then `α₂`…`α₉`(8 each), all LSB 2⁻³
TECu — Figure 6-16, Table 7-10.

The signedness is not uniform and is easy to get wrong: `α₁`, `α₃`, `α₄` and `α₅`
are unsigned, `α₂` and `α₆`…`α₉` are two's complement, and `α₅` carries the
*negative* scale factor −2⁻³.
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

The same screen `is_known_sat_type` applies to positioning, in the form
[`get_orbit_class`](@ref) answers in. B1I and B3I do not reach this: D1/D2 NAV
broadcasts no orbit-type field (see `dnav.jl`).
"""
beidou_orbit_class(sat_type) =
    sat_type == 1 ? geostationary_orbit :
    sat_type == 2 ? inclined_geosynchronous_orbit :
    sat_type == 3 ? medium_earth_orbit : nothing

# `WN_0BGTO` is 13 bits, the same width as the BDT week number it is subtracted
# from, so resolving it changes nothing except across a week-number rollover —
# which is exactly when it matters. See `resolve_reference_week`.
const BEIDOU_BGTO_WN_MODULUS = 8192

"""
    beidou_bgto_target(GNSS_ID) -> Union{Nothing,TimeSystem}

Map a broadcast BGTO `GNSS_ID` onto the time scale the offset refers to: 1 GPS,
2 Galileo, 3 GLONASS, 0 "not available" (§7.13.1 of the B1C, B2a and B2b ICDs).

GLONASS reports `nothing`, as do code 0 and an undecoded `nothing`: `GNSSSignals`
defines no GLONASS `TimeSystem`, so a GLONASS offset has no target to be asked
for and is read off the data fields directly. See [`get_time_offset`](@ref).
"""
beidou_bgto_target(GNSS_ID) = GNSS_ID == 1 ? GPST() : GNSS_ID == 2 ? GST() : nothing

"""
    beidou_bgto_offset(data, target) -> Union{Nothing,GNSSTimeOffset}

Normalise the single flat BGTO parameter set that B2a (message type 33) and B2b
(message type 40) broadcast into a [`GNSSTimeOffset`](@ref), or `nothing` when
the set is absent, marked unavailable, or refers to a system other than
`target`.

Both signals carry one set at a time, tagged with the system it refers to, so
which target is answerable is the satellite's choice at the moment of asking —
unlike B1C, which keys a set per system (see [`BeiDouB1CBGTO`](@ref)). The
`GNSS_ID == 0` "not available" sentinel is screened here rather than at decode,
because both decoders keep the raw field; B1C screens it at decode instead,
since it has a keyed store and no key to file an unavailable set under.
"""
function beidou_bgto_offset(state::GNSSDecoderState, target::TimeSystem)
    data = state.data
    beidou_bgto_target(data.GNSS_ID) === target || return nothing
    isnothing(data.A_0BGTO) && return nothing
    broadcast_time_offset(
        state,
        target,
        data.A_0BGTO,
        data.A_1BGTO,
        data.A_2BGTO;
        t_0 = data.t_0BGTO,
        WN_0 = data.WN_0BGTO,
        WN = data.WN,
        WN_0_modulus = BEIDOU_BGTO_WN_MODULUS,
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
    beidou_midi_almanac(word, word_length, start, PI) -> Union{Nothing,BeiDouMidiAlmanac}

One 156-bit midi-almanac block (Figure 6-21, Table 7-13), or `nothing` for the
`PRN_a == 0` encoding that marks an empty slot — no almanac was broadcast in this
frame, and the following blocks may still be valid.
"""
function beidou_midi_almanac(word::Unsigned, word_length::Int, start::Int, PI::Float64)
    PRN_a = Int(get_bits(word, word_length, start, 6))
    PRN_a == 0 && return nothing
    BeiDouMidiAlmanac(;
        PRN_a,
        sat_type = Int(get_bits(word, word_length, start + 6, 2)),
        WN_a = Int(get_bits(word, word_length, start + 8, 13)),
        t_0a = Int(get_bits(word, word_length, start + 21, 8)) * 2^12,
        e = Int(get_bits(word, word_length, start + 29, 11)) * 2.0^-16,
        δi = get_twos_complement_num(word, word_length, start + 40, 11) * 2.0^-14 * PI,
        sqrt_A = Int(get_bits(word, word_length, start + 51, 17)) * 2.0^-4,
        Ω_0 = get_twos_complement_num(word, word_length, start + 68, 16) * 2.0^-15 * PI,
        Ω_dot = get_twos_complement_num(word, word_length, start + 84, 11) * 2.0^-33 * PI,
        ω = get_twos_complement_num(word, word_length, start + 95, 16) * 2.0^-15 * PI,
        M_0 = get_twos_complement_num(word, word_length, start + 111, 16) * 2.0^-15 * PI,
        a_f0 = get_twos_complement_num(word, word_length, start + 127, 11) * 2.0^-20,
        a_f1 = get_twos_complement_num(word, word_length, start + 138, 10) * 2.0^-37,
        health = Int(get_bits(word, word_length, start + 148, 8)),
    )
end

"""
    beidou_reduced_almanac(word, word_length, start, WN_a, t_0a, PI)
        -> Union{Nothing,BeiDouReducedAlmanac}

One 38-bit reduced-almanac block (Figure 6-18, Table 7-16), or `nothing` for the
`PRN_a == 0` empty-slot encoding. The block carries no epoch of its own, so the
almanac reference week and time of the carrying page are passed in and copied
into the record.
"""
function beidou_reduced_almanac(
    word::Unsigned,
    word_length::Int,
    start::Int,
    WN_a::Integer,
    t_0a::Integer,
    PI::Float64,
)
    PRN_a = Int(get_bits(word, word_length, start, 6))
    PRN_a == 0 && return nothing
    BeiDouReducedAlmanac(;
        PRN_a,
        sat_type = Int(get_bits(word, word_length, start + 6, 2)),
        WN_a = Int(WN_a),
        t_0a = Int(t_0a),
        δA = get_twos_complement_num(word, word_length, start + 8, 8) * 2.0^9,
        Ω_0 = get_twos_complement_num(word, word_length, start + 16, 7) * 2.0^-6 * PI,
        Φ_0 = get_twos_complement_num(word, word_length, start + 23, 7) * 2.0^-6 * PI,
        health = Int(get_bits(word, word_length, start + 30, 8)),
    )
end
