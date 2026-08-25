# Definitions shared across BeiDou signals (B1I/B3I D1-D2 NAV and the B-CNAV
# family on B1C, B2a, and B2b): the BDCS physical constants and the GEO-orbit
# PRN partition that selects between the D1 and D2 legacy message formats.
# Each individual signal's framing, message layout, and parser live in their
# own files (`dnav.jl` + `b1i.jl`/`b3i.jl`, `b1c.jl`, `b2a.jl`, `b2b.jl`).
# This mirrors how `src/gnss.jl` holds the definitions shared across *all*
# signals and `src/galileo/galileo.jl` those shared across Galileo.

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
  - `t_oa::Int`: Almanac reference time of week (seconds, LSB 2¹²).
  - `e::Float64`: Eccentricity (dimensionless).
  - `δi::Float64`: Inclination delta from the `sat_type` reference (rad).
  - `sqrt_A::Float64`: Square root of the semi-major axis (√m).
  - `Ω_0::Float64`: Longitude of ascending node at weekly epoch (rad).
  - `Ω_dot::Float64`: Rate of right ascension (rad/s).
  - `ω::Float64`: Argument of perigee (rad).
  - `M_0::Float64`: Mean anomaly at reference time (rad).
  - `a_f0::Float64`, `a_f1::Float64`: Satellite clock bias / drift (s, s/s).
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
    t_oa::Int
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
  - `t_oa::Int`: Almanac reference time of week (seconds, from the carrying page/message).
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
    t_oa::Int
    δA::Float64
    Ω_0::Float64
    Φ_0::Float64
    health::Int
end

"""
$(SIGNATURES)

Whether a BeiDou PRN is one of the GEO satellites (PRN 1-5 and 59-63).

The distinction selects the legacy navigation message format on B1I and B3I:
MEO/IGSO satellites broadcast the D1 message (50 bps, with the NH20 secondary
code), GEO satellites the D2 message (500 bps, no secondary code)
(BDS-SIS-ICD-B1I-3.0 §5.1.1 / Table 4-1). The B-CNAV messages (B1C, B2a, B2b)
are broadcast by BDS-3 MEO/IGSO satellites only and do not use this partition.
"""
is_beidou_geo(prn::Integer) = prn in 1:5 || prn in 59:63
