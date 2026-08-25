# GPS L1C-D (CNAV-2) decoder — TOI sync + subframe 2 (issue #38).
#
# The L1C-D civil data signal carries the CNAV-2 message as 18-second frames
# of 1800 symbols at 100 sps:
#
#   - Subframe 1: 52 symbols — a BCH(51,8)+parity encoding of a 9-bit TOI
#     count (0..399), used purely for frame sync (see `src/bch_toi.jl`).
#   - Subframes 2 + 3: 1748 symbols, block-interleaved (38×46, see
#     `src/deinterleave.jl`) over a 1200-symbol rate-½ LDPC codeword for
#     subframe 2 (600 info bits) followed by a 548-symbol rate-½ LDPC
#     codeword for subframe 3 (274 info bits).
#
# This slice decodes subframe 2 fully (clock + ephemeris, IS-GPS-800J
# Figure 3.5-1 / Table 3.5-1) and LDPC-decodes + CRC-checks subframe 3 but
# only records that a page was received; the per-page subframe-3 field
# parsing lands in issue #39.
#
# Frame sync follows the generic streaming framework in `src/gnss.jl`: the
# soft-symbol `CircularDeque` is sized to one full subframe (1800) plus the
# 52-symbol subframe-1 segment of the *next* frame (1852 total) so the BCH
# match can be confirmed at both ends of the window before locking.

# `import Aff3ct` for the qualified `Aff3ct.decode` call (it would clash with
# `GNSSDecoder.decode`); bring only the LDPC constructors into scope by name.
import Aff3ct
using Aff3ct: LDPCMatrix, LDPCBPDecoder

"""
Total symbols accumulated before a sync attempt: 1800 (one frame) + 52 (next subframe-1).
"""
const L1C_D_FRAME_LENGTH = 1800
const L1C_D_SUBFRAME1_LENGTH = 52
const L1C_D_WINDOW_LENGTH = L1C_D_FRAME_LENGTH + L1C_D_SUBFRAME1_LENGTH  # 1852

# Subframe 2/3 channel-symbol counts (IS-GPS-800J §3.2.3).
const L1C_D_SF2_SYMBOLS = 1200
const L1C_D_SF3_SYMBOLS = 548
const L1C_D_PAYLOAD_SYMBOLS = L1C_D_SF2_SYMBOLS + L1C_D_SF3_SYMBOLS  # 1748
const L1C_D_INTERLEAVER_ROWS = 38
const L1C_D_INTERLEAVER_COLS = 46  # 38 * 46 == 1748

# LDPC info-block sizes (rate ½): N channel symbols -> K info bits.
const L1C_D_SF2_INFO_BITS = 600
const L1C_D_SF3_INFO_BITS = 274
const L1C_D_CRC_BITS = 24

# Packed-word types for bit-field extraction. After the LDPC decode + CRC check
# the info bits are packed MSB-first (LSB = last bit) into a wide integer and
# parsed with the shared `get_bits`/`get_twos_complement_num`/`get_bit` helpers
# (src/bit_fiddling.jl), matching the GPS L1 C/A and Galileo E1B decoders (#48).
#
# SF2's 600 info bits go in a `UInt600` and SF3's 274 in a `UInt288`, both
# defined in `gnss.jl` (BeiDou B1C reaches for the same two widths, so neither
# belongs to a signal file): `get_bits(word, 274, …)` addresses the 274 logical
# bits regardless of the wider storage as long as they are right-aligned.

# Semi-major axis reference (IS-GPS-800J Table 3.5-1 footnote, meters).
const L1C_D_A_REF = 26_559_710.0
# Rate-of-right-ascension reference (IS-GPS-800J Table 3.5-1, semi-circles/sec).
const L1C_D_OMEGA_DOT_REF = -2.6e-9

# Subframe-3 page numbers (IS-GPS-800J §3.5.4, figure-by-figure). The 6-bit
# page-number field lives in bits 9-14 of every SF3 page (bits 1-8 are PRN).
const L1C_D_SF3_PAGE_UTC_IONO = 1   # Figure 3.5-2: UTC + Klobuchar iono + ISC
const L1C_D_SF3_PAGE_GGTO_EOP = 2   # Figure 3.5-3: GGTO + Earth orientation
const L1C_D_SF3_PAGE_REDUCED_ALMANAC = 3  # Figure 3.5-4: 6 reduced-almanac packets
const L1C_D_SF3_PAGE_MIDI_ALMANAC = 4     # Figure 3.5-5: one Midi almanac
const L1C_D_SF3_PAGE_DIFF_CORRECTION = 5  # Figure 3.5-6: clock+ephemeris DC
const L1C_D_SF3_PAGE_TEXT = 6             # Figure 3.5-7: 29 ASCII characters

# Reduced-almanac references (IS-GPS-800J Table 3.5-6 footnotes): δi relative to
# i₀ = 0.30 semi-circles (54°); with δi the almanac's inclination is 55°.
# δA is relative to A_ref.
const L1C_D_REDUCED_ALMANAC_DELTA_I_REF = 0.0056  # semi-circles
const L1C_D_REDUCED_ALMANAC_I_REF = 0.30          # semi-circles
# Midi-almanac inclination reference (IS-GPS-800J Table 3.5-7 footnote): δi
# relative to i₀ = 0.30 semi-circles.
const L1C_D_MIDI_ALMANAC_I_REF = 0.30  # semi-circles

"""
    GPSL1C_DConstants

WGS 84 constants and CNAV-2 message structure parameters for GPS L1C-D decoding.

The frame is modelled through the generic streaming framework: `preamble_length`
is the 52-symbol subframe-1 BCH segment of the *next* frame retained at the tail
of the sync window, and `syncro_sequence_length` is the 1800-symbol frame that is
drained once a subframe is decoded.

# Fields

$(TYPEDFIELDS)

# Reference

IS-GPS-800J, Sections 3.2 and 3.5, Table 3.5-1.
"""
Base.@kwdef struct GPSL1C_DConstants <: AbstractGNSSConstants
    """
    Frame length drained after each decoded subframe (1800 symbols)
    """
    syncro_sequence_length::Int = L1C_D_FRAME_LENGTH
    """
    Trailing next-frame subframe-1 BCH segment retained for sync (52 symbols)
    """
    preamble_length::Int = L1C_D_SUBFRAME1_LENGTH
    """
    Mathematical constant π (IS-GPS-800J)
    """
    PI::Float64 = GNSS_PI
    """
    WGS 84 Earth rotation rate (rad/s)
    """
    Ω_dot_e::Float64 = EARTH_ROTATION_RATE
    """
    Speed of light (m/s)
    """
    c::Float64 = SPEED_OF_LIGHT
    """
    WGS 84 Earth gravitational parameter (m³/s²)
    """
    μ::Float64 = 3.986005e14
    """
    Relativistic correction constant (s/√m)
    """
    F::Float64 = -4.442807633e-10
end

"""
    GPSL1C_DReducedAlmanac

One satellite's reduced-almanac packet from subframe 3, page 3 (IS-GPS-800J
Figure 3.5-9, Table 3.5-6).

The reduced almanac gives a very coarse ephemeris for satellite selection. Each
page-3 carries six 33-bit packets; this struct holds one decoded packet plus the
page-level almanac reference week/time. A reduced almanac is *complete in a
single page* — there is no IOD-driven multi-page chaining like Galileo's word
types 7-10 — so `GPSL1C_DData.reduced_almanacs` entries are inserted
whole, keyed by `PRN_a`. Reduced and Midi almanacs use *separate* structs (their
field sets barely overlap); they share the `Dictionary` pattern.

Reference values to apply (Table 3.5-6 footnotes): `e = 0`,
`δi = +0.0056 semi-circles` (so `i₀ = 0.30 sc = 54°` and `i₀ + δi = 55°`),
`Ω̇ = -2.6e-9 semi-circles/s`, `A = A_ref + δA` with `A_ref = 26 559 710 m`,
`Φ₀ = M₀ + ω`. Semi-circle fields are converted to radians on decode.

# Fields

  - `PRN_a::Int`: Almanac satellite PRN (8-bit field; `0` marks an empty
    packet, IS-GPS-800J §3.5.4.3.5.1.1). The CNAV reduced almanac carries the
    same quantity in 6 bits — see [`GPSCNAVReducedAlmanac`](@ref).
  - `WN_a::Int`: Almanac reference week number (mod 8192).
  - `t_0a::Int`: Almanac reference time of week (seconds).
  - `δA::Float64`: Semi-major-axis delta from `A_ref` (meters).
  - `Ω_0::Float64`: Longitude of ascending node at weekly epoch (rad).
  - `Φ_0::Float64`: Argument of latitude at reference time, `M₀+ω` (rad).
  - `l1_health::Bool`, `l2_health::Bool`, `l5_health::Bool`: per-band health
    (false = OK, true = some/all signals bad).

# Reference

IS-GPS-800J, Figure 3.5-4 / Figure 3.5-9 / Table 3.5-6.
"""
Base.@kwdef struct GPSL1C_DReducedAlmanac
    PRN_a::Int
    WN_a::Int
    t_0a::Int
    δA::Float64
    Ω_0::Float64
    Φ_0::Float64
    l1_health::Bool
    l2_health::Bool
    l5_health::Bool
end

"""
    GPSL1C_DMidiAlmanac

One satellite's Midi almanac from subframe 3, page 4 (IS-GPS-800J Figure 3.5-5,
Table 3.5-7).

The Midi almanac is a medium-precision single-SV almanac. Each page-4 carries
exactly one SV's almanac, complete in that single page (no multi-page chaining),
so `GPSL1C_DData.midi_almanacs` entries are inserted whole, keyed by
`PRN_a`. Inclination is `δi` relative to `i₀ = 0.30 semi-circles` (54°);
semi-circle fields are converted to radians on decode.

# Fields

  - `PRN_a::Int`: Almanac satellite PRN.
  - `WN_a::Int`: Almanac reference week number (mod 8192).
  - `t_0a::Int`: Almanac reference time of week (seconds).
  - `e::Float64`: Eccentricity (dimensionless).
  - `δi::Float64`: Inclination delta from `i₀ = 0.30 sc` (rad); add the reference.
  - `Ω_dot::Float64`: Rate of right ascension (rad/s).
  - `sqrt_A::Float64`: Square root of the semi-major axis (√m).
  - `Ω_0::Float64`: Longitude of ascending node at weekly epoch (rad).
  - `ω::Float64`: Argument of perigee (rad).
  - `M_0::Float64`: Mean anomaly at reference time (rad).
  - `a_f0::Float64`, `a_f1::Float64`: Clock bias / drift (s, s/s).
  - `l1_health::Bool`, `l2_health::Bool`, `l5_health::Bool`: per-band health.

# Reference

IS-GPS-800J, Figure 3.5-5 / Table 3.5-7.
"""
Base.@kwdef struct GPSL1C_DMidiAlmanac
    PRN_a::Int
    WN_a::Int
    t_0a::Int
    e::Float64
    δi::Float64
    Ω_dot::Float64
    sqrt_A::Float64
    Ω_0::Float64
    ω::Float64
    M_0::Float64
    a_f0::Float64
    a_f1::Float64
    l1_health::Bool
    l2_health::Bool
    l5_health::Bool
end

"""
    GPSL1C_DDifferentialCorrection

One satellite's clock+ephemeris differential-correction packet from subframe 3,
page 5 (IS-GPS-800J Figure 3.5-6 / Figure 3.5-10 / Table 3.5-8).

A page-5 carries the predict/reference times plus exactly one DC packet (a
34-bit CDC segment and a 92-bit EDC segment that form an indivisible pair) for
*another* SV, keyed by `PRN_a`. An all-ones PRN ID (`0xFF` = 255) in any PRN ID
field marks an empty packet — the remainder of the data block is then filler
(IS-GPS-800J §3.5.4.4.4.1). `dc_data_type` selects the data the corrections
apply to: `false` ⇒
CNAV-2 (`D_L1C`), `true` ⇒ legacy NAV (`D`). Semi-circle fields → radians.

# Fields

  - `PRN_a::Int`: PRN the corrections apply to.
  - `t_op_D::Int`: DC data predict time of week (seconds).
  - `t_OD::Int`: Time of DC data (seconds).
  - `dc_data_type::Bool`: false ⇒ CNAV-2, true ⇒ legacy NAV.
  - `δa_f0::Float64`, `δa_f1::Float64`: Clock bias / drift corrections (s, s/s).
  - `UDRA_index::Int`, `UDRA_dot_index::Int`: (rate-of-)UDRA indices (signed).
  - `Δα::Float64`, `Δβ::Float64`: Ephemeris α/β corrections (dimensionless).
  - `Δγ::Float64`: Ephemeris γ correction (rad).
  - `Δi::Float64`, `ΔΩ::Float64`: Inclination / right-ascension corrections (rad).
  - `ΔA::Float64`: Semi-major-axis correction (meters).

# Reference

IS-GPS-800J, Figure 3.5-6 / Figure 3.5-10 / Table 3.5-8.
"""
Base.@kwdef struct GPSL1C_DDifferentialCorrection
    PRN_a::Int
    t_op_D::Int
    t_OD::Int
    dc_data_type::Bool
    δa_f0::Float64
    δa_f1::Float64
    UDRA_index::Int
    UDRA_dot_index::Int
    Δα::Float64
    Δβ::Float64
    Δγ::Float64
    Δi::Float64
    ΔΩ::Float64
    ΔA::Float64
end

"""
    GPSL1C_DData

Decoded GPS L1C-D (CNAV-2) navigation message data.

Holds the subframe-2 clock, ephemeris, and accuracy parameters (IS-GPS-800J
Figure 3.5-1 / Table 3.5-1). Subframe-3 page contents are not parsed in this
slice (issue #39); only the count of CRC-valid subframe-3 pages received is
tracked. Field-naming follows [`GPSL1CAData`](@ref): semi-circle quantities are
converted to radians on decode (multiplied by π), all `Union{Nothing,…}` until
first decoded.

# Sync / timing

  - `toi::Int`: Last validated Time-Of-Interval count (0..399), or `nothing`.
  - `ITOW::Int64`: Interval time of week — number of two-hour epochs since the
    start of the week (subframe 2 bits 14-21).
  - `WN::Int64`: Transmission week number, modulo-8192 (subframe 2 bits 1-13).
  - `t_op::Int64`: CEI data sequence propagation time of week (seconds). Not
    the same field as `t_op_D`, which *is* the DC Data Predict Time of Week.

# Health / accuracy

  - `l1c_health::Bool`: L1C signal health bit (false = OK, true = bad/unavailable).
  - `ura_ed_index::Int64`: Ephemeris URA index (signed).
  - `integrity_status_flag::Bool`: Level of integrity assurance the URA carries
    (bit 566; 3.5.3.10.1). `false` is the legacy level, `true` the enhanced one
    — "URA is integrity assured to the enhanced level only when the integrity
    status flag is '1'" (3.5.3.10).
  - `WN_op::Int64`: CEI data sequence propagation week number (bits 567-574,
    modulo 256). Pairs with `t_op` in the integrity-assured URA_NED expression
    of 3.5.3.10, which needs both.
  - `ura_ned0_index::Int64`, `ura_ned1_index::Int64`, `ura_ned2_index::Int64`:
    Non-elevation-dependent accuracy indices — NED Accuracy, NED Accuracy
    Change, and NED Accuracy Change Rate. Not three interchangeable clock URAs:
    they are the terms of one time-dependent NED bound.

# Ephemeris (Table 3.5-1)

  - `t_0e::Int64`: Ephemeris/clock data reference time of week (seconds; the
    ICD writes `toe`).
  - `ΔA::Float64`: Semi-major axis difference at reference time (meters).
  - `A_dot::Float64`: Change rate in semi-major axis (m/s).
  - `Δn_0::Float64`: Mean motion difference from computed value (rad/s).
  - `Δn_0_dot::Float64`: Rate of mean motion difference (rad/s²).
  - `M_0::Float64`: Mean anomaly at reference time (rad).
  - `e::Float64`: Eccentricity (dimensionless).
  - `ω::Float64`: Argument of perigee (rad).
  - `Ω_0::Float64`: Reference right ascension angle (rad).
  - `i_0::Float64`: Inclination angle at reference time (rad).
  - `ΔΩ_dot::Float64`: Rate of right ascension difference (rad/s).
  - `i_dot::Float64`: Rate of inclination angle (rad/s).
  - `C_is::Float64`, `C_ic::Float64`: Sine/cosine inclination harmonic corrections (rad).
  - `C_rs::Float64`, `C_rc::Float64`: Sine/cosine orbit-radius harmonic corrections (m).
  - `C_us::Float64`, `C_uc::Float64`: Sine/cosine argument-of-latitude harmonic corrections (rad).

# Clock (Table 3.5-1)

  - `t_0c::Int64`: Clock data reference time of week (seconds; the ICD writes
    `toc`). Equals `t_0e` in CNAV-2.
  - `a_f0::Float64`, `a_f1::Float64`, `a_f2::Float64`: Clock bias / drift / drift-rate.
  - `T_GD::Float64`: L1/L2 P(Y) inter-signal correction (seconds).
  - `ISC_L1CP::Float64`, `ISC_L1CD::Float64`: L1CP / L1CD inter-signal corrections (seconds).

# Subframe 3 (IS-GPS-800J §3.5.4 — IRN-IS-800J layout)

Subframe-3 pages are parsed after their CRC passes, dispatching on the 6-bit
page number (bits 9-14; bits 1-8 are the transmitting PRN). Layouts follow
IS-GPS-800J as amended by IRN-IS-800J-003. `num_sf3_pages_received` counts every
CRC-valid SF3 page regardless of whether its page format is parsed.

## Page 1 — UTC + Klobuchar iono + ISC

  - `A_0UTC,A_1UTC,A_2UTC::Float64`: UTC polynomial (s, s/s, s/s²; the ICD
    names them `A0-n`/`A1-n`/`A2-n`, Table 3.5-3).
  - `Δt_LS,Δt_LSF::Int64`: current/past and future leap-second counts (s).
  - `t_0t::Int64`: UTC reference time of week (s; the ICD writes `t_ot`).
  - `WN_0t,WN_LSF::Int64`: UTC and leap-second reference week numbers (the ICD
    writes `WN_ot` for the first).
  - `DN::Int64`: leap-second reference day number (1-7).
  - `α_0,α_1,α_2,α_3,β_0,β_1,β_2,β_3::Float64`: Klobuchar ionospheric coefficients.
  - `ISC_L1CA,ISC_L2C,ISC_L5I5,ISC_L5Q5::Float64`: inter-signal corrections (s).

## Page 2 — GGTO + EOP

  - `A_0GGTO,A_1GGTO,A_2GGTO::Float64`: GPS/GNSS time-offset polynomial (the ICD
    writes `A0GGTO`/`A1GGTO`/`A2GGTO`).
  - `t_GGTO::Int64`, `WN_GGTO::Int64`: GGTO reference time/week.
  - `GGTO_ID::Int64`: GNSS the time offset refers to — 0 none, 1 Galileo,
    2 GLONASS, 3-7 reserved. (Named "GNSS ID" in IS-GPS-800 ≤ Rev J; renamed
    "GGTO ID" by IRN-IS-800J-003.)
  - `t_EOP::Int64`: EOP reference time of week (s).
  - `PM_X,PM_X_dot,PM_Y,PM_Y_dot::Float64`: polar-motion values/rates.
  - `ΔUT_GPS,ΔUT_GPS_dot::Float64`: UT1-GPS (UT1−GPST) difference and rate.

## Pages 3/4/5 — keyed dictionaries (`nothing` until first decoded)

  - `reduced_almanacs::Dictionary{Int,GPSL1C_DReducedAlmanac}` (page 3).
  - `midi_almanacs::Dictionary{Int,GPSL1C_DMidiAlmanac}` (page 4).
  - `differential_corrections::Dictionary{Int,GPSL1C_DDifferentialCorrection}` (page 5).

## Page 6 — Text

  - `text_message::String`: 29 ASCII characters (control chars stripped).

## Counters

  - `num_sf3_pages_received::Int`: Count of CRC-valid subframe-3 pages received.

# Reference

IS-GPS-800J, Figures 3.5-1 through 3.5-9 and Tables 3.5-1, 3.5-3 … 3.5-8.
"""
Base.@kwdef struct GPSL1C_DData <: AbstractGPSData
    toi::Union{Nothing,Int} = nothing
    ITOW::Union{Nothing,Int64} = nothing
    WN::Union{Nothing,Int64} = nothing
    t_op::Union{Nothing,Int64} = nothing

    l1c_health::Union{Nothing,Bool} = nothing
    ura_ed_index::Union{Nothing,Int64} = nothing
    integrity_status_flag::Union{Nothing,Bool} = nothing
    WN_op::Union{Nothing,Int64} = nothing
    ura_ned0_index::Union{Nothing,Int64} = nothing
    ura_ned1_index::Union{Nothing,Int64} = nothing
    ura_ned2_index::Union{Nothing,Int64} = nothing

    t_0e::Union{Nothing,Int64} = nothing
    ΔA::Union{Nothing,Float64} = nothing
    A_dot::Union{Nothing,Float64} = nothing
    Δn_0::Union{Nothing,Float64} = nothing
    Δn_0_dot::Union{Nothing,Float64} = nothing
    M_0::Union{Nothing,Float64} = nothing
    e::Union{Nothing,Float64} = nothing
    ω::Union{Nothing,Float64} = nothing
    Ω_0::Union{Nothing,Float64} = nothing
    i_0::Union{Nothing,Float64} = nothing
    ΔΩ_dot::Union{Nothing,Float64} = nothing
    i_dot::Union{Nothing,Float64} = nothing
    C_is::Union{Nothing,Float64} = nothing
    C_ic::Union{Nothing,Float64} = nothing
    C_rs::Union{Nothing,Float64} = nothing
    C_rc::Union{Nothing,Float64} = nothing
    C_us::Union{Nothing,Float64} = nothing
    C_uc::Union{Nothing,Float64} = nothing

    t_0c::Union{Nothing,Int64} = nothing
    a_f0::Union{Nothing,Float64} = nothing
    a_f1::Union{Nothing,Float64} = nothing
    a_f2::Union{Nothing,Float64} = nothing
    T_GD::Union{Nothing,Float64} = nothing
    ISC_L1CP::Union{Nothing,Float64} = nothing
    ISC_L1CD::Union{Nothing,Float64} = nothing

    # --- Subframe 3, page 1: UTC + iono + ISC (IS-GPS-800J Fig 3.5-2) ---
    A_0UTC::Union{Nothing,Float64} = nothing
    A_1UTC::Union{Nothing,Float64} = nothing
    A_2UTC::Union{Nothing,Float64} = nothing
    Δt_LS::Union{Nothing,Int64} = nothing
    t_0t::Union{Nothing,Int64} = nothing
    WN_0t::Union{Nothing,Int64} = nothing
    WN_LSF::Union{Nothing,Int64} = nothing
    DN::Union{Nothing,Int64} = nothing
    Δt_LSF::Union{Nothing,Int64} = nothing
    α_0::Union{Nothing,Float64} = nothing
    α_1::Union{Nothing,Float64} = nothing
    α_2::Union{Nothing,Float64} = nothing
    α_3::Union{Nothing,Float64} = nothing
    β_0::Union{Nothing,Float64} = nothing
    β_1::Union{Nothing,Float64} = nothing
    β_2::Union{Nothing,Float64} = nothing
    β_3::Union{Nothing,Float64} = nothing
    ISC_L1CA::Union{Nothing,Float64} = nothing
    ISC_L2C::Union{Nothing,Float64} = nothing
    ISC_L5I5::Union{Nothing,Float64} = nothing
    ISC_L5Q5::Union{Nothing,Float64} = nothing

    # --- Subframe 3, page 2: GGTO + EOP (IS-GPS-800J Fig 3.5-3) ---
    A_0GGTO::Union{Nothing,Float64} = nothing
    A_1GGTO::Union{Nothing,Float64} = nothing
    A_2GGTO::Union{Nothing,Float64} = nothing
    t_GGTO::Union{Nothing,Int64} = nothing
    WN_GGTO::Union{Nothing,Int64} = nothing
    GGTO_ID::Union{Nothing,Int64} = nothing
    t_EOP::Union{Nothing,Int64} = nothing
    PM_X::Union{Nothing,Float64} = nothing
    PM_X_dot::Union{Nothing,Float64} = nothing
    PM_Y::Union{Nothing,Float64} = nothing
    PM_Y_dot::Union{Nothing,Float64} = nothing
    ΔUT_GPS::Union{Nothing,Float64} = nothing
    ΔUT_GPS_dot::Union{Nothing,Float64} = nothing

    # --- Subframe 3, pages 3/4/5: per-SV keyed dictionaries ---
    reduced_almanacs::Union{Nothing,Dictionary{Int,GPSL1C_DReducedAlmanac}} = nothing
    midi_almanacs::Union{Nothing,Dictionary{Int,GPSL1C_DMidiAlmanac}} = nothing
    differential_corrections::Union{
        Nothing,
        Dictionary{Int,GPSL1C_DDifferentialCorrection},
    } = nothing

    # --- Subframe 3, page 6: text (IS-GPS-800J Fig 3.5-7) ---
    text_message::Union{Nothing,String} = nothing

    num_sf3_pages_received::Int = 0
end

function GPSL1C_DData(
    data::GPSL1C_DData;
    toi = data.toi,
    ITOW = data.ITOW,
    WN = data.WN,
    t_op = data.t_op,
    l1c_health = data.l1c_health,
    ura_ed_index = data.ura_ed_index,
    integrity_status_flag = data.integrity_status_flag,
    WN_op = data.WN_op,
    ura_ned0_index = data.ura_ned0_index,
    ura_ned1_index = data.ura_ned1_index,
    ura_ned2_index = data.ura_ned2_index,
    t_0e = data.t_0e,
    ΔA = data.ΔA,
    A_dot = data.A_dot,
    Δn_0 = data.Δn_0,
    Δn_0_dot = data.Δn_0_dot,
    M_0 = data.M_0,
    e = data.e,
    ω = data.ω,
    Ω_0 = data.Ω_0,
    i_0 = data.i_0,
    ΔΩ_dot = data.ΔΩ_dot,
    i_dot = data.i_dot,
    C_is = data.C_is,
    C_ic = data.C_ic,
    C_rs = data.C_rs,
    C_rc = data.C_rc,
    C_us = data.C_us,
    C_uc = data.C_uc,
    t_0c = data.t_0c,
    a_f0 = data.a_f0,
    a_f1 = data.a_f1,
    a_f2 = data.a_f2,
    T_GD = data.T_GD,
    ISC_L1CP = data.ISC_L1CP,
    ISC_L1CD = data.ISC_L1CD,
    A_0UTC = data.A_0UTC,
    A_1UTC = data.A_1UTC,
    A_2UTC = data.A_2UTC,
    Δt_LS = data.Δt_LS,
    t_0t = data.t_0t,
    WN_0t = data.WN_0t,
    WN_LSF = data.WN_LSF,
    DN = data.DN,
    Δt_LSF = data.Δt_LSF,
    α_0 = data.α_0,
    α_1 = data.α_1,
    α_2 = data.α_2,
    α_3 = data.α_3,
    β_0 = data.β_0,
    β_1 = data.β_1,
    β_2 = data.β_2,
    β_3 = data.β_3,
    ISC_L1CA = data.ISC_L1CA,
    ISC_L2C = data.ISC_L2C,
    ISC_L5I5 = data.ISC_L5I5,
    ISC_L5Q5 = data.ISC_L5Q5,
    A_0GGTO = data.A_0GGTO,
    A_1GGTO = data.A_1GGTO,
    A_2GGTO = data.A_2GGTO,
    t_GGTO = data.t_GGTO,
    WN_GGTO = data.WN_GGTO,
    GGTO_ID = data.GGTO_ID,
    t_EOP = data.t_EOP,
    PM_X = data.PM_X,
    PM_X_dot = data.PM_X_dot,
    PM_Y = data.PM_Y,
    PM_Y_dot = data.PM_Y_dot,
    ΔUT_GPS = data.ΔUT_GPS,
    ΔUT_GPS_dot = data.ΔUT_GPS_dot,
    reduced_almanacs = data.reduced_almanacs,
    midi_almanacs = data.midi_almanacs,
    differential_corrections = data.differential_corrections,
    text_message = data.text_message,
    num_sf3_pages_received = data.num_sf3_pages_received,
)
    GPSL1C_DData(
        toi,
        ITOW,
        WN,
        t_op,
        l1c_health,
        ura_ed_index,
        integrity_status_flag,
        WN_op,
        ura_ned0_index,
        ura_ned1_index,
        ura_ned2_index,
        t_0e,
        ΔA,
        A_dot,
        Δn_0,
        Δn_0_dot,
        M_0,
        e,
        ω,
        Ω_0,
        i_0,
        ΔΩ_dot,
        i_dot,
        C_is,
        C_ic,
        C_rs,
        C_rc,
        C_us,
        C_uc,
        t_0c,
        a_f0,
        a_f1,
        a_f2,
        T_GD,
        ISC_L1CP,
        ISC_L1CD,
        A_0UTC,
        A_1UTC,
        A_2UTC,
        Δt_LS,
        t_0t,
        WN_0t,
        WN_LSF,
        DN,
        Δt_LSF,
        α_0,
        α_1,
        α_2,
        α_3,
        β_0,
        β_1,
        β_2,
        β_3,
        ISC_L1CA,
        ISC_L2C,
        ISC_L5I5,
        ISC_L5Q5,
        A_0GGTO,
        A_1GGTO,
        A_2GGTO,
        t_GGTO,
        WN_GGTO,
        GGTO_ID,
        t_EOP,
        PM_X,
        PM_X_dot,
        PM_Y,
        PM_Y_dot,
        ΔUT_GPS,
        ΔUT_GPS_dot,
        reduced_almanacs,
        midi_almanacs,
        differential_corrections,
        text_message,
        num_sf3_pages_received,
    )
end

# The default struct `==` falls back to `===` (reference equality), which fails
# for the mutable `Dictionary` fields even when their contents match. Compare
# field-by-field (mirrors `GalileoINAVData`).
Base.:(==)(a::GPSL1C_DData, b::GPSL1C_DData) = fields_equal(a, b)

"""
$(TYPEDEF)

Per-decoder cache for the GPS L1C-D signal.

Holds the soft-symbol `CircularDeque{Float32}` (capacity = 1852 = 1800 frame +
52 next-frame BCH segment) and the two Aff3ct LDPC belief-propagation decoders
(subframe 2: K=600, N=1200; subframe 3: K=274, N=548), built lazily from the
committed `.alist` parity matrices in `data/`. The decoder objects are mutable
Aff3ct handles; they are shared by reference through the otherwise-immutable
[`GNSSDecoderState`](@ref).

# Fields

$(TYPEDFIELDS)
"""
struct GPSL1C_DCache <: AbstractGNSSCache
    """
    Soft-symbol buffer (1852 = 1800 frame + 52 next subframe-1)
    """
    soft_buffer::CircularDeque{Float32}
    """
    Aff3ct LDPC BP decoder for subframe 2 (K=600, N=1200)
    """
    sf2_decoder::LDPCScratch
    """
    Aff3ct LDPC BP decoder for subframe 3 (K=274, N=548)
    """
    sf3_decoder::LDPCScratch
    """
    Polarity-resolved 1748-symbol interleaved payload, copied out per frame
    """
    payload::Vector{Float32}
    """
    The same payload deinterleaved: subframe 2 then subframe 3
    """
    deinterleaved::Vector{Float32}
end

function GPSL1C_DCache()
    # The CNAV-2 LDPC codes are systematic with the codeword laid out as
    # [info | parity] (IS-GPS-800J §3.2.3.4: H = [[A B T] [C D E]] acts on
    # [info; parity]); `load_ldpc_decoder` forces that layout onto Aff3ct.
    # Verified against a Spirent GSS post-FEC L1C recording.
    GPSL1C_DCache(
        CircularDeque{Float32}(L1C_D_WINDOW_LENGTH),
        LDPCScratch(alist_path("cnv2_sf2.alist")),
        LDPCScratch(alist_path("cnv2_sf3.alist")),
        Vector{Float32}(undef, L1C_D_PAYLOAD_SYMBOLS),
        Vector{Float32}(undef, L1C_D_PAYLOAD_SYMBOLS),
    )
end

# The LDPC decoder handles are stateless w.r.t. equality (they are runtime
# Aff3ct objects); two L1C-D caches are equal when their soft buffers match.
function Base.:(==)(a::GPSL1C_DCache, b::GPSL1C_DCache)
    deques_equal(a.soft_buffer, b.soft_buffer)
end

function is_subframe2_decoded(data::GPSL1C_DData)
    !isnothing(data.WN) &&
        !isnothing(data.ITOW) &&
        !isnothing(data.t_0e) &&
        !isnothing(data.ΔA) &&
        !isnothing(data.M_0) &&
        !isnothing(data.e) &&
        !isnothing(data.ω) &&
        !isnothing(data.Ω_0) &&
        !isnothing(data.i_0) &&
        !isnothing(data.a_f0) &&
        !isnothing(data.a_f1) &&
        !isnothing(data.l1c_health)
end

function is_decoding_completed_for_positioning(data::GPSL1C_DData)
    !isnothing(data.toi) && is_subframe2_decoded(data)
end

"""
$(TYPEDSIGNATURES)

Create a decoder state for GPS L1C-D (CNAV-2) navigation messages.

Wires up a [`GNSSDecoderState`](@ref) with a 1852-symbol soft-symbol buffer,
the 400-entry BCH(51,8) TOI codeword table (`src/bch_toi.jl`), and two Aff3ct
LDPC belief-propagation decoders loaded lazily from the committed `.alist`
parity matrices in `data/`.

# Arguments

  - `prn::Int`: Pseudo-Random Noise code identifier (1-63 for L1C).

# Returns

  - `GNSSDecoderState{GPSL1C_DData}`: Initialized decoder state for GPS L1C-D.

# Example

```julia
state = GPSL1C_DDecoderState(1)            # PRN 1
state = decode(state, soft_symbols, num_symbols)
if is_sat_healthy(state)
    # Use state.data for positioning
end
```

# See Also

  - [`GNSSDecoderState`](@ref): The underlying state structure
  - [`decode`](@ref): Decode soft symbols using this state
  - [`reset_decoder_state`](@ref): Reset after signal loss
"""
function GPSL1C_DDecoderState(prn)
    GNSSDecoderState(
        prn,
        GPSL1C_DData(),
        GPSL1C_DData(),
        GPSL1C_DConstants(),
        GPSL1C_DCache(),
        nothing,
        false,
    )
end

function GNSSDecoderState(system::GPSL1C_D, prn)
    GPSL1C_DDecoderState(prn)
end

# The signal this decoder demodulates; signal metadata is forwarded through it
# (see `src/gps/l1ca.jl`). CNAV-2 rides on the L1C data component — `GPSL1C_D` —
# while `GPSL1C_P` is the dataless pilot.
get_signal_type(::GPSL1C_DConstants) = GPSL1C_D

"""
$(TYPEDSIGNATURES)

Reset the GPS L1C-D decoder state after a signal loss or reacquisition.

Clears the in-flight sync state (soft-symbol buffer and TOI) and the validated
data, while preserving the long-lived CED/clock fields in `raw_data` so a
[`GNSSReceiver`] can re-use the satellite after reacquisition without re-decoding
all of subframe 2. Mirrors the semantics of the GPS L1 C/A and Galileo E1B
implementations.

# Arguments

  - `state::GNSSDecoderState{<:GPSL1C_DData}`: Current GPS L1C-D decoder state.

# Returns

  - `GNSSDecoderState{<:GPSL1C_DData}`: Reset decoder state with cleared buffers.

# See Also

  - [`GPSL1C_DDecoderState`](@ref): Create a fresh decoder state
  - [`decode`](@ref): Continue decoding after reset
"""
function reset_decoder_state(state::GNSSDecoderState{<:GPSL1C_DData})
    empty!(state.cache.soft_buffer)
    GNSSDecoderState(
        state;
        raw_data = GPSL1C_DData(state.raw_data; toi = nothing),
        data = GPSL1C_DData(),
        num_bits_after_valid_syncro_sequence = nothing,
        is_shifted_by_180_degrees = false,
    )
end

# ---- Sync ------------------------------------------------------------------
#
# Override the generic packed-buffer sync. L1C-D sync runs the BCH(51,8) TOI
# match directly on the soft-symbol deque: the first 52 symbols must equal the
# codeword for some `toi`, and the last 52 the codeword for `(toi+1) mod 400`,
# in either polarity. `sync_bch_toi` (src/bch_toi.jl) implements exactly this.

"""
    try_sync(state::GNSSDecoderState{<:GPSL1C_DData}) -> Union{Nothing,BCHToiSync}

Per-signal sync hook for GPS L1C-D. Hard-slices the leading 52 and trailing 52
soft symbols of the 1852-symbol window and runs the multi-subframe BCH(51,8)
match (`sync_bch_toi`). Returns the [`BCHToiSync`](@ref) on a match (carrying
the TOI and the detected polarity flip) or `nothing`.
"""
function try_sync(state::GNSSDecoderState{<:GPSL1C_DData})
    deque = soft_buffer(state)
    # `pack_bits_lsb_first` (src/gnss.jl) is `soft_to_hard_codeword` read straight
    # off the deque: this runs once per symbol, so materialising two 52-element
    # slices here would allocate on every one of them.
    first52 = pack_bits_lsb_first(deque, 1, L1C_D_SUBFRAME1_LENGTH)
    next52 = pack_bits_lsb_first(deque, L1C_D_FRAME_LENGTH + 1, L1C_D_SUBFRAME1_LENGTH)
    sync_bch_toi(first52, next52)
end

"""
    complement_buffer_if_necessary(state::GNSSDecoderState{<:GPSL1C_DData}, sync)

Resolve the 180° polarity ambiguity for L1C-D. `sync::BCHToiSync` already
carries the detected polarity from `sync_bch_toi`; record it on the state and
pass the sync object through unchanged for `decode_syncro_sequence`.
"""
function complement_buffer_if_necessary(
    state::GNSSDecoderState{<:GPSL1C_DData},
    sync::BCHToiSync,
)
    GNSSDecoderState(state; is_shifted_by_180_degrees = sync.polarity_flipped), sync
end

# ---- Subframe pipeline -----------------------------------------------------

"""
    decode_syncro_sequence(state::GNSSDecoderState{<:GPSL1C_DData}, sync::BCHToiSync)

Process one locked CNAV-2 frame. Enforces a monotonic TOI (each frame's TOI must
be `(previous + 1) mod 400`; a break resets the decoder to searching), then
deinterleaves and LDPC-decodes subframes 2 and 3 from the 1748-symbol payload,
CRC-checking each. Subframe 2 fields are parsed into `raw_data`; subframe 3 is
recorded as a received page (issue #39 will parse it).
"""
function decode_syncro_sequence(state::GNSSDecoderState{<:GPSL1C_DData}, sync::BCHToiSync)
    # Monotonic-TOI validation: a locked frame whose TOI does not follow the
    # previous one by exactly +1 (mod 400) indicates loss of frame lock.
    prev_toi = state.raw_data.toi
    if !isnothing(prev_toi) && (prev_toi + 1) % TOI_RANGE != sync.toi
        return reset_decoder_state(state)
    end
    state = GNSSDecoderState(state; raw_data = GPSL1C_DData(state.raw_data; toi = sync.toi))

    # Extract the 1748-symbol interleaved SF2+SF3 payload (symbols 53..1800),
    # applying the polarity flip by negating soft symbols up front.
    interleaved = copy_soft_window!(
        state.cache.payload,
        soft_buffer(state),
        L1C_D_SUBFRAME1_LENGTH,
        L1C_D_PAYLOAD_SYMBOLS,
        state.is_shifted_by_180_degrees,
    )
    deinterleaved = deinterleave!(
        state.cache.deinterleaved,
        interleaved,
        L1C_D_INTERLEAVER_ROWS,
        L1C_D_INTERLEAVER_COLS,
    )
    sf2_symbols = @view deinterleaved[1:L1C_D_SF2_SYMBOLS]
    sf3_symbols = @view deinterleaved[(L1C_D_SF2_SYMBOLS+1):L1C_D_PAYLOAD_SYMBOLS]

    state = decode_subframe2(state, sf2_symbols)
    state = decode_subframe3(state, sf3_symbols)
    return state
end

# ---- Subframe 2 bit-field extraction (IS-GPS-800J Figure 3.5-1) ------------
#
# `word` is the 600-bit subframe-2 info block packed MSB-first into a `UInt600`
# (bit 1 = MSB). Fields are read by 1-based start bit and length through the
# shared `get_bits` / `get_twos_complement_num` / `get_bit` helpers.

function decode_subframe2(state::GNSSDecoderState{<:GPSL1C_DData}, sf2_symbols)
    word = ldpc_decode_word(state.cache.sf2_decoder, sf2_symbols, UInt600)
    isnothing(word) && return state  # silently drop on CRC failure
    word_length = L1C_D_SF2_INFO_BITS

    PI = state.constants.PI

    WN = Int(get_bits(word, word_length, 1, 13))
    ITOW = Int(get_bits(word, word_length, 14, 8))
    t_op = Int(get_bits(word, word_length, 22, 11)) * 300
    l1c_health = get_bit(word, word_length, 33)
    ura_ed_index = get_twos_complement_num(word, word_length, 34, 5)
    t_0e = Int(get_bits(word, word_length, 39, 11)) * 300
    ΔA = get_twos_complement_num(word, word_length, 50, 26) * 2.0^-9
    A_dot = get_twos_complement_num(word, word_length, 76, 25) * 2.0^-21
    Δn_0 = get_twos_complement_num(word, word_length, 101, 17) * 2.0^-44 * PI
    Δn_0_dot = get_twos_complement_num(word, word_length, 118, 23) * 2.0^-57 * PI
    M_0 = get_twos_complement_num(word, word_length, 141, 33) * 2.0^-32 * PI
    # 33-bit field: convert through Int64 explicitly (a platform Int would
    # overflow on 32-bit), matching the BeiDou B-CNAV parsers.
    e = Int64(get_bits(word, word_length, 174, 33)) * 2.0^-34
    ω = get_twos_complement_num(word, word_length, 207, 33) * 2.0^-32 * PI
    Ω_0 = get_twos_complement_num(word, word_length, 240, 33) * 2.0^-32 * PI
    i_0 = get_twos_complement_num(word, word_length, 273, 33) * 2.0^-32 * PI
    ΔΩ_dot = get_twos_complement_num(word, word_length, 306, 17) * 2.0^-44 * PI
    i_dot = get_twos_complement_num(word, word_length, 323, 15) * 2.0^-44 * PI
    C_is = get_twos_complement_num(word, word_length, 338, 16) * 2.0^-30
    C_ic = get_twos_complement_num(word, word_length, 354, 16) * 2.0^-30
    C_rs = get_twos_complement_num(word, word_length, 370, 24) * 2.0^-8
    C_rc = get_twos_complement_num(word, word_length, 394, 24) * 2.0^-8
    C_us = get_twos_complement_num(word, word_length, 418, 21) * 2.0^-30
    C_uc = get_twos_complement_num(word, word_length, 439, 21) * 2.0^-30
    ura_ned0_index = get_twos_complement_num(word, word_length, 460, 5)
    ura_ned1_index = Int(get_bits(word, word_length, 465, 3))
    ura_ned2_index = Int(get_bits(word, word_length, 468, 3))
    a_f0 = get_twos_complement_num(word, word_length, 471, 26) * 2.0^-35
    a_f1 = get_twos_complement_num(word, word_length, 497, 20) * 2.0^-48
    a_f2 = get_twos_complement_num(word, word_length, 517, 10) * 2.0^-60
    T_GD = get_twos_complement_num(word, word_length, 527, 13) * 2.0^-35
    ISC_L1CP = get_twos_complement_num(word, word_length, 540, 13) * 2.0^-35
    ISC_L1CD = get_twos_complement_num(word, word_length, 553, 13) * 2.0^-35
    # Bits 566 and 567-574 close subframe 2 before its 2 reserved bits and the
    # CRC (Figure 3.5-1). Figure 3.5-1 drops `ISF - 1 BIT` and `RESERVED - 2
    # BITS` onto a second label row below the main one, which is how they came
    # to be skipped when this parser was first transcribed from the figure.
    integrity_status_flag = get_bit(word, word_length, 566)
    WN_op = Int(get_bits(word, word_length, 567, 8))

    raw = GPSL1C_DData(
        state.raw_data;
        WN,
        ITOW,
        t_op,
        l1c_health,
        ura_ed_index,
        integrity_status_flag,
        WN_op,
        t_0e,
        ΔA,
        A_dot,
        Δn_0,
        Δn_0_dot,
        M_0,
        e,
        ω,
        Ω_0,
        i_0,
        ΔΩ_dot,
        i_dot,
        C_is,
        C_ic,
        C_rs,
        C_rc,
        C_us,
        C_uc,
        ura_ned0_index,
        ura_ned1_index,
        ura_ned2_index,
        t_0c = t_0e,  # CNAV-2 shares one reference time for clock and ephemeris
        a_f0,
        a_f1,
        a_f2,
        T_GD,
        ISC_L1CP,
        ISC_L1CD,
    )
    GNSSDecoderState(state; raw_data = raw)
end

# ---- Subframe 3 page parsing (IS-GPS-800J §3.5.4) --------------------------
#
# Every SF3 page is a 274-bit info block: bits 1-8 are the transmitting PRN,
# bits 9-14 the 6-bit page number, then page-specific fields, with the trailing
# 24 bits a CRC-24Q. After the CRC passes the 274 bits are packed MSB-first into
# a `UInt288` (`get_bits(word, 274, …)` addresses the right-aligned 274 logical
# bits); we dispatch on the page number and merge the parsed fields into
# `raw_data` immutably (same style as SF2). Layouts follow IS-GPS-800J as
# amended by IRN-IS-800J-003 — which for page 1 only restores Greek letters the
# base PDF rendered as question marks; the four ISC fields at bits 177/190/203/216
# are in base Rev J already (Figure 3.5-2, Table 6.2-18).

function decode_subframe3(state::GNSSDecoderState{<:GPSL1C_DData}, sf3_symbols)
    word = ldpc_decode_word(state.cache.sf3_decoder, sf3_symbols, UInt288)
    isnothing(word) && return state  # silently drop on CRC failure

    # CRC-valid page received: count it, then dispatch on the page number.
    raw = GPSL1C_DData(
        state.raw_data;
        num_sf3_pages_received = state.raw_data.num_sf3_pages_received + 1,
    )

    word_length = L1C_D_SF3_INFO_BITS
    page = Int(get_bits(word, word_length, 9, 6))  # bits 1-8 PRN, bits 9-14 page (IS-GPS-800J §3.5.4)
    raw = if page == L1C_D_SF3_PAGE_UTC_IONO
        parse_sf3_page1(raw, word, state.constants.PI)
    elseif page == L1C_D_SF3_PAGE_GGTO_EOP
        parse_sf3_page2(raw, word, state.constants.PI)
    elseif page == L1C_D_SF3_PAGE_REDUCED_ALMANAC
        parse_sf3_page3(raw, word, state.constants.PI)
    elseif page == L1C_D_SF3_PAGE_MIDI_ALMANAC
        parse_sf3_page4(raw, word, state.constants.PI)
    elseif page == L1C_D_SF3_PAGE_DIFF_CORRECTION
        parse_sf3_page5(raw, word, state.constants.PI)
    elseif page == L1C_D_SF3_PAGE_TEXT
        parse_sf3_page6(raw, word)
    else
        raw  # unsupported/reserved page (e.g. 7 SV-config, 8 ISM): counted, ignored
    end

    GNSSDecoderState(state; raw_data = raw)
end

"""
Subframe 3, page 1 — UTC + Klobuchar iono + ISC (IS-GPS-800J Fig 3.5-2, Table 3.5-3).
"""
function parse_sf3_page1(raw::GPSL1C_DData, word::UInt288, PI::Float64)
    word_length = L1C_D_SF3_INFO_BITS
    GPSL1C_DData(
        raw;
        # UTC polynomial (Table 3.5-3).
        A_0UTC = get_twos_complement_num(word, word_length, 15, 16) * 2.0^-35,
        A_1UTC = get_twos_complement_num(word, word_length, 31, 13) * 2.0^-51,
        A_2UTC = get_twos_complement_num(word, word_length, 44, 7) * 2.0^-68,
        Δt_LS = get_twos_complement_num(word, word_length, 51, 8),
        t_0t = Int(get_bits(word, word_length, 59, 16)) * 2^4,
        WN_0t = Int(get_bits(word, word_length, 75, 13)),
        WN_LSF = Int(get_bits(word, word_length, 88, 13)),
        DN = Int(get_bits(word, word_length, 101, 4)),
        Δt_LSF = get_twos_complement_num(word, word_length, 105, 8),
        # Klobuchar ionospheric coefficients (IS-GPS-200 Table 20-X; all 8-bit
        # two's-complement, scaled in seconds / seconds-per-semicircle^n).
        α_0 = get_twos_complement_num(word, word_length, 113, 8) * 2.0^-30,
        α_1 = get_twos_complement_num(word, word_length, 121, 8) * 2.0^-27,
        α_2 = get_twos_complement_num(word, word_length, 129, 8) * 2.0^-24,
        α_3 = get_twos_complement_num(word, word_length, 137, 8) * 2.0^-24,
        β_0 = get_twos_complement_num(word, word_length, 145, 8) * 2.0^11,
        β_1 = get_twos_complement_num(word, word_length, 153, 8) * 2.0^14,
        β_2 = get_twos_complement_num(word, word_length, 161, 8) * 2.0^16,
        β_3 = get_twos_complement_num(word, word_length, 169, 8) * 2.0^16,
        # Inter-signal corrections (Fig 3.5-2; 13-bit two's complement, 2^-35 s).
        ISC_L1CA = get_twos_complement_num(word, word_length, 177, 13) * 2.0^-35,
        ISC_L2C = get_twos_complement_num(word, word_length, 190, 13) * 2.0^-35,
        ISC_L5I5 = get_twos_complement_num(word, word_length, 203, 13) * 2.0^-35,
        ISC_L5Q5 = get_twos_complement_num(word, word_length, 216, 13) * 2.0^-35,
    )
end

"""
Subframe 3, page 2 — GGTO + EOP (IS-GPS-800J Fig 3.5-3, Tables 3.5-4/3.5-5).
"""
function parse_sf3_page2(raw::GPSL1C_DData, word::UInt288, PI::Float64)
    word_length = L1C_D_SF3_INFO_BITS
    GPSL1C_DData(
        raw;
        # GGTO (Table 3.5-4). Field order in Fig 3.5-3: tGGTO, WNGGTO, A0, A1, A2.
        t_GGTO = Int(get_bits(word, word_length, 18, 16)) * 2^4,
        WN_GGTO = Int(get_bits(word, word_length, 34, 13)),
        A_0GGTO = get_twos_complement_num(word, word_length, 47, 16) * 2.0^-35,
        A_1GGTO = get_twos_complement_num(word, word_length, 63, 13) * 2.0^-51,
        A_2GGTO = get_twos_complement_num(word, word_length, 76, 7) * 2.0^-68,
        GGTO_ID = Int(get_bits(word, word_length, 15, 3)),
        # EOP (Table 3.5-5). All fields are contiguous in the info block; Figure
        # 3.5-3 only *draws* PM_X across its 100-bit row boundary — the 2 MSBs end
        # row 1 (bits 99-100) and the 19 LSBs begin row 2 (bits 101-119) — so the
        # plain 21-bit read at bit 99 is correct (likewise ΔUT_GPS spans the next row).
        t_EOP = Int(get_bits(word, word_length, 83, 16)) * 2^4,
        PM_X = get_twos_complement_num(word, word_length, 99, 21) * 2.0^-20,
        PM_X_dot = get_twos_complement_num(word, word_length, 120, 15) * 2.0^-21,
        PM_Y = get_twos_complement_num(word, word_length, 135, 21) * 2.0^-20,
        PM_Y_dot = get_twos_complement_num(word, word_length, 156, 15) * 2.0^-21,
        ΔUT_GPS = get_twos_complement_num(word, word_length, 171, 31) * 2.0^-23,
        ΔUT_GPS_dot = get_twos_complement_num(word, word_length, 202, 19) * 2.0^-25,
    )
end

"""
Decode one 33-bit reduced-almanac packet starting at 1-based bit `start` (IS-GPS-800J Fig 3.5-9).
"""
function _reduced_almanac_packet(
    word::UInt288,
    start::Int,
    WN_a::Int,
    t_0a::Int,
    PI::Float64,
)
    word_length = L1C_D_SF3_INFO_BITS
    PRN_a = Int(get_bits(word, word_length, start, 8))
    # PRNa == 0 marks an empty packet; per IS-GPS-800J §3.5.4.3.5.1.1 all
    # subsequent bits through the last packet are then filler, so the caller
    # stops here rather than parsing the remaining packets.
    PRN_a == 0 && return nothing
    GPSL1C_DReducedAlmanac(;
        PRN_a,
        WN_a,
        t_0a,
        δA = get_twos_complement_num(word, word_length, start + 8, 8) * 2.0^9,
        Ω_0 = get_twos_complement_num(word, word_length, start + 16, 7) * 2.0^-6 * PI,
        Φ_0 = get_twos_complement_num(word, word_length, start + 23, 7) * 2.0^-6 * PI,
        l1_health = get_bit(word, word_length, start + 30),
        l2_health = get_bit(word, word_length, start + 31),
        l5_health = get_bit(word, word_length, start + 32),
    )
end

"""
Subframe 3, page 3 — six reduced-almanac packets (IS-GPS-800J Fig 3.5-4).
"""
function parse_sf3_page3(raw::GPSL1C_DData, word::UInt288, PI::Float64)
    word_length = L1C_D_SF3_INFO_BITS
    WN_a = Int(get_bits(word, word_length, 15, 13))
    t_0a = Int(get_bits(word, word_length, 28, 8)) * 2^12
    almanacs = raw.reduced_almanacs
    # Six 33-bit packets at bits 36, 69, 102, 135, 168, 201.
    for start in (36, 69, 102, 135, 168, 201)
        packet = _reduced_almanac_packet(word, start, WN_a, t_0a, PI)
        isnothing(packet) && break  # PRNa==0 ⇒ rest of page is filler (§3.5.4.3.5.1.1)
        almanacs = _merge_keyed(almanacs, packet.PRN_a, packet)
    end
    GPSL1C_DData(raw; reduced_almanacs = almanacs)
end

"""
Subframe 3, page 4 — one Midi almanac (IS-GPS-800J Fig 3.5-5, Table 3.5-7).
"""
function parse_sf3_page4(raw::GPSL1C_DData, word::UInt288, PI::Float64)
    word_length = L1C_D_SF3_INFO_BITS
    PRN_a = Int(get_bits(word, word_length, 36, 8))
    PRN_a == 0 && return raw  # empty almanac
    alm = GPSL1C_DMidiAlmanac(;
        PRN_a,
        WN_a = Int(get_bits(word, word_length, 15, 13)),
        t_0a = Int(get_bits(word, word_length, 28, 8)) * 2^12,
        l1_health = get_bit(word, word_length, 44),
        l2_health = get_bit(word, word_length, 45),
        l5_health = get_bit(word, word_length, 46),
        e = Int(get_bits(word, word_length, 47, 11)) * 2.0^-16,
        δi = get_twos_complement_num(word, word_length, 58, 11) * 2.0^-14 * PI,
        Ω_dot = get_twos_complement_num(word, word_length, 69, 11) * 2.0^-33 * PI,
        sqrt_A = Int(get_bits(word, word_length, 80, 17)) * 2.0^-4,
        Ω_0 = get_twos_complement_num(word, word_length, 97, 16) * 2.0^-15 * PI,
        ω = get_twos_complement_num(word, word_length, 113, 16) * 2.0^-15 * PI,
        M_0 = get_twos_complement_num(word, word_length, 129, 16) * 2.0^-15 * PI,
        a_f0 = get_twos_complement_num(word, word_length, 145, 11) * 2.0^-20,
        a_f1 = get_twos_complement_num(word, word_length, 156, 10) * 2.0^-37,
    )
    GPSL1C_DData(raw; midi_almanacs = _merge_keyed(raw.midi_almanacs, PRN_a, alm))
end

"""
Subframe 3, page 5 — one differential-correction packet (IS-GPS-800J Fig 3.5-6/3.5-10, Table 3.5-8).
"""
function parse_sf3_page5(raw::GPSL1C_DData, word::UInt288, PI::Float64)
    word_length = L1C_D_SF3_INFO_BITS
    # Page-level fields precede the 126-bit CDC+EDC packet. Layout (Fig 3.5-6):
    # bit 15 t_op-D (11, scale 300), bit 26 t_OD (11, scale 300),
    # bit 37 DC data type (1), then the CDC segment (bit 38) and EDC segment.
    t_op_D = Int(get_bits(word, word_length, 15, 11)) * 300
    t_OD = Int(get_bits(word, word_length, 26, 11)) * 300
    dc_data_type = get_bit(word, word_length, 37)
    # CDC segment (Fig 3.5-10): PRN ID(8) δaf0(13) δaf1(8) UDRA(5) — starts bit 38.
    cdc_start_bit = 38
    PRN_a = Int(get_bits(word, word_length, cdc_start_bit, 8))
    PRN_a == 0xff && return raw  # all-ones PRN ⇒ no DC data in this packet
    δa_f0 = get_twos_complement_num(word, word_length, cdc_start_bit + 8, 13) * 2.0^-35
    δa_f1 = get_twos_complement_num(word, word_length, cdc_start_bit + 21, 8) * 2.0^-51
    UDRA_index = get_twos_complement_num(word, word_length, cdc_start_bit + 29, 5)
    # EDC segment (Fig 3.5-10): PRN ID(8) Δα(14) Δβ(14) Δγ(15) Δi(12) ΔΩ(12)
    # ΔA(12) UDRA-dot(5) — starts at bit 72 (38 + 34).
    edc_start_bit = 72
    Δα = get_twos_complement_num(word, word_length, edc_start_bit + 8, 14) * 2.0^-34
    Δβ = get_twos_complement_num(word, word_length, edc_start_bit + 22, 14) * 2.0^-34
    Δγ = get_twos_complement_num(word, word_length, edc_start_bit + 36, 15) * 2.0^-32 * PI
    Δi = get_twos_complement_num(word, word_length, edc_start_bit + 51, 12) * 2.0^-32 * PI
    ΔΩ = get_twos_complement_num(word, word_length, edc_start_bit + 63, 12) * 2.0^-32 * PI
    ΔA = get_twos_complement_num(word, word_length, edc_start_bit + 75, 12) * 2.0^-9
    UDRA_dot_index = get_twos_complement_num(word, word_length, edc_start_bit + 87, 5)
    differential_correction = GPSL1C_DDifferentialCorrection(;
        PRN_a,
        t_op_D,
        t_OD,
        dc_data_type,
        δa_f0,
        δa_f1,
        UDRA_index,
        UDRA_dot_index,
        Δα,
        Δβ,
        Δγ,
        Δi,
        ΔΩ,
        ΔA,
    )
    GPSL1C_DData(
        raw;
        differential_corrections = _merge_keyed(
            raw.differential_corrections,
            PRN_a,
            differential_correction,
        ),
    )
end

"""
Subframe 3, page 6 — 29 ASCII characters at bits 19-250 (IS-GPS-800J Fig 3.5-7).
"""
function parse_sf3_page6(raw::GPSL1C_DData, word::UInt288)
    word_length = L1C_D_SF3_INFO_BITS
    chars = Char[]
    for k = 0:28
        code = Int(get_bits(word, word_length, 19 + 8k, 8))
        # Keep printable ASCII; skip NUL/control padding so the message is clean.
        (code >= 0x20 && code < 0x7f) && push!(chars, Char(code))
    end
    GPSL1C_DData(raw; text_message = String(chars))
end

"""
    validate_data(state::GNSSDecoderState{<:GPSL1C_DData})

Promote `raw_data` to `data` once subframe 2 has been fully decoded. The
per-frame BCH re-check and monotonic-TOI enforcement are performed inline by
[`decode_syncro_sequence`](@ref); this hook only publishes validated data and
arms the streaming counter.
"""
function validate_data(state::GNSSDecoderState{<:GPSL1C_DData})
    if is_decoding_completed_for_positioning(state.raw_data)
        return GNSSDecoderState(
            state;
            data = state.raw_data,
            num_bits_after_valid_syncro_sequence = state.constants.preamble_length,
        )
    end
    return state
end

"""
$(TYPEDSIGNATURES)

Check if the GPS L1C-D satellite is healthy and usable for positioning.

Examines the 1-bit L1C signal health flag from subframe 2 (IS-GPS-800J
§3.5.3.4): a satellite is healthy iff the health bit is 0 (Signal OK).

!!! warning

    Requires subframe 2 to have been decoded; returns `false` until then.

# Arguments

  - `state::GNSSDecoderState{<:GPSL1C_DData}`: GPS L1C-D decoder state.

# Returns

  - `Bool`: `true` iff the L1C signal-health bit indicates OK.
"""
function is_sat_healthy(state::GNSSDecoderState{<:GPSL1C_DData})
    state.data.l1c_health === false
end
