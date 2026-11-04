# GPS L1C-D (CNAV-2) decoder — TOI sync, subframe 2 and subframe 3.
#
# L1C-D carries CNAV-2 as 18-second frames of 1800 symbols at 100 sps:
#
#   - Subframe 1: 52 symbols, a BCH(51,8)+parity encoding of a 9-bit TOI count
#     (0..399), used only for frame sync (see `src/coding/bch_toi.jl`).
#   - Subframes 2 + 3: 1748 symbols, block-interleaved (38×46, see
#     `src/coding/deinterleave.jl`): a 1200-symbol rate-½ LDPC codeword for
#     subframe 2 (600 info bits), then a 548-symbol one for subframe 3 (274 bits).
#
# Subframe 2 carries clock and ephemeris (IS-GPS-800J Figure 3.5-1 / Table 3.5-1).
# Subframe 3 is paged: pages 1-6 (UTC + ionosphere + ISC, GGTO + EOP, reduced almanac,
# midi almanac, differential correction, text; Figures 3.5-2 to 3.5-7) and page 8 (the
# ISM Packet, IRN-IS-800J-004 Figure 3.5-8a) are parsed; other pages are counted only.
#
# The soft-symbol `CircularDeque` holds one frame (1800) plus the 52-symbol subframe-1
# segment of the *next* frame (1852), so the BCH match is confirmed at both ends of
# the window before locking (generic streaming framework, `src/gnss.jl`).

# `import Aff3ct` for the qualified `Aff3ct.decode` call (clashes with
# `GNSSDecoder.decode`); only the LDPC constructors are brought in by name.
import Aff3ct
using Aff3ct: LDPCMatrix, LDPCBPDecoder

const L1C_D_FRAME_LENGTH = 1800
const L1C_D_SUBFRAME1_LENGTH = 52
# Total symbols accumulated before a sync attempt: one frame plus the next subframe 1.
const L1C_D_WINDOW_LENGTH = L1C_D_FRAME_LENGTH + L1C_D_SUBFRAME1_LENGTH  # 1852

# Time of week (IS-GPS-800J §3.2.3.1, Figure 3.5-1): 84 two-hour intervals counted by
# subframe-2 `ITOW`, each of 400 18-second frames counted by subframe-1 `TOI`.
const L1C_D_ITOW_SECONDS = 7200
const L1C_D_ITOW_RANGE = 84  # two-hour intervals per week
const L1C_D_FRAME_SECONDS = 18

# Subframe 2/3 channel-symbol counts (IS-GPS-800J §3.2.3).
const L1C_D_SF2_SYMBOLS = 1200
const L1C_D_SF3_SYMBOLS = 548
const L1C_D_PAYLOAD_SYMBOLS = L1C_D_SF2_SYMBOLS + L1C_D_SF3_SYMBOLS  # 1748
const L1C_D_INTERLEAVER_ROWS = 38
const L1C_D_INTERLEAVER_COLS = 46  # 38 * 46 == 1748

# LDPC info-block sizes (rate ½): N channel symbols -> K info bits.
const L1C_D_SF2_INFO_BITS = 600
const L1C_D_SF3_INFO_BITS = 274

# Bit-field extraction: after LDPC + CRC the info bits are packed MSB-first into a wide
# integer and parsed with the shared `get_bits`/`get_twos_complement_num`/`get_bit`
# helpers (src/bit_fiddling.jl), as in the GPS L1 C/A and Galileo E1B decoders (#48).
#
# SF2's 600 bits go in a `UInt600`, SF3's 274 in a `UInt288` (both defined in
# `gnss.jl`, shared with BeiDou B1C): `get_bits(word, 274, …)` addresses the 274
# logical bits as long as they are right-aligned.

# Subframe-3 page numbers (IS-GPS-800J §3.5.4, figure-by-figure). The 6-bit
# page-number field lives in bits 9-14 of every SF3 page (bits 1-8 are PRN).
const L1C_D_SF3_PAGE_UTC_IONO = 1   # Figure 3.5-2: UTC + Klobuchar iono + ISC
const L1C_D_SF3_PAGE_GGTO_EOP = 2   # Figure 3.5-3: GGTO + Earth orientation
const L1C_D_SF3_PAGE_REDUCED_ALMANAC = 3  # Figure 3.5-4: 6 reduced-almanac packets
const L1C_D_SF3_PAGE_MIDI_ALMANAC = 4     # Figure 3.5-5: one Midi almanac
const L1C_D_SF3_PAGE_DIFF_CORRECTION = 5  # Figure 3.5-6: clock+ephemeris DC
const L1C_D_SF3_PAGE_TEXT = 6             # Figure 3.5-7: 29 ASCII characters
const L1C_D_SF3_PAGE_ISM = 8  # IRN-IS-800J-004 Figure 3.5-8a: ISM Packet at bits 15-250

# Characters in a page-6 text message (IS-GPS-800J Figure 3.5-7).
const L1C_D_TEXT_LENGTH = 29

# Slots of the PRN-keyed subframe-3 stores (reduced/Midi almanacs, differential
# corrections): `PRN_a` and `PRN_ID` are 8 bits wide, so one slot per value 0..255 covers
# any PRN.
const L1C_D_PRN_SLOTS = 256

"""
    GPSL1C_DConstants

CNAV-2 message structure parameters for GPS L1C-D decoding.

In the generic streaming framework, `preamble_length` is the 52-symbol subframe-1 BCH
segment of the *next* frame kept at the tail of the sync window, and
`sync_unit_length` is the 1800-symbol frame drained once a subframe is decoded.

# Fields

$(TYPEDFIELDS)

# Reference

IS-GPS-800J, Sections 3.2 and 3.5, Table 3.5-1.
"""
Base.@kwdef struct GPSL1C_DConstants <: AbstractGNSSConstants
    """
    Frame length drained after each decoded subframe (1800 symbols)
    """
    sync_unit_length::Int = L1C_D_FRAME_LENGTH
    """
    Trailing next-frame subframe-1 BCH segment retained for sync (52 symbols)
    """
    preamble_length::Int = L1C_D_SUBFRAME1_LENGTH
end

"""
    GPSL1C_DReducedAlmanac

One satellite's reduced-almanac packet from subframe 3, page 3 (IS-GPS-800J
Figure 3.5-9, Table 3.5-6).

A very coarse ephemeris for satellite selection. Each page 3 carries six 33-bit
packets; this struct holds one decoded packet plus the page-level reference
week/time. It is complete in a single page (no multi-page chaining as in Galileo's
word types 7-10), so `GPSCNAV2Data.reduced_almanacs` entries are inserted whole,
keyed by `PRN_a`.

Reference values to apply (Table 3.5-6 footnotes): `e = 0`,
`δ_i = +0.0056 semi-circles` (so `i₀ = 0.30 sc = 54°` and `i₀ + δ_i = 55°`),
`Ω̇ = -2.6e-9 semi-circles/s`, `A = A_ref + δ_A` with `A_ref = 26 559 710 m`,
`Φ₀ = M₀ + ω`. Semi-circle fields are converted to radians on decode.

# Fields

  - `PRN_a::Int`: Almanac satellite PRN (8-bit field; `0` marks an empty packet,
    IS-GPS-800J §3.5.4.3.5.1.1). CNAV carries it in 6 bits, see
    [`GPSCNAVReducedAlmanac`](@ref).
  - `WN_a::Int`: Almanac reference week number (mod 8192).
  - `t_oa::Int`: Almanac reference time of week (seconds).
  - `δ_A::Float64`: Semi-major-axis delta from `A_ref` (meters).
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
    t_oa::Int
    δ_A::Float64
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

A medium-precision single-SV almanac. Each page 4 carries exactly one, complete in
that page, so `GPSCNAV2Data.midi_almanacs` entries are inserted whole, keyed by
`PRN_a`. Inclination is `δ_i` relative to `i₀ = 0.30 semi-circles` (54°);
semi-circle fields are converted to radians on decode.

# Fields

  - `PRN_a::Int`: Almanac satellite PRN.
  - `WN_a::Int`: Almanac reference week number (mod 8192).
  - `t_oa::Int`: Almanac reference time of week (seconds).
  - `e::Float64`: Eccentricity (dimensionless).
  - `δ_i::Float64`: Inclination delta from `i₀ = 0.30 sc` (rad); add the reference.
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
    t_oa::Int
    e::Float64
    δ_i::Float64
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

A page 5 carries the predict/reference times plus one DC packet (a 34-bit CDC and a
92-bit EDC segment, an indivisible pair) for *another* SV, keyed by `PRN_ID`. An
all-ones PRN ID (`0xFF`) in any PRN ID field marks an empty packet; the rest of the
block is filler (IS-GPS-800J §3.5.4.4.4.1). `dc_data_type` selects the data the
corrections apply to: `false` ⇒ CNAV-2 (`D_L1C`), `true` ⇒ legacy NAV (`D`).
Semi-circle fields → radians.

# Fields

  - `PRN_ID::Int`: PRN the corrections apply to.
  - `t_op_D::Int`: DC data predict time of week (seconds).
  - `t_OD::Int`: Time of DC data (seconds).
  - `dc_data_type::Bool`: false ⇒ CNAV-2, true ⇒ legacy NAV.
  - `δa_f0::Float64`, `δa_f1::Float64`: Clock bias / drift corrections (s, s/s).
  - `UDRA::Int`, `UDRA_dot::Int`: (rate-of-)UDRA indices (signed).
  - `Δα::Float64`, `Δβ::Float64`: Ephemeris α/β corrections (dimensionless).
  - `Δγ::Float64`: Ephemeris γ correction (rad).
  - `Δi::Float64`, `ΔΩ::Float64`: Inclination / right-ascension corrections (rad).
  - `ΔA::Float64`: Semi-major-axis correction (meters).

# Reference

IS-GPS-800J, Figure 3.5-6 / Figure 3.5-10 / Table 3.5-8.
"""
Base.@kwdef struct GPSL1C_DDifferentialCorrection
    PRN_ID::Int
    t_op_D::Int
    t_OD::Int
    dc_data_type::Bool
    δa_f0::Float64
    δa_f1::Float64
    UDRA::Int
    UDRA_dot::Int
    Δα::Float64
    Δβ::Float64
    Δγ::Float64
    Δi::Float64
    ΔΩ::Float64
    ΔA::Float64
end

"""
    GPSCNAV2Data

Decoded GPS L1C-D (CNAV-2) navigation message data.

Subframe-2 clock, ephemeris and accuracy parameters (IS-GPS-800J Figure 3.5-1 /
Table 3.5-1), subframe-3 pages 1-6 (Figures 3.5-2 to 3.5-7) and 8 (IRN-IS-800J-004
Figure 3.5-8a) and the count of CRC-valid subframe-3 pages. Field names are those of
IS-GPS-800J (see CONTEXT.md, Field naming). The ICD's satellite suffix `-n` (`M0-n`,
`af0-n`, `WNa-n`, …) is dropped throughout. Semi-circle quantities are converted to
radians (times π); all fields are `Union{Nothing,…}` until first decoded.

# Sync / timing

  - `TOI::Int`: Last validated Time-Of-Interval count (0..399), or `nothing`.
  - `ITOW::Int64`: Interval time of week — number of two-hour epochs since the
    start of the week (subframe 2 bits 14-21).
  - `WN::Int64`: Transmission week number, modulo-8192 (subframe 2 bits 1-13).
  - `t_op::Int64`: CEI data sequence propagation time of week (seconds). Distinct
    from `t_op_D` of `GPSL1C_DDifferentialCorrection`, the DC data predict time of week.

# Health / accuracy

  - `l1c_health::Bool`: L1C signal health bit (false = OK, true = bad/unavailable).
  - `URA_ED_index::Int64`: Ephemeris URA index (signed).
  - `integrity_status_flag::Bool`: Integrity level of the URA (bit 566; 3.5.3.10.1).
    `false` is legacy, `true` enhanced: "URA is integrity assured to the enhanced
    level only when the integrity status flag is '1'" (3.5.3.10).
  - `WN_op::Int64`: CEI data sequence propagation week number (bits 567-574,
    modulo 256). Pairs with `t_op` in the integrity-assured URA_NED expression of
    3.5.3.10.
  - `URA_NED0_index::Int64`, `URA_NED1_index::Int64`, `URA_NED2_index::Int64`:
    Non-elevation-dependent accuracy indices (NED Accuracy, Change, Change Rate).
    They are terms of one time-dependent NED bound, not interchangeable URAs.

# Ephemeris (Table 3.5-1)

  - `t_oe::Int64`: Ephemeris/clock data reference time of week (seconds).
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
  - `IDOT::Float64`: Rate of inclination angle (rad/s).
  - `C_is::Float64`, `C_ic::Float64`: Sine/cosine inclination harmonic corrections (rad).
  - `C_rs::Float64`, `C_rc::Float64`: Sine/cosine orbit-radius harmonic corrections (m).
  - `C_us::Float64`, `C_uc::Float64`: Sine/cosine argument-of-latitude harmonic corrections (rad).

# Clock (Table 3.5-1)

  - `t_oc::Int64`: Clock data reference time of week (seconds). CNAV-2 has no `t_oc`
    symbol (IS-GPS-800J §3.5.3.7.1: users "shall use toe, provided in bits 39 through 49
    of subframe 2, to replace toc"), so this is a copy of `t_oe`, named as on the other
    GPS messages.
  - `a_f0::Float64`, `a_f1::Float64`, `a_f2::Float64`: Clock bias / drift / drift-rate.
  - `T_GD::Float64`: L1/L2 P(Y) inter-signal correction (seconds).
  - `ISC_L1CP::Float64`, `ISC_L1CD::Float64`: L1CP / L1CD inter-signal corrections (seconds).

# Subframe 3 (IS-GPS-800J §3.5.4 — IRN-IS-800J layout)

Pages are parsed after their CRC passes, dispatching on the 6-bit page number (bits
9-14; bits 1-8 are the transmitting PRN). Layouts follow IS-GPS-800J as amended by
IRN-IS-800J-003, page 8 as amended by IRN-IS-800J-004. `num_sf3_pages_received` counts
every CRC-valid page, parsed or not.

## Page 1 — UTC + Klobuchar iono + ISC

  - `A_0,A_1,A_2::Float64`: UTC polynomial (s, s/s, s/s²; the ICD's `A_0-n`, `A_1-n`
    and `A_2-n` of Table 3.5-3, here without the satellite index `-n`).
  - `Δt_LS,Δt_LSF::Int64`: current/past and future leap-second counts (s).
  - `t_ot::Int64`: UTC reference time of week (s).
  - `WN_ot,WN_LSF::Int64`: UTC and leap-second reference week numbers.
  - `DN::Int64`: leap-second reference day number (1-7).
  - `α_0,α_1,α_2,α_3,β_0,β_1,β_2,β_3::Float64`: Klobuchar ionospheric coefficients.
  - `ISC_L1CA,ISC_L2C,ISC_L5I5,ISC_L5Q5::Float64`: inter-signal corrections (s).

## Page 2 — GGTO + EOP

  - `A_0GGTO,A_1GGTO,A_2GGTO::Float64`: GPS/GNSS time-offset polynomial.
  - `t_GGTO::Int64`, `WN_GGTO::Int64`: GGTO reference time/week.
  - `GGTO_ID::Int64`: GGTO Type ID: 0 none, 1 Galileo, 2 GLONASS,
    3-7 reserved. Called "GGTO ID" in the ICD since IRN-IS-800J-003 (before, "GNSS
    ID").
  - `t_EOP::Int64`: EOP reference time of week (s).
  - `PM_X,PM_X_dot,PM_Y,PM_Y_dot::Float64`: polar-motion values/rates.
  - `ΔUTGPS,ΔUTGPS_dot::Float64`: UT1-GPS (UT1−GPST) difference and rate.

## Pages 3/4/5 — keyed dictionaries (`nothing` until first decoded)

Keyed by the 8-bit `PRN_a` / `PRN_ID`, one preallocated slot per value
([`SlotDictionary`](@ref), iterated in ascending PRN order). A decoded packet
**overwrites** its PRN's slot in place (see [`decode!`](@ref)).

  - `reduced_almanacs::SlotDictionary{GPSL1C_DReducedAlmanac,256}` (page 3).
  - `midi_almanacs::SlotDictionary{GPSL1C_DMidiAlmanac,256}` (page 4).
  - `differential_corrections::SlotDictionary{GPSL1C_DDifferentialCorrection,256}` (page 5).

## Page 6 — Text

  - `text_message::CStaticString{29}`: up to 29 ASCII characters (control chars
    stripped), stored inline; compares equal to the matching `String`.

## Page 8 — Integrity Support Message (IRN-IS-800J-004)

IRN-IS-800J-004 replaces the page's own ISM fields (IS-GPS-800J Table 3.5-9) with the
236-bit ISM Packet of IS-GPS-200 §30.3.3.10 at bits 15-250, the same packet CNAV message
type 40 carries (see `parse_sf3_page8`).

  - `ism::GPSIntegritySupportMessage`: the most recent GPS ISM Packet whose ISM
    CRC matched; a test packet, a reserved GNSS ID or a CRC mismatch keeps the
    previous one.

## Counters

  - `num_sf3_pages_received::Int`: Count of CRC-valid subframe-3 pages received.

# Reference

IS-GPS-800J, Figures 3.5-1 through 3.5-10 and Tables 3.5-1, 3.5-3 … 3.5-8; subframe 3
page 8 as amended by IRN-IS-800J-004, its ISM Packet in IS-GPS-200 §30.3.3.10 as
amended by IRN-IS-200N-004.
"""
Base.@kwdef struct GPSCNAV2Data <: AbstractGPSCNAVData
    TOI::Union{Nothing,Int} = nothing
    ITOW::Union{Nothing,Int64} = nothing
    WN::Union{Nothing,Int64} = nothing
    t_op::Union{Nothing,Int64} = nothing

    l1c_health::Union{Nothing,Bool} = nothing
    URA_ED_index::Union{Nothing,Int64} = nothing
    integrity_status_flag::Union{Nothing,Bool} = nothing
    WN_op::Union{Nothing,Int64} = nothing
    URA_NED0_index::Union{Nothing,Int64} = nothing
    URA_NED1_index::Union{Nothing,Int64} = nothing
    URA_NED2_index::Union{Nothing,Int64} = nothing

    t_oe::Union{Nothing,Int64} = nothing
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
    IDOT::Union{Nothing,Float64} = nothing
    C_is::Union{Nothing,Float64} = nothing
    C_ic::Union{Nothing,Float64} = nothing
    C_rs::Union{Nothing,Float64} = nothing
    C_rc::Union{Nothing,Float64} = nothing
    C_us::Union{Nothing,Float64} = nothing
    C_uc::Union{Nothing,Float64} = nothing

    t_oc::Union{Nothing,Int64} = nothing
    a_f0::Union{Nothing,Float64} = nothing
    a_f1::Union{Nothing,Float64} = nothing
    a_f2::Union{Nothing,Float64} = nothing
    T_GD::Union{Nothing,Float64} = nothing
    ISC_L1CP::Union{Nothing,Float64} = nothing
    ISC_L1CD::Union{Nothing,Float64} = nothing

    # --- Subframe 3, page 1: UTC + iono + ISC (IS-GPS-800J Fig 3.5-2) ---
    A_0::Union{Nothing,Float64} = nothing
    A_1::Union{Nothing,Float64} = nothing
    A_2::Union{Nothing,Float64} = nothing
    Δt_LS::Union{Nothing,Int64} = nothing
    t_ot::Union{Nothing,Int64} = nothing
    WN_ot::Union{Nothing,Int64} = nothing
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
    ΔUTGPS::Union{Nothing,Float64} = nothing
    ΔUTGPS_dot::Union{Nothing,Float64} = nothing

    # --- Subframe 3, pages 3/4/5: per-SV keyed dictionaries ---
    reduced_almanacs::Union{
        Nothing,
        SlotDictionary{GPSL1C_DReducedAlmanac,L1C_D_PRN_SLOTS},
    } = nothing
    midi_almanacs::Union{Nothing,SlotDictionary{GPSL1C_DMidiAlmanac,L1C_D_PRN_SLOTS}} =
        nothing
    differential_corrections::Union{
        Nothing,
        SlotDictionary{GPSL1C_DDifferentialCorrection,L1C_D_PRN_SLOTS},
    } = nothing

    # --- Subframe 3, page 6: text (IS-GPS-800J Fig 3.5-7) ---
    text_message::Union{Nothing,CStaticString{L1C_D_TEXT_LENGTH}} = nothing

    # --- Subframe 3, page 8: ISM Packet (IRN-IS-800J-004 Fig 3.5-8a) ---
    ism::Union{Nothing,GPSIntegritySupportMessage} = nothing

    num_sf3_pages_received::Int = 0
end

# `WN_ot` and `WN_LSF` are 13 bits (IS-GPS-800J Figure 3.5-2, subframe 3 page 1), like
# `WN`: see `get_utc_parameters`.
const CNAV2_UTC_WN_MODULUS = 8192
const CNAV2_LSF_WN_MODULUS = 8192

@inline function GPSCNAV2Data(
    data::GPSCNAV2Data;
    TOI = data.TOI,
    ITOW = data.ITOW,
    WN = data.WN,
    t_op = data.t_op,
    l1c_health = data.l1c_health,
    URA_ED_index = data.URA_ED_index,
    integrity_status_flag = data.integrity_status_flag,
    WN_op = data.WN_op,
    URA_NED0_index = data.URA_NED0_index,
    URA_NED1_index = data.URA_NED1_index,
    URA_NED2_index = data.URA_NED2_index,
    t_oe = data.t_oe,
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
    IDOT = data.IDOT,
    C_is = data.C_is,
    C_ic = data.C_ic,
    C_rs = data.C_rs,
    C_rc = data.C_rc,
    C_us = data.C_us,
    C_uc = data.C_uc,
    t_oc = data.t_oc,
    a_f0 = data.a_f0,
    a_f1 = data.a_f1,
    a_f2 = data.a_f2,
    T_GD = data.T_GD,
    ISC_L1CP = data.ISC_L1CP,
    ISC_L1CD = data.ISC_L1CD,
    A_0 = data.A_0,
    A_1 = data.A_1,
    A_2 = data.A_2,
    Δt_LS = data.Δt_LS,
    t_ot = data.t_ot,
    WN_ot = data.WN_ot,
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
    ΔUTGPS = data.ΔUTGPS,
    ΔUTGPS_dot = data.ΔUTGPS_dot,
    reduced_almanacs = data.reduced_almanacs,
    midi_almanacs = data.midi_almanacs,
    differential_corrections = data.differential_corrections,
    text_message = data.text_message,
    ism = data.ism,
    num_sf3_pages_received = data.num_sf3_pages_received,
)
    GPSCNAV2Data(
        TOI,
        ITOW,
        WN,
        t_op,
        l1c_health,
        URA_ED_index,
        integrity_status_flag,
        WN_op,
        URA_NED0_index,
        URA_NED1_index,
        URA_NED2_index,
        t_oe,
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
        IDOT,
        C_is,
        C_ic,
        C_rs,
        C_rc,
        C_us,
        C_uc,
        t_oc,
        a_f0,
        a_f1,
        a_f2,
        T_GD,
        ISC_L1CP,
        ISC_L1CD,
        A_0,
        A_1,
        A_2,
        Δt_LS,
        t_ot,
        WN_ot,
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
        ΔUTGPS,
        ΔUTGPS_dot,
        reduced_almanacs,
        midi_almanacs,
        differential_corrections,
        text_message,
        ism,
        num_sf3_pages_received,
    )
end

# The default `==` is `===`, which fails for the mutable `SlotDictionary` fields even
# when contents match. Compare field-by-field (mirrors `GalileoINAVData`).
Base.:(==)(a::GPSCNAV2Data, b::GPSCNAV2Data) = fields_equal(a, b)

"""
$(TYPEDEF)

Per-decoder cache for the GPS L1C-D signal.

Holds the soft-symbol `CircularDeque{Float32}` (capacity 1852 = 1800 frame + 52 next-frame
BCH segment) and two Aff3ct LDPC belief-propagation decoders (subframe 2: K=600, N=1200;
subframe 3: K=274, N=548), loaded at construction from the committed `.alist` matrices
(`data/`, else the embedded copies). The decoders are mutable Aff3ct handles shared by
reference through the otherwise-immutable [`GNSSDecoderState`](@ref).

Each decoder comes in an `LDPCScratch` holding the buffers its decode would otherwise
allocate per frame; the two payload buffers do the same for the deinterleave, so a
decoded frame allocates nothing.

# Fields

$(TYPEDFIELDS)
"""
struct GPSL1C_DCache <: AbstractGNSSCache
    """
    Soft-symbol buffer (1852 = 1800 frame + 52 next subframe-1)
    """
    soft_buffer::CircularDeque{Float32}
    """
    Subframe-2 LDPC decoder and its scratch buffers (K=600, N=1200)
    """
    sf2_ldpc::LDPCScratch
    """
    Subframe-3 LDPC decoder and its scratch buffers (K=274, N=548)
    """
    sf3_ldpc::LDPCScratch
    """
    Polarity-resolved 1748-symbol interleaved payload, copied out per frame
    """
    payload::Vector{Float32}
    """
    The same payload deinterleaved: subframe 2 then subframe 3
    """
    deinterleaved::Vector{Float32}
    """
    Preallocated containers `raw_data` and `data` are decoded into (the
    PRN-keyed subframe-3 stores), overwritten in place
    """
    storage::DataStorage{GPSCNAV2Data}
end

# Nothing satellite-specific (see `reset_decoder_state!(state, prn)`): the buffers are
# FEC scratch, refilled every frame.
clear_satellite_cache(cache::GPSL1C_DCache) = cache

function GPSL1C_DCache()
    # The CNAV-2 LDPC codes are systematic with codeword layout [info | parity]
    # (IS-GPS-800J §3.2.3.4: H = [[A B T] [C D E]] acts on [info; parity]);
    # `load_ldpc_decoder` forces that layout onto Aff3ct. Verified against a Spirent
    # GSS post-FEC L1C recording.
    GPSL1C_DCache(
        CircularDeque{Float32}(L1C_D_WINDOW_LENGTH),
        committed_ldpc_scratch("cnv2_sf2.alist"),
        committed_ldpc_scratch("cnv2_sf3.alist"),
        Vector{Float32}(undef, L1C_D_PAYLOAD_SYMBOLS),
        Vector{Float32}(undef, L1C_D_PAYLOAD_SYMBOLS),
        DataStorage{GPSCNAV2Data}(),
    )
end

# Container fields at ICD size: one slot per 8-bit `PRN_a` / `PRN_ID` value.
preallocated_data(::Type{GPSCNAV2Data}) = GPSCNAV2Data(;
    reduced_almanacs = SlotDictionary{GPSL1C_DReducedAlmanac,L1C_D_PRN_SLOTS}(),
    midi_almanacs = SlotDictionary{GPSL1C_DMidiAlmanac,L1C_D_PRN_SLOTS}(),
    differential_corrections = SlotDictionary{
        GPSL1C_DDifferentialCorrection,
        L1C_D_PRN_SLOTS,
    }(),
)

# The LDPC handles are runtime Aff3ct objects: caches are equal when their soft
# buffers match.
function Base.:(==)(a::GPSL1C_DCache, b::GPSL1C_DCache)
    deques_equal(a.soft_buffer, b.soft_buffer)
end

function is_subframe2_decoded(data::GPSCNAV2Data)
    !isnothing(data.WN) &&
        !isnothing(data.ITOW) &&
        !isnothing(data.t_oe) &&
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

"""
$(TYPEDSIGNATURES)

Seconds of week of the epoch the TOI stamps: `ITOW·7200 + TOI·18`, with `TOI = 0`
counted into the *next* two-hour interval.

The time of week is split across two parts: subframe-2 `ITOW` counts two-hour
intervals since the start of the week (IS-GPS-800J §3.5.3.2) and subframe-1 `TOI` counts
18-second frames within the interval (§3.5.2). Neither is a time of week alone; `nothing`
until both are present.

The `TOI = 0` special case is the ICD's. §3.5.2 defines TOI relative to "the two-hour
period represented by ITOW count in the subframe 2 of the *next* 18-second frame"
(range 0–399): the frame broadcasting `TOI = 0` is the last of the old interval, so
its own `ITOW` is one interval behind. Without the carry, one frame in 400 would stamp
a time two hours in the past. At the week crossover (`ITOW = 83`, `TOI = 0`) the result
wraps to second 0 of the *following* week while the paired `WN` is still the old week; a
consumer combining it with `WN` inherits that one-frame-per-week
ambiguity, also the ICD's.
"""
function get_broadcast_time_of_week(data::GPSCNAV2Data)
    (isnothing(data.ITOW) || isnothing(data.TOI)) && return nothing
    itow = Int(data.ITOW) + (data.TOI == 0)
    return mod(
        itow * L1C_D_ITOW_SECONDS + Int(data.TOI) * L1C_D_FRAME_SECONDS,
        SECONDS_PER_WEEK,
    )
end

# CNAV-2's week number is 13 bits, absolute until 2137.
get_broadcast_week(data::GPSCNAV2Data; approximate_year::Integer) = data.WN

# CNAV-2 carries the full L1 ISC set: L1C-D, L1C-P or C/A, each against T_GD.
get_group_delay(state::GNSSDecoderState{<:GPSCNAV2Data}, ::GPSL1C_D) =
    cnav_group_delay(state, state.data.ISC_L1CD)
get_group_delay(state::GNSSDecoderState{<:GPSCNAV2Data}, ::GPSL1C_P) =
    cnav_group_delay(state, state.data.ISC_L1CP)
get_group_delay(state::GNSSDecoderState{<:GPSCNAV2Data}, ::GPSL1CA) =
    cnav_group_delay(state, state.data.ISC_L1CA)

"""
$(TYPEDSIGNATURES)

The GPS/GNSS time offset from subframe 3 page 2, if this satellite is currently
broadcasting one for `target`.

`GGTO_ID` names the system the single set refers to (0 none, 1 Galileo, 2 GLONASS,
3-7 reserved; IS-GPS-800J §3.5.4.2.1), as in CNAV message type 35; see
`gps_ggto_offset`.
"""
time_offset_parameters(state::GNSSDecoderState{<:GPSCNAV2Data}, target::TimeSystem) =
    gps_ggto_offset(state, target)

function has_positioning_messages(data::GPSCNAV2Data)
    !isnothing(data.TOI) && is_subframe2_decoded(data)
end

"""
$(TYPEDSIGNATURES)

Create a decoder state for GPS L1C-D (CNAV-2) navigation messages.

Wires up a [`GNSSDecoderState`](@ref) with a 1852-symbol soft-symbol buffer, the
400-entry BCH(51,8) TOI codeword table (`src/coding/bch_toi.jl`) and two Aff3ct LDPC
decoders loaded at construction from the committed `.alist` matrices (`data/`, else the
embedded copies).

# Arguments

  - `prn::Int`: Pseudo-Random Noise code identifier (1-63 for L1C).

# Returns

  - `GNSSDecoderState{GPSCNAV2Data}`: Initialized decoder state for GPS L1C-D.

# Example

```julia
state = GPSL1C_DDecoderState(1)            # PRN 1
state = decode!(state, soft_bits, num_bits)
if is_sat_healthy(state)
    # Use state.data for positioning
end
```

# See Also

  - [`GNSSDecoderState`](@ref): The underlying state structure
  - [`decode!`](@ref): Decode soft symbols using this state
  - [`reset_decoder_state!`](@ref): Reset after signal loss
"""
function GPSL1C_DDecoderState(prn)
    GNSSDecoderState(
        prn,
        GPSCNAV2Data(),
        GPSCNAV2Data(),
        GPSL1C_DConstants(),
        GPSL1C_DCache(),
        nothing,
        false,
    )
end

function GNSSDecoderState(system::GPSL1C_D, prn)
    GPSL1C_DDecoderState(prn)
end

# Signal metadata is forwarded through this (see `src/gps/l1ca.jl`). CNAV-2 is on the
# data component `GPSL1C_D`; `GPSL1C_P` is the dataless pilot.
get_signal_type(::GPSL1C_DConstants) = GPSL1C_D

"""
$(TYPEDSIGNATURES)

Reset the GPS L1C-D decoder state after a signal loss or reacquisition.

Clears the sync state (soft buffer emptied in place, `TOI` and `ITOW` dropped) and the
validated data, keeping the long-lived ephemeris/clock fields in `raw_data` so a
`GNSSReceiver` can reuse the satellite after reacquisition without re-decoding
subframe 2. Same semantics as the GPS L1 C/A and Galileo E1B decoders.

# Arguments

  - `state::GNSSDecoderState{<:GPSCNAV2Data}`: Current GPS L1C-D decoder state.

# Returns

  - `GNSSDecoderState{<:GPSCNAV2Data}`: Reset decoder state with cleared buffers.

# See Also

  - [`GPSL1C_DDecoderState`](@ref): Create a fresh decoder state
  - [`decode!`](@ref): Continue decoding after reset
"""
function reset_decoder_state!(state::GNSSDecoderState{<:GPSCNAV2Data})
    empty!(state.cache.soft_buffer)
    GNSSDecoderState(
        state;
        # `ITOW` is cleared with `TOI`: they only mean a time as a pair, and an outage
        # can span any number of interval boundaries, so a kept `ITOW` would date the
        # first re-locked frames whole intervals in the past. Re-learning costs one
        # decoded subframe 2.
        raw_data = GPSCNAV2Data(state.raw_data; TOI = nothing, ITOW = nothing),
        data = GPSCNAV2Data(),
        num_bits_since_time_reference = nothing,
        is_shifted_by_180_degrees = false,
    )
end

# ---- Sync ------------------------------------------------------------------
#
# Overrides the generic packed-buffer sync: the BCH(51,8) TOI match runs directly on
# the soft-symbol deque. The first 52 symbols must equal the codeword for some `TOI`
# and the last 52 that for `(TOI+1) mod 400`, in either polarity (`sync_bch_toi`,
# src/coding/bch_toi.jl).

"""
    try_sync(state::GNSSDecoderState{<:GPSCNAV2Data}) -> Union{Nothing,BCHToiSync}

Sync hook for GPS L1C-D. Hard-slices the leading and trailing 52 soft symbols of the
1852-symbol window and runs `sync_bch_toi`. Returns the [`BCHToiSync`](@ref) (TOI and
polarity flip) on a match, else `nothing`.
"""
function try_sync(state::GNSSDecoderState{<:GPSCNAV2Data})
    deque = soft_buffer(state)
    # `pack_bits_lsb_first` (src/gnss.jl) hard-slices straight off the deque in
    # `pack_hard_codeword`'s bit order; this runs per symbol, so materialising two
    # 52-element slices would allocate every time.
    first52 = pack_bits_lsb_first(deque, 1, L1C_D_SUBFRAME1_LENGTH)
    next52 = pack_bits_lsb_first(deque, L1C_D_FRAME_LENGTH + 1, L1C_D_SUBFRAME1_LENGTH)
    sync_bch_toi(first52, next52)
end

"""
    complement_buffer_if_necessary(state::GNSSDecoderState{<:GPSCNAV2Data}, sync)

Resolve the 180° polarity ambiguity for L1C-D: record the polarity already detected
by `sync_bch_toi` in `sync::BCHToiSync` on the state and pass `sync` through
unchanged to `decode_sync_unit`.
"""
function complement_buffer_if_necessary(
    state::GNSSDecoderState{<:GPSCNAV2Data},
    sync::BCHToiSync,
)
    GNSSDecoderState(state; is_shifted_by_180_degrees = sync.polarity_flipped), sync
end

# ---- Subframe pipeline -----------------------------------------------------

"""
    decode_sync_unit(state::GNSSDecoderState{<:GPSCNAV2Data}, sync::BCHToiSync)

Process one locked CNAV-2 frame. Enforces a monotonic TOI (`(previous + 1) mod 400`;
a break resets the decoder to searching) and resolves the TOI/polarity complement
ambiguity by continuity, then deinterleaves and LDPC-decodes
subframes 2 and 3 from the 1748-symbol payload, CRC-checking each, and parses them
into `raw_data`.
"""
function decode_sync_unit(state::GNSSDecoderState{<:GPSCNAV2Data}, sync::BCHToiSync)
    # `sync` carries the *low* branch of the BCH complement ambiguity (`sync_bch_toi`
    # cannot tell (TOI, polarity) from (TOI + 256, !polarity) when both are in range
    # and reports the lower TOI). Resolve it by TOI continuity when a previous frame is
    # held (a locked frame's TOI follows the previous by +1 mod 400 on whichever branch
    # it settled); any other step is loss of frame lock.
    TOI = sync.TOI
    flipped = sync.polarity_flipped
    prev_toi = state.raw_data.TOI
    if !isnothing(prev_toi)
        expected = (prev_toi + 1) % TOI_RANGE
        if expected == TOI + TOI_COMPLEMENT_OFFSET
            TOI += TOI_COMPLEMENT_OFFSET
            flipped = !flipped
        elseif expected != TOI
            return reset_decoder_state!(state)
        end
    end
    # `TOI` advances every locked frame, but `ITOW` only refreshes when a subframe 2
    # clears LDPC+CRC. The frame at `TOI = 1` is the first to broadcast the new `ITOW`
    # (`TOI = 0` still carries the old one; see `get_broadcast_time_of_week`), so carry
    # the interval across here, or a noisy subframe 2 from then on promotes a time of
    # week exactly 7200 s low. A later decoded subframe 2 overwrites it either way.
    # Continuity above makes `prev_toi == 0` here. The week crossover (ITOW 83 → 0 while
    # `WN` awaits its subframe 2) has the same one-frame ambiguity as every WN signal.
    ITOW = state.raw_data.ITOW
    if !isnothing(prev_toi) && TOI == 1 && !isnothing(ITOW)
        ITOW = (ITOW + 1) % L1C_D_ITOW_RANGE
    end
    # split: a Union keyword value takes the allocating kw path on Julia 1.10
    raw_data = @split_nothing ITOW GPSCNAV2Data(state.raw_data; TOI, ITOW)
    state = GNSSDecoderState(state; is_shifted_by_180_degrees = flipped, raw_data)

    # The 1748-symbol interleaved SF2+SF3 payload (symbols 53..1800), polarity-flipped
    # by negating soft symbols up front.
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

    # The explicit success flag matters: `GNSSDecoderState` is immutable, so
    # `decoded === state` would also hold when subframe 2 decoded but repeated what
    # `raw_data` already held (e.g. after a reset), wrongly triggering the complement retry.
    decoded, sf2_ok = decode_subframe2(state, sf2_symbols)
    complement =
        TOI < TOI_COMPLEMENT_OFFSET ? TOI + TOI_COMPLEMENT_OFFSET :
        TOI - TOI_COMPLEMENT_OFFSET
    if !sf2_ok && complement < TOI_RANGE
        # Subframe 2 failed on this branch. Its CRC is the only in-band oracle, so try
        # the complement branch (negating the payload is exactly the 180° flip, and LLR
        # negation commutes with the deinterleaver). This recovers a cold start with
        # true TOI ≥ 256 and un-sticks a chain that settled on the wrong branch while
        # its subframe 2 was noisy.
        deinterleaved .= .-deinterleaved
        alt_state = GNSSDecoderState(
            state;
            is_shifted_by_180_degrees = !flipped,
            raw_data = GPSCNAV2Data(state.raw_data; TOI = complement),
        )
        alt_decoded, alt_ok = decode_subframe2(alt_state, sf2_symbols)
        if alt_ok
            return decode_subframe3(alt_decoded, sf3_symbols)
        end
        # Both failed: noise. Restore the polarity and keep the original hypothesis.
        deinterleaved .= .-deinterleaved
    end
    state = decoded
    state = decode_subframe3(state, sf3_symbols)
    return state
end

# ---- Subframe 2 bit-field extraction (IS-GPS-800J Figure 3.5-1) ------------
#
# `word` is the 600-bit subframe-2 info block packed MSB-first into a `UInt600`
# (bit 1 = MSB). Fields are read by 1-based start bit and length through the
# shared `get_bits` / `get_twos_complement_num` / `get_bit` helpers.

# Returns `(state, decoded::Bool)`; `decoded` is `false` on a CRC failure.
function decode_subframe2(state::GNSSDecoderState{<:GPSCNAV2Data}, sf2_symbols)
    word = ldpc_decode_word(state.cache.sf2_ldpc, sf2_symbols, UInt600)
    isnothing(word) && return state, false  # silently drop on CRC failure
    word_length = L1C_D_SF2_INFO_BITS
    WN = Int(get_bits(word, word_length, 1, 13))
    ITOW = Int(get_bits(word, word_length, 14, 8))
    t_op = Int(get_bits(word, word_length, 22, 11)) * 300
    l1c_health = get_bit(word, word_length, 33)
    URA_ED_index = get_twos_complement_num(word, word_length, 34, 5)
    t_oe = Int(get_bits(word, word_length, 39, 11)) * 300
    ΔA = get_twos_complement_num(word, word_length, 50, 26) * 2.0^-9
    A_dot = get_twos_complement_num(word, word_length, 76, 25) * 2.0^-21
    Δn_0 = get_twos_complement_num(word, word_length, 101, 17) * 2.0^-44 * GNSS_PI
    Δn_0_dot = get_twos_complement_num(word, word_length, 118, 23) * 2.0^-57 * GNSS_PI
    M_0 = get_twos_complement_num(word, word_length, 141, 33) * 2.0^-32 * GNSS_PI
    # 33-bit field: convert through Int64 explicitly (a platform Int would
    # overflow on 32-bit), matching the BeiDou B-CNAV parsers.
    e = Int64(get_bits(word, word_length, 174, 33)) * 2.0^-34
    ω = get_twos_complement_num(word, word_length, 207, 33) * 2.0^-32 * GNSS_PI
    Ω_0 = get_twos_complement_num(word, word_length, 240, 33) * 2.0^-32 * GNSS_PI
    i_0 = get_twos_complement_num(word, word_length, 273, 33) * 2.0^-32 * GNSS_PI
    ΔΩ_dot = get_twos_complement_num(word, word_length, 306, 17) * 2.0^-44 * GNSS_PI
    IDOT = get_twos_complement_num(word, word_length, 323, 15) * 2.0^-44 * GNSS_PI
    C_is = get_twos_complement_num(word, word_length, 338, 16) * 2.0^-30
    C_ic = get_twos_complement_num(word, word_length, 354, 16) * 2.0^-30
    C_rs = get_twos_complement_num(word, word_length, 370, 24) * 2.0^-8
    C_rc = get_twos_complement_num(word, word_length, 394, 24) * 2.0^-8
    C_us = get_twos_complement_num(word, word_length, 418, 21) * 2.0^-30
    C_uc = get_twos_complement_num(word, word_length, 439, 21) * 2.0^-30
    URA_NED0_index = get_twos_complement_num(word, word_length, 460, 5)
    URA_NED1_index = Int(get_bits(word, word_length, 465, 3))
    URA_NED2_index = Int(get_bits(word, word_length, 468, 3))
    a_f0 = get_twos_complement_num(word, word_length, 471, 26) * 2.0^-35
    a_f1 = get_twos_complement_num(word, word_length, 497, 20) * 2.0^-48
    a_f2 = get_twos_complement_num(word, word_length, 517, 10) * 2.0^-60
    T_GD = get_twos_complement_num(word, word_length, 527, 13) * 2.0^-35
    ISC_L1CP = get_twos_complement_num(word, word_length, 540, 13) * 2.0^-35
    ISC_L1CD = get_twos_complement_num(word, word_length, 553, 13) * 2.0^-35
    # Bits 566 and 567-574 close subframe 2 before its 2 reserved bits and the CRC.
    # Easy to miss: Figure 3.5-1 puts `ISF - 1 BIT` and `RESERVED - 2 BITS` on a
    # second label row.
    integrity_status_flag = get_bit(word, word_length, 566)
    WN_op = Int(get_bits(word, word_length, 567, 8))

    raw = GPSCNAV2Data(
        state.raw_data;
        WN,
        ITOW,
        t_op,
        l1c_health,
        URA_ED_index,
        integrity_status_flag,
        WN_op,
        t_oe,
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
        IDOT,
        C_is,
        C_ic,
        C_rs,
        C_rc,
        C_us,
        C_uc,
        URA_NED0_index,
        URA_NED1_index,
        URA_NED2_index,
        t_oc = t_oe,  # CNAV-2 shares one reference time for clock and ephemeris
        a_f0,
        a_f1,
        a_f2,
        T_GD,
        ISC_L1CP,
        ISC_L1CD,
    )
    GNSSDecoderState(state; raw_data = raw), true
end

# ---- Subframe 3 page parsing (IS-GPS-800J §3.5.4) --------------------------
#
# Every SF3 page is a 274-bit info block: bits 1-8 the transmitting PRN, bits 9-14 the
# page number, then page fields, with a trailing CRC-24Q. After the CRC passes, the bits
# are packed MSB-first into a `UInt288` and the parsed fields merged into `raw_data`
# (as for SF2), except for the PRN-keyed pages 3/4/5, whose packets overwrite their slot
# of a preallocated store in place. Layouts follow IS-GPS-800J as amended by
# IRN-IS-800J-003, which for page 1 only restores Greek letters the base PDF rendered as
# question marks; the four ISC fields at bits 177/190/203/216 are already in base Rev J
# (Figure 3.5-2, Table 6.2-18). Page 8 follows IRN-IS-800J-004 (`parse_sf3_page8`).

function decode_subframe3(state::GNSSDecoderState{<:GPSCNAV2Data}, sf3_symbols)
    word = ldpc_decode_word(state.cache.sf3_ldpc, sf3_symbols, UInt288)
    isnothing(word) && return state  # silently drop on CRC failure

    # CRC-valid page received: count it, then dispatch on the page number.
    raw = GPSCNAV2Data(
        state.raw_data;
        num_sf3_pages_received = state.raw_data.num_sf3_pages_received + 1,
    )

    word_length = L1C_D_SF3_INFO_BITS
    page = Int(get_bits(word, word_length, 9, 6))  # bits 1-8 PRN, bits 9-14 page (IS-GPS-800J §3.5.4)
    raw = if page == L1C_D_SF3_PAGE_UTC_IONO
        parse_sf3_page1(raw, word)
    elseif page == L1C_D_SF3_PAGE_GGTO_EOP
        parse_sf3_page2(raw, word)
    elseif page == L1C_D_SF3_PAGE_REDUCED_ALMANAC
        parse_sf3_page3(raw, state.cache.storage.raw, word)
    elseif page == L1C_D_SF3_PAGE_MIDI_ALMANAC
        parse_sf3_page4(raw, state.cache.storage.raw, word)
    elseif page == L1C_D_SF3_PAGE_DIFF_CORRECTION
        parse_sf3_page5(raw, state.cache.storage.raw, word)
    elseif page == L1C_D_SF3_PAGE_TEXT
        parse_sf3_page6(raw, word)
    elseif page == L1C_D_SF3_PAGE_ISM
        parse_sf3_page8(raw, word)
    else
        raw  # unsupported/reserved page (e.g. 7 SV-config): counted, ignored
    end

    GNSSDecoderState(state; raw_data = raw)
end

"""
Subframe 3, page 1 — UTC + Klobuchar iono + ISC (IS-GPS-800J Fig 3.5-2, Table 3.5-3).
"""
function parse_sf3_page1(raw::GPSCNAV2Data, word::UInt288)
    word_length = L1C_D_SF3_INFO_BITS
    GPSCNAV2Data(
        raw;
        # UTC polynomial (Table 3.5-3).
        A_0 = get_twos_complement_num(word, word_length, 15, 16) * 2.0^-35,
        A_1 = get_twos_complement_num(word, word_length, 31, 13) * 2.0^-51,
        A_2 = get_twos_complement_num(word, word_length, 44, 7) * 2.0^-68,
        Δt_LS = get_twos_complement_num(word, word_length, 51, 8),
        t_ot = Int(get_bits(word, word_length, 59, 16)) * 2^4,
        WN_ot = Int(get_bits(word, word_length, 75, 13)),
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
function parse_sf3_page2(raw::GPSCNAV2Data, word::UInt288)
    word_length = L1C_D_SF3_INFO_BITS
    GPSCNAV2Data(
        raw;
        # GGTO (Table 3.5-4). Field order in Fig 3.5-3: tGGTO, WNGGTO, A0, A1, A2.
        t_GGTO = Int(get_bits(word, word_length, 18, 16)) * 2^4,
        WN_GGTO = Int(get_bits(word, word_length, 34, 13)),
        A_0GGTO = get_twos_complement_num(word, word_length, 47, 16) * 2.0^-35,
        A_1GGTO = get_twos_complement_num(word, word_length, 63, 13) * 2.0^-51,
        A_2GGTO = get_twos_complement_num(word, word_length, 76, 7) * 2.0^-68,
        GGTO_ID = Int(get_bits(word, word_length, 15, 3)),
        # EOP (Table 3.5-5). Fields are contiguous; Figure 3.5-3 only *draws* PM_X
        # across its 100-bit row boundary (2 MSBs end row 1, 19 LSBs begin row 2), so
        # the plain 21-bit read at bit 99 is correct (likewise ΔUTGPS).
        t_EOP = Int(get_bits(word, word_length, 83, 16)) * 2^4,
        PM_X = get_twos_complement_num(word, word_length, 99, 21) * 2.0^-20,
        PM_X_dot = get_twos_complement_num(word, word_length, 120, 15) * 2.0^-21,
        PM_Y = get_twos_complement_num(word, word_length, 135, 21) * 2.0^-20,
        PM_Y_dot = get_twos_complement_num(word, word_length, 156, 15) * 2.0^-21,
        ΔUTGPS = get_twos_complement_num(word, word_length, 171, 31) * 2.0^-23,
        ΔUTGPS_dot = get_twos_complement_num(word, word_length, 202, 19) * 2.0^-25,
    )
end

"""
Decode one 33-bit reduced-almanac packet starting at 1-based bit `start` (IS-GPS-800J Fig 3.5-9).
"""
function _reduced_almanac_packet(word::UInt288, start::Int, WN_a::Int, t_oa::Int)
    word_length = L1C_D_SF3_INFO_BITS
    PRN_a = Int(get_bits(word, word_length, start, 8))
    # PRNa == 0 marks an empty packet; per IS-GPS-800J §3.5.4.3.5.1.1 all later bits are
    # then filler, so the caller
    # stops here rather than parsing the remaining packets.
    PRN_a == 0 && return nothing
    GPSL1C_DReducedAlmanac(;
        PRN_a,
        WN_a,
        t_oa,
        δ_A = get_twos_complement_num(word, word_length, start + 8, 8) * 2.0^9,
        Ω_0 = get_twos_complement_num(word, word_length, start + 16, 7) * 2.0^-6 * GNSS_PI,
        Φ_0 = get_twos_complement_num(word, word_length, start + 23, 7) * 2.0^-6 * GNSS_PI,
        l1_health = get_bit(word, word_length, start + 30),
        l2_health = get_bit(word, word_length, start + 31),
        l5_health = get_bit(word, word_length, start + 32),
    )
end

"""
Subframe 3, page 3 — six reduced-almanac packets (IS-GPS-800J Fig 3.5-4).

Each packet **overwrites** its PRN's slot of `raw.reduced_almanacs` in place, or, on the
first page 3, of the preallocated `storage.reduced_almanacs`.
"""
function parse_sf3_page3(raw::GPSCNAV2Data, storage::GPSCNAV2Data, word::UInt288)
    almanacs = writable_container(raw.reduced_almanacs, storage.reduced_almanacs)
    word_length = L1C_D_SF3_INFO_BITS
    WN_a = Int(get_bits(word, word_length, 15, 13))
    t_oa = Int(get_bits(word, word_length, 28, 8)) * 2^12
    # Six 33-bit packets at bits 36, 69, 102, 135, 168, 201.
    for start in (36, 69, 102, 135, 168, 201)
        packet = _reduced_almanac_packet(word, start, WN_a, t_oa)
        isnothing(packet) && break  # PRNa==0 ⇒ rest of page is filler (§3.5.4.3.5.1.1)
        set!(almanacs, packet.PRN_a, packet)
        raw = GPSCNAV2Data(raw; reduced_almanacs = almanacs)
    end
    return raw
end

"""
Subframe 3, page 4 — one Midi almanac (IS-GPS-800J Fig 3.5-5, Table 3.5-7).

The almanac **overwrites** its PRN's slot of `raw.midi_almanacs` in place, or, on the
first page 4, of the preallocated `storage.midi_almanacs`.
"""
function parse_sf3_page4(raw::GPSCNAV2Data, storage::GPSCNAV2Data, word::UInt288)
    word_length = L1C_D_SF3_INFO_BITS
    PRN_a = Int(get_bits(word, word_length, 36, 8))
    PRN_a == 0 && return raw  # empty almanac
    alm = GPSL1C_DMidiAlmanac(;
        PRN_a,
        WN_a = Int(get_bits(word, word_length, 15, 13)),
        t_oa = Int(get_bits(word, word_length, 28, 8)) * 2^12,
        l1_health = get_bit(word, word_length, 44),
        l2_health = get_bit(word, word_length, 45),
        l5_health = get_bit(word, word_length, 46),
        e = Int(get_bits(word, word_length, 47, 11)) * 2.0^-16,
        δ_i = get_twos_complement_num(word, word_length, 58, 11) * 2.0^-14 * GNSS_PI,
        Ω_dot = get_twos_complement_num(word, word_length, 69, 11) * 2.0^-33 * GNSS_PI,
        sqrt_A = Int(get_bits(word, word_length, 80, 17)) * 2.0^-4,
        Ω_0 = get_twos_complement_num(word, word_length, 97, 16) * 2.0^-15 * GNSS_PI,
        ω = get_twos_complement_num(word, word_length, 113, 16) * 2.0^-15 * GNSS_PI,
        M_0 = get_twos_complement_num(word, word_length, 129, 16) * 2.0^-15 * GNSS_PI,
        a_f0 = get_twos_complement_num(word, word_length, 145, 11) * 2.0^-20,
        a_f1 = get_twos_complement_num(word, word_length, 156, 10) * 2.0^-37,
    )
    midi_almanacs = writable_container(raw.midi_almanacs, storage.midi_almanacs)
    set!(midi_almanacs, PRN_a, alm)
    GPSCNAV2Data(raw; midi_almanacs)
end

"""
Subframe 3, page 5 — one differential-correction packet (IS-GPS-800J Fig 3.5-6/3.5-10, Table 3.5-8).

The packet **overwrites** its PRN's slot of `raw.differential_corrections` in
place, or, on the first page 5, of the preallocated `storage.differential_corrections`.
"""
function parse_sf3_page5(raw::GPSCNAV2Data, storage::GPSCNAV2Data, word::UInt288)
    word_length = L1C_D_SF3_INFO_BITS
    # Page-level fields precede the 126-bit CDC+EDC packet. Layout (Fig 3.5-6):
    # bit 15 t_op-D (11, scale 300), bit 26 t_OD (11, scale 300),
    # bit 37 DC data type (1), then the CDC segment (bit 38) and EDC segment.
    t_op_D = Int(get_bits(word, word_length, 15, 11)) * 300
    t_OD = Int(get_bits(word, word_length, 26, 11)) * 300
    dc_data_type = get_bit(word, word_length, 37)
    # CDC segment (Fig 3.5-10): PRN ID(8) δaf0(13) δaf1(8) UDRA(5) — starts bit 38.
    cdc_start_bit = 38
    PRN_ID = Int(get_bits(word, word_length, cdc_start_bit, 8))
    PRN_ID == 0xff && return raw  # all-ones PRN ⇒ no DC data in this packet
    δa_f0 = get_twos_complement_num(word, word_length, cdc_start_bit + 8, 13) * 2.0^-35
    δa_f1 = get_twos_complement_num(word, word_length, cdc_start_bit + 21, 8) * 2.0^-51
    UDRA = get_twos_complement_num(word, word_length, cdc_start_bit + 29, 5)
    # EDC segment (Fig 3.5-10): PRN ID(8) Δα(14) Δβ(14) Δγ(15) Δi(12) ΔΩ(12)
    # ΔA(12) UDRA-dot(5) — starts at bit 72 (38 + 34).
    edc_start_bit = 72
    Δα = get_twos_complement_num(word, word_length, edc_start_bit + 8, 14) * 2.0^-34
    Δβ = get_twos_complement_num(word, word_length, edc_start_bit + 22, 14) * 2.0^-34
    Δγ =
        get_twos_complement_num(word, word_length, edc_start_bit + 36, 15) *
        2.0^-32 *
        GNSS_PI
    Δi =
        get_twos_complement_num(word, word_length, edc_start_bit + 51, 12) *
        2.0^-32 *
        GNSS_PI
    ΔΩ =
        get_twos_complement_num(word, word_length, edc_start_bit + 63, 12) *
        2.0^-32 *
        GNSS_PI
    ΔA = get_twos_complement_num(word, word_length, edc_start_bit + 75, 12) * 2.0^-9
    UDRA_dot = get_twos_complement_num(word, word_length, edc_start_bit + 87, 5)
    differential_correction = GPSL1C_DDifferentialCorrection(;
        PRN_ID,
        t_op_D,
        t_OD,
        dc_data_type,
        δa_f0,
        δa_f1,
        UDRA,
        UDRA_dot,
        Δα,
        Δβ,
        Δγ,
        Δi,
        ΔΩ,
        ΔA,
    )
    differential_corrections =
        writable_container(raw.differential_corrections, storage.differential_corrections)
    set!(differential_corrections, PRN_ID, differential_correction)
    GPSCNAV2Data(raw; differential_corrections)
end

# Character `k` (0-based) of a page-6 text message, and whether it is printable.
_text_char(word::UInt288, k::Int) = UInt8(get_bits(word, L1C_D_SF3_INFO_BITS, 19 + 8k, 8))
_is_printable(code::UInt8) = code >= 0x20 && code < 0x7f

# The `i`-th printable character of a page-6 text message, or `0x00` past the last.
function _text_printable_char(word::UInt288, i::Int)
    seen = 0
    for k = 0:(L1C_D_TEXT_LENGTH-1)
        code = _text_char(word, k)
        if _is_printable(code)
            seen += 1
            seen == i && return code
        end
    end
    return 0x00
end

"""
Subframe 3, page 6 — 29 ASCII characters at bits 19-250 (IS-GPS-800J Fig 3.5-7).
"""
function parse_sf3_page6(raw::GPSCNAV2Data, word::UInt288)
    # Built inline (a `CStaticString`, not a heap `String`): no allocation. Printable
    # ASCII is kept and NUL/control padding skipped, so each kept character moves to the
    # `i`-th printable position; the tail stays NUL, where a `CStaticString` ends.
    units = ntuple(i -> _text_printable_char(word, i), Val(L1C_D_TEXT_LENGTH))
    GPSCNAV2Data(raw; text_message = CStaticString(units))
end

"""
Subframe 3, page 8 — Integrity Support Message: PRN (bits 1-8), page number (9-14), the
236-bit ISM Packet at bits 15-250 and the CRC (IRN-IS-800J-004 Figure 3.5-8a).

IRN-IS-800J-004 (16 June 2026) replaces the page's GNSS ID through ISM CRC (IS-GPS-800J
Table 3.5-9, bits 15-250) with the ISM Packet of IS-GPS-200 §30.3.3.10, the one CNAV
message type 40 carries, so it is decoded by the same `gps_ism_packet` under the same
rules: a test packet, a reserved GNSS ID or an ISM CRC mismatch leaves `ism` as it was.

Only the amended layout is decoded. Unlike message type 40, the earlier page 8 cannot
be told from it: its ISM CRC (DO-246E) sits at bits 219-250 over bits 15-218, exactly
where the packet's ISM CRC sits and over the same bits, and only its GNSS ID, `WN_ISM`,
`TOW_ISM` and `t_correl` keep their place and meaning in the packet. If that CRC is the
packet's CRC-32Q, an earlier-layout page with GNSS ID `0100` passes and is misread from
`b_nom` on; the notice gives no format indicator to tell them apart.
"""
function parse_sf3_page8(raw::GPSCNAV2Data, word::UInt288)
    ism = gps_ism_packet(word, L1C_D_SF3_INFO_BITS, 15)
    isnothing(ism) ? raw : GPSCNAV2Data(raw; ism)
end

"""
    validate_data(state::GNSSDecoderState{<:GPSCNAV2Data})

Promote `raw_data` to `data` once subframe 2 is decoded. The BCH match happens in `try_sync`
and TOI continuity in [`decode_sync_unit`](@ref); this hook only publishes validated
data and arms the streaming counter.
"""
function validate_data(state::GNSSDecoderState{<:GPSCNAV2Data})
    if has_positioning_messages(state.raw_data)
        # `data` gets its own container copies (`publish_data`), as `raw_data`'s keep
        # being written by later pages.
        return GNSSDecoderState(
            state;
            data = publish_data(state.cache.storage, state.raw_data),
            num_bits_since_time_reference = state.constants.preamble_length,
        )
    end
    return state
end

"""
$(TYPEDSIGNATURES)

Check if the GPS L1C-D satellite is healthy and usable for positioning.

Healthy iff the 1-bit L1C health flag of subframe 2 is 0 (Signal OK; IS-GPS-800J
§3.5.3.4).

!!! warning

    Requires subframe 2 to have been decoded; returns `false` until then.

# Arguments

  - `state::GNSSDecoderState{<:GPSCNAV2Data}`: GPS L1C-D decoder state.

# Returns

  - `Bool`: `true` iff the L1C signal-health bit indicates OK.
"""
function is_sat_healthy(state::GNSSDecoderState{<:GPSCNAV2Data})
    state.data.l1c_health === false
end
