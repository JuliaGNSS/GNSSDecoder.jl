# BeiDou B1C (B-CNAV1) decoder — BDS-SIS-ICD-B1C-1.0 (2017-12).
#
# The B1C data component (`BeiDouB1C_D`) carries the B-CNAV1 message as
# 18-second frames of 1800 symbols at 100 sps (ICD §6.2.1):
#
#   - Subframe 1: 72 symbols — the 6-bit PRN encoded BCH(21,6) followed by
#     the 8-bit Seconds-Of-Hour (SOH) count encoded BCH(51,8) (ICD §6.2.2.1,
#     §6.2.3.1). Carries system time as well as framing: SOH is a BDT parameter
#     counting 18 s units within the hour (§7.3, Table 7-2, range 0..3582 s), and
#     the decoder uses it for both — the sync match and the time stamp, exactly
#     like GPS L1C-D's TOI.
#   - Subframes 2 + 3: 1728 symbols, block-interleaved (36×48 with a
#     staggered write, ICD §6.2.2.4) over a 1200-symbol 64-ary LDPC(200,100)
#     codeword for subframe 2 (600 info bits) and a 528-symbol 64-ary
#     LDPC(88,44) codeword for subframe 3 (264 info bits). Decoded through
#     the binary-image `.alist` matrices (see `scripts/generate_beidou_alist.jl`)
#     with the shared `load_ldpc_decoder`/`ldpc_decode_word` helpers.
#
# Subframe 2 carries WN/HOW, IODC/IODE, the full ephemeris, the clock
# correction, and the group-delay parameters (ICD Figure 6-6); subframe 3 is
# paged (PageID 1-4: iono + BDT-UTC, reduced almanacs, EOP + BGTO, midi
# almanac; ICD Figures 6-8 .. 6-11). Every layout and scale below was
# transcribed from the ICD figures/tables cited inline.
#
# Frame sync follows the generic streaming framework in `src/gnss.jl`: the
# soft-symbol `CircularDeque` is sized to one full frame (1800) plus the
# 72-symbol subframe-1 segment of the *next* frame (1872 total) so the BCH
# match can be confirmed at both ends of the window before locking — the
# same two-subframe pattern as GPS L1C-D, with SOH in the role of TOI.

# ---- Subframe-1 BCH codewords (ICD §6.2.2.1, Table 6-1, Figure 6-2) ---------
#
# Subframe 1 encodes its 14 info bits as two concatenated cyclic codewords:
# the 6 PRN bits into 21 symbols and the 8 SOH bits into 51 symbols (21 + 51
# = 72; transmitted PRN-codeword first, matching the subframe-1 MSB-first
# layout of Figure 6-5 — cross-checked against PocketSDR's `sync_BCNV1_frame`).
#
# Both codewords come from `bch_lfsr_codeword` (`src/bch_toi.jl`), the same
# Fibonacci-LFSR construction the GPS L1C-D TOI table uses — only the symbol
# count, register width and tap mask differ.
#
# The ICD's "generator polynomials" g21,6(x) = x⁶+x⁴+x²+x+1 and
# g51,8(x) = x⁸+x⁷+x⁴+x³+x²+x+1 have degree k (not n-k): they are the
# *characteristic polynomials* of the k-stage encoder registers of Figure
# 6-2. Each codeword is the k-stage Fibonacci LFSR output clocked n times
# from an initial state holding the info value: seeding with the
# *bit-reversed* info value and emitting the register LSB makes the code
# systematic — the first k emitted symbols spell the info value MSB-first —
# which is the same construction (and the same verified convention) as the
# GPS L1C-D TOI table in `src/bch_toi.jl`, and matches PocketSDR's
# `LFSR(n, rev_reg(v, k), taps, k)`. The tap masks below are the
# non-leading coefficients of the ICD polynomials, which is what
# `bch_lfsr_codeword` takes.

const B1C_PRN_BCH_WIDTH = 6
const B1C_PRN_BCH_TAP = 0b010111       # g21,6 = x⁶ + x⁴ + x² + x + 1
const B1C_PRN_CODEWORD_LEN = 21
const B1C_SOH_BCH_WIDTH = 8
const B1C_SOH_BCH_TAP = 0b10011111     # g51,8 = x⁸ + x⁷ + x⁴ + x³ + x² + x + 1
const B1C_SOH_CODEWORD_LEN = 51
const B1C_SOH_RANGE = 200              # SOH counts 18 s units within the hour (ICD §7.3)
const B1C_SOH_SECONDS = 18             # seconds per SOH step, i.e. per B-CNAV1 frame

const B1C_PRN_MASK21 = (UInt64(1) << B1C_PRN_CODEWORD_LEN) - UInt64(1)
const B1C_SOH_MASK51 = (UInt64(1) << B1C_SOH_CODEWORD_LEN) - UInt64(1)

# Correction radii for the two subframe-1 BCH codes. Both are the classical
# bounded-distance radius `t = floor((d_min - 1) / 2)`, and both minimum
# distances are properties of the codebooks this file generates rather than
# numbers taken on trust — `test/beidou_b1c.jl` recomputes them exhaustively.
#
#   BCH(21,6) over the 63 assigned PRNs: d_min = 8 -> t = 3
#   BCH(51,8) over the 200 assigned SOH counts: d_min = 24 -> t = 11
#
# Matching a codeword exactly, as this decoder used to, throws that away: at
# 28 dB-Hz a 72-symbol subframe 1 arrives error-free about 6 % of the time and
# within the correction radius about 99 % of it, which is roughly 6 dB of
# acquisition sensitivity on a field the ICD protected precisely so it would
# survive the margin. The cost is a false-sync rate that rises from
# unmeasurable to still unmeasurable: `try_sync` demands two windows a frame
# apart carrying consecutive SOH counts, which puts a random-symbol false sync
# at ~1e-13 per window pair, well under one per thousand years at 100 sps —
# and the LDPC CRC on subframes 2 and 3 discards it in any case.
const B1C_PRN_MAX_ERRORS = 3
const B1C_SOH_MAX_ERRORS = 11

"""
21-symbol BCH(21,6) codeword for a PRN (1..63), bit 0 = first transmitted symbol.
"""
b1c_prn_codeword(prn::Int) =
    bch_lfsr_codeword(prn, B1C_PRN_CODEWORD_LEN, B1C_PRN_BCH_WIDTH, B1C_PRN_BCH_TAP)

"""
51-symbol BCH(51,8) codeword for an SOH count (0..199), bit 0 = first transmitted symbol.
"""
b1c_soh_codeword(soh::Int) =
    bch_lfsr_codeword(soh, B1C_SOH_CODEWORD_LEN, B1C_SOH_BCH_WIDTH, B1C_SOH_BCH_TAP)

# 200-entry SOH codeword table, materialised once at precompile time.
const B1C_SOH_CODEWORDS = ntuple(i -> b1c_soh_codeword(i - 1), B1C_SOH_RANGE)

# ---- Frame geometry (ICD §6.2.1, §6.2.2.4) ----------------------------------

const B1C_FRAME_LENGTH = 1800
const B1C_SUBFRAME1_LENGTH = 72
const B1C_WINDOW_LENGTH = B1C_FRAME_LENGTH + B1C_SUBFRAME1_LENGTH  # 1872

const B1C_SF2_SYMBOLS = 1200
const B1C_SF3_SYMBOLS = 528
const B1C_PAYLOAD_SYMBOLS = B1C_SF2_SYMBOLS + B1C_SF3_SYMBOLS  # 1728
const B1C_INTERLEAVER_ROWS = 36
const B1C_INTERLEAVER_COLS = 48

const B1C_SF2_INFO_BITS = 600
const B1C_SF3_INFO_BITS = 264

# Staggered interleaver row assignment (ICD §6.2.2.4): the transmitter writes
# rows 1,2 from subframe 2, row 3 from subframe 3, repeating until the 528
# subframe-3 symbols fill 11 rows (row 33), then the last 3 rows with the
# remaining 144 subframe-2 symbols; it reads the 36×48 array out by column.
# After the receiver undoes the column read (`deinterleave(…, 36, 48)` gives
# the array back in row-major order) these constants pick the rows back
# apart: subframe 2 = rows {3k+1, 3k+2 : k=0..10} ∪ {34,35,36} in that
# order, subframe 3 = rows {3k+3 : k=0..10}.
const B1C_SF2_ROW_ORDER = (vcat([[3k + 1, 3k + 2] for k = 0:10]...)..., 34, 35, 36)
const B1C_SF3_ROW_ORDER = ntuple(k -> 3k, 11)

# ---- Physical constants and orbit references (ICD §7.7) ---------------------

# Semi-major axis reference values (ICD Table 7-8 footnote ***, meters).
const B1C_A_REF_MEO = 27_906_100.0
const B1C_A_REF_IGSO_GEO = 42_162_200.0

# Reduced-almanac reference values (ICD Table 7-16 footnote ****): e = 0,
# δi = 0 relative to i = 55° (MEO/IGSO) or i = 0° (GEO).

"""
    BeiDouB1CConstants

BDCS constants and B-CNAV1 message structure parameters for BeiDou B1C
decoding.

The frame is modelled through the generic streaming framework:
`preamble_length` is the 72-symbol subframe-1 BCH segment of the *next* frame
retained at the tail of the sync window, and `syncro_sequence_length` is the
1800-symbol frame that is drained once a frame is decoded.

# Fields

$(TYPEDFIELDS)

# Reference

BDS-SIS-ICD-B1C-1.0, Sections 6.2 and 7.5-7.7 (Tables 7-8, 7-9).
"""
Base.@kwdef struct BeiDouB1CConstants <: AbstractGNSSConstants
    """
    Frame length drained after each decoded frame (1800 symbols)
    """
    syncro_sequence_length::Int = B1C_FRAME_LENGTH
    """
    Trailing next-frame subframe-1 BCH segment retained for sync (72 symbols)
    """
    preamble_length::Int = B1C_SUBFRAME1_LENGTH
    """
    Mathematical constant π (BDS-SIS-ICD-B1C-1.0 Table 7-9)
    """
    PI::Float64 = GNSS_PI
    """
    BDCS Earth rotation rate (rad/s) — differs from the WGS-84 value
    """
    Ω_dot_e::Float64 = BEIDOU_EARTH_ROTATION_RATE
    """
    Speed of light (m/s)
    """
    c::Float64 = SPEED_OF_LIGHT
    """
    BDCS Earth gravitational parameter (m³/s²)
    """
    μ::Float64 = BEIDOU_μ
    """
    Relativistic correction constant F = -2√μ/c² (s/√m, ICD §7.5.2)
    """
    F::Float64 = -4.442807309e-10
end

"""
    BeiDouB1CBGTO

One BDT-GNSS time offset (BGTO) parameter set from B-CNAV1 subframe 3, page
type 3 (BDS-SIS-ICD-B1C-1.0 Figure 6-20, Table 7-21).

`Δt = t_BD - t_GNSS = A_0BGTO + A_1BGTO·Δτ + A_2BGTO·Δτ²` with
`Δτ = t_BD - t_0BGTO + 604800(WN - WN_0BGTO)` (ICD Eq. 7-30). Different
frames may broadcast offsets for different systems, so sets are keyed by
`GNSS_ID` in [`BeiDouB1CData`](@ref).

# Fields

  - `GNSS_ID::Int`: GNSS type the offset refers to (1 GPS, 2 Galileo,
    3 GLONASS; 0 marks the parameters as unavailable and is never stored).
  - `WN_0BGTO::Int`: Reference week number.
  - `t_0BGTO::Int`: Reference time of week (seconds).
  - `A_0BGTO::Float64`, `A_1BGTO::Float64`, `A_2BGTO::Float64`: Bias /
    drift / drift-rate coefficients (s, s/s, s/s²).

# Reference

BDS-SIS-ICD-B1C-1.0, Figure 6-20, Table 7-21, §7.13.2.
"""
Base.@kwdef struct BeiDouB1CBGTO
    GNSS_ID::Int
    WN_0BGTO::Int
    t_0BGTO::Int
    A_0BGTO::Float64
    A_1BGTO::Float64
    A_2BGTO::Float64
end

"""
    BeiDouB1CData

Decoded BeiDou B1C (B-CNAV1) navigation message data.

Holds the subframe-2 system time, ephemeris, clock, and group-delay
parameters (BDS-SIS-ICD-B1C-1.0 Figure 6-6, Tables 7-5 .. 7-8) plus the
paged subframe-3 contents (Figures 6-8 .. 6-11). Semi-circle quantities are
converted to radians on decode (multiplied by π); all fields are
`Union{Nothing,…}` until first decoded.

# Sync / timing (ICD §7.3)

  - `soh::Int`: Last validated Seconds-Of-Hour count (0..199, in 18 s units);
    the epoch it denotes is the leading edge of the *current* frame's
    subframe 1. Seconds of week at that epoch = `HOW·3600 + soh·18`.
  - `HOW::Int64`: Hours of week (0..167, subframe 2).
  - `WN::Int64`: BDT week number (0..8191, subframe 2).

# Issue of data (ICD §7.4)

  - `IODC::Int64`: Issue of data, clock (10 bits).
  - `IODE::Int64`: Issue of data, ephemeris (8 bits). Not necessarily equal to
    the 8 LSBs of `IODC`, even inside one CRC-protected subframe 2: §7.4.3 says
    the two "may be different … during the update of the ephemeris and clock
    correction data", and a user "shall use the preceding matched pair" until
    they agree. `is_decoding_completed_for_positioning` enforces that match.

# Ephemeris (Figure 6-12/6-13, Table 7-8)

  - `t_0e::Int64`: Ephemeris reference time of week (seconds; the ICD writes `toe`).
  - `sat_type::Int64`: Satellite orbit type (raw 2 bits: 1 GEO, 2 IGSO,
    3 MEO, 0 reserved). Selects the semi-major-axis reference
    `A_ref = 27 906 100 m` (MEO) or `42 162 200 m` (IGSO/GEO).
  - `ΔA::Float64`: Semi-major axis difference at reference time (meters).
  - `A_dot::Float64`: Change rate in semi-major axis (m/s).
  - `Δn_0::Float64`: Mean motion difference at reference time (rad/s).
  - `Δn_0_dot::Float64`: Rate of mean motion difference (rad/s²).
  - `M_0::Float64`: Mean anomaly at reference time (rad).
  - `e::Float64`: Eccentricity (dimensionless).
  - `ω::Float64`: Argument of perigee (rad).
  - `Ω_0::Float64`: Longitude of ascending node at weekly epoch (rad).
  - `i_0::Float64`: Inclination angle at reference time (rad).
  - `Ω_dot::Float64`: Rate of right ascension (rad/s).
  - `i_dot::Float64`: Rate of inclination angle (rad/s).
  - `C_is::Float64`, `C_ic::Float64`: Sine/cosine inclination harmonic corrections (rad).
  - `C_rs::Float64`, `C_rc::Float64`: Sine/cosine orbit-radius harmonic corrections (m).
  - `C_us::Float64`, `C_uc::Float64`: Sine/cosine argument-of-latitude harmonic corrections (rad).

# Clock and group delay (Figure 6-14, Tables 7-5, 7-6)

  - `t_0c::Int64`: Clock reference time of week (seconds; the ICD writes `toc`).

  - `a_f0::Float64`, `a_f1::Float64`, `a_f2::Float64`: Clock bias / drift /
    drift-rate (s, s/s, s/s²; the ICD names them `a_0`/`a_1`/`a_2`). The
    broadcast clock is referenced to the B3I signal (ICD §7.6.1).

  - `T_GD_B2ap::Float64`: Group delay differential of the B2a pilot (s).

  - `ISC_B1Cd::Float64`: Group delay differential between the B1C data and
    pilot components (s).

  - `T_GD_B1Cp::Float64`: Group delay differential of the B1C pilot, relative to
    B3I (s).

    A B1C receiver needs both: `a_f0` is referenced to B3I (§7.6.1), so reaching
    the B1C data component means adding `T_GD_B1Cp` (B1C pilot vs B3I) and then
    `ISC_B1Cd` (B1C data vs B1C pilot). They ride in the same CRC-protected
    subframe 2 as the clock, so `is_decoding_completed_for_positioning` requires
    both and a consumer never has to treat either as optional. `T_GD_B2ap` is
    the B2a pilot's delta and is not required here.

# Subframe 3 page header (Figures 6-8 .. 6-11, §7.14-7.17)

Present on every defined page type, refreshed whenever any page decodes:

  - `HS::Int64`: Satellite health status (Table 7-22: 0 healthy,
    1 unhealthy or in test, 2-3 reserved).
  - `DIF::Bool`, `SIF::Bool`, `AIF::Bool`: Data / signal / accuracy
    integrity flags for B1C (Table 7-23; `true` flags a problem).
  - `SISMAI::Int64`: Signal-in-space monitoring accuracy index (4 bits).
  - `SISAI_oe::Int64`: Orbit along-track/cross-track accuracy index
    (5 bits, pages 1 and 3).
  - `t_op::Int64`: Time of week for data prediction (raw 11-bit count; B1C-1.0
    §7.16 defers the SISAI *definitions* and gives no parameter table at all, so
    no scale factor is applied), and
    `SISAI_ocb::Int64`, `SISAI_oc1::Int64`, `SISAI_oc2::Int64`: orbit
    radius / clock accuracy indices (pages 1, 2, and 4).

# Subframe 3, page type 1 — ionosphere + BDT-UTC (Figures 6-16, 6-17)

  - `α_bdgim_1 .. α_bdgim_9::Float64`: BDGIM ionospheric coefficients (TECu; Table 7-10 —
    the broadcast `α_bdgim_5` carries scale factor −2⁻³, applied on decode). The ICD
    names them `α_1` … `α_9`; the `bdgim` qualifier is this package's, because the
    Klobuchar `α_0` … `α_3` on `GPSL1CAData`, `GPSCNAVData`, `GPSL1C_DData` and
    `BeiDouDNAVData` are a different model with different units, and one name for
    both would let generic code apply the wrong one without complaint.
  - `A_0UTC,A_1UTC,A_2UTC::Float64`: BDT-UTC polynomial (s, s/s, s/s²).
  - `Δt_LS,Δt_LSF::Int64`: current *or past* and future leap-second counts (s;
    Table 7-20 — the past case is what selects Eq. (7-29) over (7-25)).
  - `t_0t::Int64`: UTC reference time of week (s); `WN_0t::Int64`: reference
    week. (The ICD writes `t_ot`/`WN_ot`.)
    (The ICD writes `t_ot`/`WN_ot`.)
  - `WN_LSF::Int64`, `DN::Int64`: leap-second reference week/day (day 0..6).

# Subframe 3, page type 3 — EOP + BGTO (Figures 6-19, 6-20)

  - `t_EOP::Int64`: EOP reference time of week (s).
  - `PM_X,PM_X_dot,PM_Y,PM_Y_dot::Float64`: polar motion (arc-seconds, arc-seconds/day).
  - `ΔUT1::Float64`, `ΔUT1_dot::Float64`: UT1-UTC difference (s) and rate (s/day).
  - `bgtos::Dictionary{Int,BeiDouB1CBGTO}`: BDT-GNSS time offsets keyed by GNSS ID.

# Subframe 3, page types 2/4 — keyed almanac dictionaries

  - `reduced_almanacs::Dictionary{Int,BeiDouReducedAlmanac}` (page type 2).
  - `midi_almanacs::Dictionary{Int,BeiDouMidiAlmanac}` (page type 4).

# Counters

  - `num_sf3_pages_received::Int`: Count of CRC-valid subframe-3 pages received.

# Reference

BDS-SIS-ICD-B1C-1.0, Figures 6-5 .. 6-21 and Tables 7-2 .. 7-23.
"""
Base.@kwdef struct BeiDouB1CData <: AbstractBeiDouData
    soh::Union{Nothing,Int} = nothing
    HOW::Union{Nothing,Int64} = nothing
    WN::Union{Nothing,Int64} = nothing

    IODC::Union{Nothing,Int64} = nothing
    IODE::Union{Nothing,Int64} = nothing

    t_0e::Union{Nothing,Int64} = nothing
    sat_type::Union{Nothing,Int64} = nothing
    ΔA::Union{Nothing,Float64} = nothing
    A_dot::Union{Nothing,Float64} = nothing
    Δn_0::Union{Nothing,Float64} = nothing
    Δn_0_dot::Union{Nothing,Float64} = nothing
    M_0::Union{Nothing,Float64} = nothing
    e::Union{Nothing,Float64} = nothing
    ω::Union{Nothing,Float64} = nothing
    Ω_0::Union{Nothing,Float64} = nothing
    i_0::Union{Nothing,Float64} = nothing
    Ω_dot::Union{Nothing,Float64} = nothing
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
    T_GD_B2ap::Union{Nothing,Float64} = nothing
    ISC_B1Cd::Union{Nothing,Float64} = nothing
    T_GD_B1Cp::Union{Nothing,Float64} = nothing

    # --- Subframe 3 page header (every defined page type) ---
    HS::Union{Nothing,Int64} = nothing
    DIF::Union{Nothing,Bool} = nothing
    SIF::Union{Nothing,Bool} = nothing
    AIF::Union{Nothing,Bool} = nothing
    SISMAI::Union{Nothing,Int64} = nothing
    SISAI_oe::Union{Nothing,Int64} = nothing
    t_op::Union{Nothing,Int64} = nothing
    SISAI_ocb::Union{Nothing,Int64} = nothing
    SISAI_oc1::Union{Nothing,Int64} = nothing
    SISAI_oc2::Union{Nothing,Int64} = nothing

    # --- Subframe 3, page type 1: BDGIM iono + BDT-UTC ---
    α_bdgim_1::Union{Nothing,Float64} = nothing
    α_bdgim_2::Union{Nothing,Float64} = nothing
    α_bdgim_3::Union{Nothing,Float64} = nothing
    α_bdgim_4::Union{Nothing,Float64} = nothing
    α_bdgim_5::Union{Nothing,Float64} = nothing
    α_bdgim_6::Union{Nothing,Float64} = nothing
    α_bdgim_7::Union{Nothing,Float64} = nothing
    α_bdgim_8::Union{Nothing,Float64} = nothing
    α_bdgim_9::Union{Nothing,Float64} = nothing
    A_0UTC::Union{Nothing,Float64} = nothing
    A_1UTC::Union{Nothing,Float64} = nothing
    A_2UTC::Union{Nothing,Float64} = nothing
    Δt_LS::Union{Nothing,Int64} = nothing
    t_0t::Union{Nothing,Int64} = nothing
    WN_0t::Union{Nothing,Int64} = nothing
    WN_LSF::Union{Nothing,Int64} = nothing
    DN::Union{Nothing,Int64} = nothing
    Δt_LSF::Union{Nothing,Int64} = nothing

    # --- Subframe 3, page type 3: EOP + BGTO ---
    t_EOP::Union{Nothing,Int64} = nothing
    PM_X::Union{Nothing,Float64} = nothing
    PM_X_dot::Union{Nothing,Float64} = nothing
    PM_Y::Union{Nothing,Float64} = nothing
    PM_Y_dot::Union{Nothing,Float64} = nothing
    ΔUT1::Union{Nothing,Float64} = nothing
    ΔUT1_dot::Union{Nothing,Float64} = nothing
    bgtos::Union{Nothing,Dictionary{Int,BeiDouB1CBGTO}} = nothing

    # --- Subframe 3, page types 2/4: per-SV keyed dictionaries ---
    reduced_almanacs::Union{Nothing,Dictionary{Int,BeiDouReducedAlmanac}} = nothing
    midi_almanacs::Union{Nothing,Dictionary{Int,BeiDouMidiAlmanac}} = nothing

    num_sf3_pages_received::Int = 0
end

function BeiDouB1CData(
    data::BeiDouB1CData;
    soh = data.soh,
    HOW = data.HOW,
    WN = data.WN,
    IODC = data.IODC,
    IODE = data.IODE,
    t_0e = data.t_0e,
    sat_type = data.sat_type,
    ΔA = data.ΔA,
    A_dot = data.A_dot,
    Δn_0 = data.Δn_0,
    Δn_0_dot = data.Δn_0_dot,
    M_0 = data.M_0,
    e = data.e,
    ω = data.ω,
    Ω_0 = data.Ω_0,
    i_0 = data.i_0,
    Ω_dot = data.Ω_dot,
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
    T_GD_B2ap = data.T_GD_B2ap,
    ISC_B1Cd = data.ISC_B1Cd,
    T_GD_B1Cp = data.T_GD_B1Cp,
    HS = data.HS,
    DIF = data.DIF,
    SIF = data.SIF,
    AIF = data.AIF,
    SISMAI = data.SISMAI,
    SISAI_oe = data.SISAI_oe,
    t_op = data.t_op,
    SISAI_ocb = data.SISAI_ocb,
    SISAI_oc1 = data.SISAI_oc1,
    SISAI_oc2 = data.SISAI_oc2,
    α_bdgim_1 = data.α_bdgim_1,
    α_bdgim_2 = data.α_bdgim_2,
    α_bdgim_3 = data.α_bdgim_3,
    α_bdgim_4 = data.α_bdgim_4,
    α_bdgim_5 = data.α_bdgim_5,
    α_bdgim_6 = data.α_bdgim_6,
    α_bdgim_7 = data.α_bdgim_7,
    α_bdgim_8 = data.α_bdgim_8,
    α_bdgim_9 = data.α_bdgim_9,
    A_0UTC = data.A_0UTC,
    A_1UTC = data.A_1UTC,
    A_2UTC = data.A_2UTC,
    Δt_LS = data.Δt_LS,
    t_0t = data.t_0t,
    WN_0t = data.WN_0t,
    WN_LSF = data.WN_LSF,
    DN = data.DN,
    Δt_LSF = data.Δt_LSF,
    t_EOP = data.t_EOP,
    PM_X = data.PM_X,
    PM_X_dot = data.PM_X_dot,
    PM_Y = data.PM_Y,
    PM_Y_dot = data.PM_Y_dot,
    ΔUT1 = data.ΔUT1,
    ΔUT1_dot = data.ΔUT1_dot,
    bgtos = data.bgtos,
    reduced_almanacs = data.reduced_almanacs,
    midi_almanacs = data.midi_almanacs,
    num_sf3_pages_received = data.num_sf3_pages_received,
)
    BeiDouB1CData(
        soh,
        HOW,
        WN,
        IODC,
        IODE,
        t_0e,
        sat_type,
        ΔA,
        A_dot,
        Δn_0,
        Δn_0_dot,
        M_0,
        e,
        ω,
        Ω_0,
        i_0,
        Ω_dot,
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
        T_GD_B2ap,
        ISC_B1Cd,
        T_GD_B1Cp,
        HS,
        DIF,
        SIF,
        AIF,
        SISMAI,
        SISAI_oe,
        t_op,
        SISAI_ocb,
        SISAI_oc1,
        SISAI_oc2,
        α_bdgim_1,
        α_bdgim_2,
        α_bdgim_3,
        α_bdgim_4,
        α_bdgim_5,
        α_bdgim_6,
        α_bdgim_7,
        α_bdgim_8,
        α_bdgim_9,
        A_0UTC,
        A_1UTC,
        A_2UTC,
        Δt_LS,
        t_0t,
        WN_0t,
        WN_LSF,
        DN,
        Δt_LSF,
        t_EOP,
        PM_X,
        PM_X_dot,
        PM_Y,
        PM_Y_dot,
        ΔUT1,
        ΔUT1_dot,
        bgtos,
        reduced_almanacs,
        midi_almanacs,
        num_sf3_pages_received,
    )
end

# The default struct `==` falls back to `===` for the mutable `Dictionary`
# fields; compare field-by-field (mirrors `GPSL1C_DData`).
Base.:(==)(a::BeiDouB1CData, b::BeiDouB1CData) = fields_equal(a, b)

"""
$(TYPEDEF)

Per-decoder cache for the BeiDou B1C signal.

Holds the soft-symbol `CircularDeque{Float32}` (capacity 1872 symbols: one
1800-symbol frame plus the next frame's 72-symbol subframe 1) and the two Aff3ct
LDPC belief-propagation decoders for the binary images of the 64-ary
LDPC(200,100) / LDPC(88,44) codes (subframe 2: K=600, N=1200; subframe 3: K=264,
N=528 — bit counts of the `data/bcnv1_sf{2,3}.alist` matrices). The decoder
objects are mutable Aff3ct handles; they are shared by reference through the
otherwise-immutable [`GNSSDecoderState`](@ref).

The rest is scratch, sized once so that a decoded frame allocates nothing: this
satellite's PRN codeword, which `try_sync` compares against on every symbol, and
the four symbol buffers the deinterleave and the two LDPC decodes run through.

# Fields

$(TYPEDFIELDS)
"""
struct BeiDouB1CCache <: AbstractGNSSCache
    """
    Soft-symbol buffer (1872 = 1800 frame + 72 next subframe-1)
    """
    soft_buffer::CircularDeque{Float32}
    """
    Subframe-2 LDPC decoder and its scratch buffers (binary image, K=600, N=1200)
    """
    sf2_ldpc::LDPCScratch
    """
    Subframe-3 LDPC decoder and its scratch buffers (binary image, K=264, N=528)
    """
    sf3_ldpc::LDPCScratch
    """
    This satellite's 21-symbol BCH(21,6) PRN codeword, built once — `try_sync`
    compares against it on every symbol
    """
    prn_codeword::UInt64
    """
    Polarity-resolved 1728-symbol interleaved SF2+SF3 payload, copied out per frame
    """
    payload::Vector{Float32}
    """
    The 36×48 interleaver array in row-major order, i.e. `payload` deinterleaved
    """
    array_rows::Vector{Float32}
    """
    The subframe-2 rows of `array_rows`, concatenated into its 1200-symbol codeword
    """
    sf2_symbols::Vector{Float32}
    """
    The subframe-3 rows of `array_rows`, concatenated into its 528-symbol codeword
    """
    sf3_symbols::Vector{Float32}
end

function BeiDouB1CCache(prn::Int)
    BeiDouB1CCache(
        CircularDeque{Float32}(B1C_WINDOW_LENGTH),
        LDPCScratch(alist_path("bcnv1_sf2.alist")),
        LDPCScratch(alist_path("bcnv1_sf3.alist")),
        b1c_prn_codeword(prn),
        Vector{Float32}(undef, B1C_PAYLOAD_SYMBOLS),
        Vector{Float32}(undef, B1C_PAYLOAD_SYMBOLS),
        Vector{Float32}(undef, B1C_SF2_SYMBOLS),
        Vector{Float32}(undef, B1C_SF3_SYMBOLS),
    )
end

# The LDPC decoder handles are stateless w.r.t. equality; two B1C caches are
# equal when their soft buffers match (mirrors `GPSL1C_DCache`).
function Base.:(==)(a::BeiDouB1CCache, b::BeiDouB1CCache)
    deques_equal(a.soft_buffer, b.soft_buffer)
end

function is_subframe2_decoded(data::BeiDouB1CData)
    !isnothing(data.WN) &&
        !isnothing(data.HOW) &&
        !isnothing(data.IODE) &&
        !isnothing(data.t_0e) &&
        is_known_sat_type(data.sat_type) &&
        !isnothing(data.ΔA) &&
        !isnothing(data.M_0) &&
        !isnothing(data.e) &&
        !isnothing(data.ω) &&
        !isnothing(data.Ω_0) &&
        !isnothing(data.i_0) &&
        !isnothing(data.a_f0) &&
        !isnothing(data.a_f1)
end

# Positioning readiness: a validated SOH (time), the subframe-2 ephemeris +
# clock set, and the subframe-3 health status `HS` — so `is_sat_healthy` is
# guaranteed decodable whenever this is `true` (see
# `is_decoding_completed_for_positioning` in src/gnss.jl).
#
# There is no cross-subframe IOD *stitching* to do — IODE and IODC both live
# in the CRC-protected subframe 2 — but they must still be a matched pair.
# ICD §7.4.3: "The IODE value received by the user may be different from the
# 8 LSBs of IODC during the update of the ephemeris and clock correction
# data … the user shall use the preceding matched pair … until the updated
# IODE and the 8 LSBs of IODC are the same." So the ephemeris and the clock
# set can disagree inside one valid block, and the same gate `b2a.jl` applies
# is applied here.
#
# THE GROUP DELAYS ARE REQUIRED HERE, unlike on every other signal in this
# package. `is_decoding_completed_for_positioning` (src/gnss.jl) excludes group
# delays *that would mean waiting for a message the ICD does not schedule*; the
# single-band correction is included wherever it rides in the required set
# already. On B1C it does: `T_GD_B1Cp` and `ISC_B1Cd` sit at bits 558-569 and
# 546-557 of the same 600-bit, CRC-24Q-protected subframe 2 as the ephemeris and
# the clock, decoded by the same constructor call, and subframe 2 is in every
# 18-second B-CNAV1 frame. Requiring them therefore costs a B1C receiver exactly
# zero extra time to first fix.
#
# And they are genuinely required, not a refinement. ICD §7.6: the clock
# parameter a0 carries "the equipment group delay of the B3I signal … the
# reference equipment group delay for the B1C signal". A B1C receiver is not
# using B3I, so it must add `T_GD_B1Cp` (B1C pilot vs B3I) and then `ISC_B1Cd`
# (B1C data vs B1C pilot) to reach its own component. Leaving them out is a
# metre-scale bias, not a refinement, and a consumer that treats `nothing` as
# zero — which is what the general policy tells it to do — would silently eat it.
#
# `T_GD_B2ap` is deliberately *not* required: it is the B2a pilot's delta and
# means nothing to a B1C-only user. Same rule, other side of it.
# Orbit class from the satellite's own broadcast `sat_type` (Ephemeris I,
# subframe 2); `nothing` until subframe 2 is decoded, and for the reserved code.
get_orbit_class(state::GNSSDecoderState{<:BeiDouB1CData}) =
    beidou_orbit_class(state.data.sat_type)

"""
$(TYPEDSIGNATURES)

Seconds of week of the epoch the SOH stamps: `HOW·3600 + soh·18`.

B-CNAV1 splits the time of week across two subframes and two rates — the
subframe-2 `HOW` counts whole hours of the week and the subframe-1 `soh` counts
the 18-second frames within the hour (BDS-SIS-ICD-B1C-1.0 §7.3, Table 7-2) — so
neither field is a time of week on its own. Both must be present; `nothing`
until they are.

The epoch is the leading edge of the *current* frame's subframe 1, which is what
`num_bits_after_valid_syncro_sequence` is anchored to for this signal (see
`validate_data` below); GPS L1C-D's TOI stamps the next frame instead, so its
otherwise identical-looking reconstruction is anchored differently.
"""
function get_time_of_week(data::BeiDouB1CData)
    (isnothing(data.HOW) || isnothing(data.soh)) && return nothing
    return Int(data.HOW) * 3600 + Int(data.soh) * B1C_SOH_SECONDS
end

"""
$(TYPEDSIGNATURES)

The BDT-to-`target` offset from subframe 3 page type 3, or `nothing` when this
satellite has not broadcast one for `target`.

B1C is the one shape that keys its offsets by system: successive frames may
carry sets for different systems, and each is filed under its `GNSS_ID` in
`data.bgtos` rather than overwriting the last (see [`BeiDouB1CBGTO`](@ref)), so
every system this satellite has ever broadcast an offset for stays answerable.
Sets tagged `GNSS_ID == 0` ("not available", §7.13.1) are dropped at decode and
never appear here.
"""
function get_time_offset(state::GNSSDecoderState{<:BeiDouB1CData}, target::TimeSystem)
    bgtos = state.data.bgtos
    isnothing(bgtos) && return nothing
    for (GNSS_ID, bgto) in pairs(bgtos)
        beidou_bgto_target(GNSS_ID) === target || continue
        return GNSSTimeOffset(
            target,
            bgto.A_0BGTO,
            bgto.A_1BGTO,
            bgto.A_2BGTO,
            bgto.t_0BGTO,
            bgto.WN_0BGTO,
        )
    end
    return nothing
end

function is_decoding_completed_for_positioning(data::BeiDouB1CData)
    !isnothing(data.soh) &&
        !isnothing(data.HS) &&
        is_subframe2_decoded(data) &&
        !isnothing(data.T_GD_B1Cp) &&
        !isnothing(data.ISC_B1Cd) &&
        !isnothing(data.IODC) &&
        data.IODE == data.IODC & 0xff
end

"""
$(TYPEDSIGNATURES)

Create a decoder state for BeiDou B1C (B-CNAV1) navigation messages.

Wires up a [`GNSSDecoderState`](@ref) with a 1872-symbol soft-symbol buffer,
the BCH(21,6) PRN / BCH(51,8) SOH subframe-1 codeword tables, and two Aff3ct
LDPC belief-propagation decoders loaded from the committed binary-image
`.alist` parity matrices in `data/` (see `scripts/generate_beidou_alist.jl`).

Like GPS L1C-D, the LDPC decode is flooding sum-product and therefore
scale-sensitive: feed soft symbols whose magnitudes are confidence-weighted
on a roughly LLR-like scale (`≈ 2·r/σ²`) for best performance at marginal
SNR (see the soft-symbol convention note on [`decode`](@ref)).

# Arguments

  - `prn::Int`: Pseudo-Random Noise code identifier (1-63 for BeiDou).

# Returns

  - `GNSSDecoderState{BeiDouB1CData}`: Initialized decoder state for BeiDou B1C.

# Example

```julia
state = BeiDouB1CDecoderState(30)          # PRN 30
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
function BeiDouB1CDecoderState(prn)
    1 <= prn <= 63 || throw(ArgumentError("BeiDou PRN must be in 1..63"))
    GNSSDecoderState(
        prn,
        BeiDouB1CData(),
        BeiDouB1CData(),
        BeiDouB1CConstants(),
        BeiDouB1CCache(prn),
        nothing,
        false,
    )
end

# Dispatch from the GNSSSignals type. B-CNAV1 rides on the B1C data
# component — `BeiDouB1C_D` — while `BeiDouB1C_P` is the dataless pilot.
function GNSSDecoderState(system::BeiDouB1C_D, prn)
    BeiDouB1CDecoderState(prn)
end

# The signal this decoder demodulates; signal metadata is forwarded through it
# (see `src/gps/l1ca.jl`). Named for the data component that carries B-CNAV1,
# mirroring how the Galileo E5a decoder is keyed to `GalileoE5aI`.
get_signal_type(::BeiDouB1CConstants) = BeiDouB1C_D

"""
$(TYPEDSIGNATURES)

Reset the BeiDou B1C decoder state after a signal loss or reacquisition.

Clears the in-flight sync state (soft-symbol buffer and SOH) and the
validated data, while preserving the long-lived ephemeris/clock fields in
`raw_data` so a receiver can re-use the satellite after reacquisition without
re-decoding all of subframe 2. Mirrors the GPS L1C-D implementation.

# Arguments

  - `state::GNSSDecoderState{<:BeiDouB1CData}`: Current BeiDou B1C decoder state.

# Returns

  - `GNSSDecoderState{<:BeiDouB1CData}`: Reset decoder state with cleared buffers.

# See Also

  - [`BeiDouB1CDecoderState`](@ref): Create a fresh decoder state
  - [`decode`](@ref): Continue decoding after reset
"""
function reset_decoder_state(state::GNSSDecoderState{<:BeiDouB1CData})
    empty!(state.cache.soft_buffer)
    GNSSDecoderState(
        state;
        raw_data = BeiDouB1CData(state.raw_data; soh = nothing),
        data = BeiDouB1CData(),
        num_bits_after_valid_syncro_sequence = nothing,
        is_shifted_by_180_degrees = false,
    )
end

# ---- Sync (ICD §6.2.2.1 / §7.3) ---------------------------------------------
#
# B-CNAV1 has no preamble; sync runs the two-codeword BCH match directly on
# the soft-symbol deque, with SOH in the role GPS L1C-D gives the TOI: the
# first 72 symbols must equal [PRN₂₁ | SOH₅₁(s)] and the last 72 (the next
# frame's subframe 1) [PRN₂₁ | SOH₅₁((s+1) mod 200)], in a common polarity.
#
# Unlike the GPS TOI construction there is no MSB overlay, and the all-ones
# word is not a codeword of either BCH code (both generator polynomials have
# an odd number of terms), so a codeword's complement is never itself a
# codeword: the polarity is determined unambiguously by the PRN segment —
# even for SOH = 0, whose 51-symbol codeword is all zeros.

"""
    B1CSF1Sync(soh::Int, polarity_flipped::Bool)

Result of a successful two-frame B-CNAV1 subframe-1 match. `soh` is the
Seconds-Of-Hour count (0..199, 18 s units) of the *first* of the two frames.
`polarity_flipped == true` means the receiver is Costas-locked 180° off and
every symbol must be inverted before downstream processing.
"""
struct B1CSF1Sync
    soh::Int
    polarity_flipped::Bool
end

# Decode one 72-symbol subframe-1 window as [PRN₂₁ | SOH₅₁(soh)] under a known
# polarity (`flip`), returning the SOH it carries or `nothing`. Both fields are
# read straight off the deque with the shared `pack_bits_lsb_first`
# (src/gnss.jl), which packs bit 0 = first symbol, the codeword tables' order.
#
# Both BCH fields are decoded to the nearest codeword within their correction
# radius rather than matched exactly, which is what the ICD encodes them for
# (see `B1C_PRN_MAX_ERRORS` / `B1C_SOH_MAX_ERRORS`). Neither step needs a
# tie-break: `2t < d_min` holds for both codes, so at most one codeword can lie
# within the radius, and the search can stop early only on an exact hit.
#
# The order of the two fields is what makes the SOH radius safe. Its 200
# codewords are 24 apart from each other but only 19 from some *complement* of
# another, so a search that had to consider both polarities could not use the
# full radius of 11. It does not have to: the PRN field is decoded first, and
# there its complement distance is 9 against a radius of 3, so a PRN match
# settles the polarity before any SOH symbol is de-flipped.
function _match_sf1_soh(deque, start::Int, prn_codeword::UInt64, flip::Bool)
    prn_word = pack_bits_lsb_first(deque, start, B1C_PRN_CODEWORD_LEN)
    flip && (prn_word = ~prn_word & B1C_PRN_MASK21)
    count_ones(prn_word ⊻ prn_codeword) <= B1C_PRN_MAX_ERRORS || return nothing
    soh_word =
        pack_bits_lsb_first(deque, start + B1C_PRN_CODEWORD_LEN, B1C_SOH_CODEWORD_LEN)
    flip && (soh_word = ~soh_word & B1C_SOH_MASK51)
    best_soh = nothing
    best_distance = B1C_SOH_MAX_ERRORS + 1
    @inbounds for soh = 0:(B1C_SOH_RANGE-1)
        distance = count_ones(soh_word ⊻ B1C_SOH_CODEWORDS[soh+1])
        if distance < best_distance
            best_distance = distance
            best_soh = soh
            distance == 0 && break
        end
    end
    return best_soh
end

"""
    try_sync(state::GNSSDecoderState{<:BeiDouB1CData}) -> Union{Nothing,B1CSF1Sync}

Per-signal sync hook for BeiDou B1C. Hard-slices the leading 72 and trailing
72 soft symbols of the 1872-symbol window and requires, under a common
polarity: the own-PRN BCH(21,6) codeword at both subframe-1 positions, some
SOH `s` in the leading BCH(51,8) segment, and `(s+1) mod 200` in the trailing
one. Returns the [`B1CSF1Sync`](@ref) on a match (carrying the SOH and the
detected polarity flip) or `nothing`.
"""
function try_sync(state::GNSSDecoderState{<:BeiDouB1CData})
    deque = soft_buffer(state)
    prn_codeword = state.cache.prn_codeword
    for flip in (false, true)
        soh = _match_sf1_soh(deque, 1, prn_codeword, flip)
        isnothing(soh) && continue
        soh_next = _match_sf1_soh(deque, B1C_FRAME_LENGTH + 1, prn_codeword, flip)
        if !isnothing(soh_next) && soh_next == (soh + 1) % B1C_SOH_RANGE
            return B1CSF1Sync(soh, flip)
        end
    end
    return nothing
end

"""
    complement_buffer_if_necessary(state::GNSSDecoderState{<:BeiDouB1CData}, sync)

Resolve the 180° polarity ambiguity for BeiDou B1C. `sync::B1CSF1Sync`
already carries the detected polarity from [`try_sync`](@ref); record it on
the state and pass the sync object through unchanged for
`decode_syncro_sequence`.
"""
function complement_buffer_if_necessary(
    state::GNSSDecoderState{<:BeiDouB1CData},
    sync::B1CSF1Sync,
)
    GNSSDecoderState(state; is_shifted_by_180_degrees = sync.polarity_flipped), sync
end

# ---- Frame pipeline ----------------------------------------------------------

"""
    decode_syncro_sequence(state::GNSSDecoderState{<:BeiDouB1CData}, sync::B1CSF1Sync)

Process one locked B-CNAV1 frame. Enforces a monotonic SOH (each frame's SOH
must be `(previous + 1) mod 200`; a break resets the decoder to searching),
then deinterleaves the 1728-symbol payload (36×48 staggered block interleaver,
ICD §6.2.2.4) and LDPC-decodes + CRC-checks subframes 2 and 3, parsing their
fields into `raw_data`.
"""
function decode_syncro_sequence(state::GNSSDecoderState{<:BeiDouB1CData}, sync::B1CSF1Sync)
    prev_soh = state.raw_data.soh
    if !isnothing(prev_soh) && (prev_soh + 1) % B1C_SOH_RANGE != sync.soh
        return reset_decoder_state(state)
    end
    state =
        GNSSDecoderState(state; raw_data = BeiDouB1CData(state.raw_data; soh = sync.soh))

    # Extract the 1728-symbol interleaved SF2+SF3 payload (symbols 73..1800),
    # applying the polarity flip by negating soft symbols up front. All four
    # working buffers live in the cache, so a decoded frame allocates nothing.
    interleaved = copy_soft_window!(
        state.cache.payload,
        soft_buffer(state),
        B1C_SUBFRAME1_LENGTH,
        B1C_PAYLOAD_SYMBOLS,
        state.is_shifted_by_180_degrees,
    )

    # Undo the column-wise readout: `deinterleave!(…, 36, 48)` restores the
    # 36×48 array in row-major order; the staggered row orders then split the
    # rows back into the SF2 and SF3 codewords.
    array_rows = deinterleave!(
        state.cache.array_rows,
        interleaved,
        B1C_INTERLEAVER_ROWS,
        B1C_INTERLEAVER_COLS,
    )
    sf2_symbols = state.cache.sf2_symbols
    for (i, row) in enumerate(B1C_SF2_ROW_ORDER)
        copyto!(
            sf2_symbols,
            (i - 1) * B1C_INTERLEAVER_COLS + 1,
            array_rows,
            (row - 1) * B1C_INTERLEAVER_COLS + 1,
            B1C_INTERLEAVER_COLS,
        )
    end
    sf3_symbols = state.cache.sf3_symbols
    for (i, row) in enumerate(B1C_SF3_ROW_ORDER)
        copyto!(
            sf3_symbols,
            (i - 1) * B1C_INTERLEAVER_COLS + 1,
            array_rows,
            (row - 1) * B1C_INTERLEAVER_COLS + 1,
            B1C_INTERLEAVER_COLS,
        )
    end

    state = decode_b1c_subframe2(state, sf2_symbols)
    state = decode_b1c_subframe3(state, sf3_symbols)
    return state
end

# ---- Subframe 2 bit-field extraction (ICD Figure 6-6, 6-12 .. 6-14) --------
#
# `word` is the 600-bit subframe-2 info block packed MSB-first into a
# `UInt600` (bit 1 = MSB). Fields are read by 1-based start bit and length
# through the shared `get_bits` / `get_twos_complement_num` helpers. Layout:
# WN(13) HOW(8) IODC(10) IODE(8) EphemerisI(203) EphemerisII(222)
# Clock(69) T_GD_B2ap(12) ISC_B1Cd(12) T_GD_B1Cp(12) Rev(7) CRC(24).

function decode_b1c_subframe2(state::GNSSDecoderState{<:BeiDouB1CData}, sf2_symbols)
    word = ldpc_decode_word(state.cache.sf2_ldpc, sf2_symbols, UInt600)
    isnothing(word) && return state  # silently drop on CRC failure
    word_length = B1C_SF2_INFO_BITS

    PI = state.constants.PI

    WN = Int(get_bits(word, word_length, 1, 13))
    HOW = Int(get_bits(word, word_length, 14, 8))
    IODC = Int(get_bits(word, word_length, 22, 10))
    IODE = Int(get_bits(word, word_length, 32, 8))
    # Ephemeris I (Figure 6-12, Table 7-8), bits 40-242.
    t_0e = Int(get_bits(word, word_length, 40, 11)) * 300
    sat_type = Int(get_bits(word, word_length, 51, 2))
    ΔA = get_twos_complement_num(word, word_length, 53, 26) * 2.0^-9
    A_dot = get_twos_complement_num(word, word_length, 79, 25) * 2.0^-21
    Δn_0 = get_twos_complement_num(word, word_length, 104, 17) * 2.0^-44 * PI
    Δn_0_dot = get_twos_complement_num(word, word_length, 121, 23) * 2.0^-57 * PI
    M_0 = get_twos_complement_num(word, word_length, 144, 33) * 2.0^-32 * PI
    # 33-bit field: convert through Int64 explicitly (a platform `Int` would
    # overflow on 32-bit), matching the sibling B2a/B2b parsers.
    e = Int64(get_bits(word, word_length, 177, 33)) * 2.0^-34
    ω = get_twos_complement_num(word, word_length, 210, 33) * 2.0^-32 * PI
    # Ephemeris II (Figure 6-13, Table 7-8), bits 243-464.
    Ω_0 = get_twos_complement_num(word, word_length, 243, 33) * 2.0^-32 * PI
    i_0 = get_twos_complement_num(word, word_length, 276, 33) * 2.0^-32 * PI
    Ω_dot = get_twos_complement_num(word, word_length, 309, 19) * 2.0^-44 * PI
    i_dot = get_twos_complement_num(word, word_length, 328, 15) * 2.0^-44 * PI
    C_is = get_twos_complement_num(word, word_length, 343, 16) * 2.0^-30
    C_ic = get_twos_complement_num(word, word_length, 359, 16) * 2.0^-30
    C_rs = get_twos_complement_num(word, word_length, 375, 24) * 2.0^-8
    C_rc = get_twos_complement_num(word, word_length, 399, 24) * 2.0^-8
    C_us = get_twos_complement_num(word, word_length, 423, 21) * 2.0^-30
    C_uc = get_twos_complement_num(word, word_length, 444, 21) * 2.0^-30
    # Clock correction (Figure 6-14, Table 7-5), bits 465-533.
    (; t_0c, a_f0, a_f1, a_f2) = beidou_clock_block(word, word_length, 465)
    # Group delay differentials (Table 7-6), bits 534-569.
    T_GD_B2ap = get_twos_complement_num(word, word_length, 534, 12) * 2.0^-34
    ISC_B1Cd = get_twos_complement_num(word, word_length, 546, 12) * 2.0^-34
    T_GD_B1Cp = get_twos_complement_num(word, word_length, 558, 12) * 2.0^-34

    raw = BeiDouB1CData(
        state.raw_data;
        WN,
        HOW,
        IODC,
        IODE,
        t_0e,
        sat_type,
        ΔA,
        A_dot,
        Δn_0,
        Δn_0_dot,
        M_0,
        e,
        ω,
        Ω_0,
        i_0,
        Ω_dot,
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
        T_GD_B2ap,
        ISC_B1Cd,
        T_GD_B1Cp,
    )
    GNSSDecoderState(state; raw_data = raw)
end

# ---- Subframe 3 page parsing (ICD Figures 6-7 .. 6-11) ----------------------
#
# Every SF3 page is a 264-bit info block: bits 1-6 the PageID, then (for the
# defined page types 1-4) a common header HS(2) DIF(1) SIF(1) AIF(1)
# SISMAI(4) at bits 7-15, page-specific fields, and a trailing CRC-24Q. After
# the CRC passes the 264 bits are packed MSB-first into a `UInt288`
# (`get_bits(word, 264, …)` addresses the right-aligned logical bits); we
# dispatch on the PageID and merge parsed fields into `raw_data` immutably.

function decode_b1c_subframe3(state::GNSSDecoderState{<:BeiDouB1CData}, sf3_symbols)
    word = ldpc_decode_word(state.cache.sf3_ldpc, sf3_symbols, UInt288)
    isnothing(word) && return state  # silently drop on CRC failure

    raw = BeiDouB1CData(
        state.raw_data;
        num_sf3_pages_received = state.raw_data.num_sf3_pages_received + 1,
    )

    word_length = B1C_SF3_INFO_BITS
    page = Int(get_bits(word, word_length, 1, 6))
    # PageID 0 is invalid and > 4 reserved (Table 7-1); their layouts are
    # undefined, so only the header of the defined page types is parsed.
    if 1 <= page <= 4
        raw = BeiDouB1CData(
            raw;
            HS = Int(get_bits(word, word_length, 7, 2)),
            DIF = get_bit(word, word_length, 9),
            SIF = get_bit(word, word_length, 10),
            AIF = get_bit(word, word_length, 11),
            SISMAI = Int(get_bits(word, word_length, 12, 4)),
        )
    end
    raw = if page == 1
        parse_b1c_sf3_page1(raw, word)
    elseif page == 2
        parse_b1c_sf3_page2(raw, word, state.constants.PI)
    elseif page == 3
        parse_b1c_sf3_page3(raw, word)
    elseif page == 4
        parse_b1c_sf3_page4(raw, word, state.constants.PI)
    else
        raw  # invalid/reserved page: counted, ignored
    end

    GNSSDecoderState(state; raw_data = raw)
end

# SISAI_oc data block (Figure 6-15), 22 bits at 1-based `start`; the shared
# parser is in `beidou/beidou.jl` because B2a and B2b carry the same block.
_parse_b1c_sisai_oc(raw::BeiDouB1CData, word::UInt288, start::Int) =
    BeiDouB1CData(raw; beidou_sisai_oc_block(word, B1C_SF3_INFO_BITS, start)...)

"""
Subframe 3, page type 1 — SISAI + BDGIM iono + BDT-UTC (ICD Figures 6-8, 6-16, 6-17).
"""
function parse_b1c_sf3_page1(raw::BeiDouB1CData, word::UInt288)
    word_length = B1C_SF3_INFO_BITS
    raw = BeiDouB1CData(raw; SISAI_oe = Int(get_bits(word, word_length, 16, 5)))
    raw = _parse_b1c_sisai_oc(raw, word, 21)
    BeiDouB1CData(
        raw;
        # BDGIM ionospheric coefficients (Figure 6-16, Table 7-10), bits 43-116:
        beidou_bdgim_block(word, word_length, 43)...,
        # BDT-UTC time offset (Figure 6-17, Table 7-20).
        # BDT-UTC (Figure 6-17, Table 7-20), bits 117-213:
        beidou_bdt_utc_block(word, word_length, 117)...,
    )
end

"""
Subframe 3, page type 2 — SISAI + four reduced almanacs (ICD Figures 6-9, 6-18).
"""
function parse_b1c_sf3_page2(raw::BeiDouB1CData, word::UInt288, PI::Float64)
    word_length = B1C_SF3_INFO_BITS
    raw = _parse_b1c_sisai_oc(raw, word, 16)
    WN_a = Int(get_bits(word, word_length, 38, 13))
    t_0a = Int(get_bits(word, word_length, 51, 8)) * 2^12
    almanacs = raw.reduced_almanacs
    # Four 38-bit packets at bits 59, 97, 135, 173 (Figure 6-9).
    for start in (59, 97, 135, 173)
        packet = beidou_reduced_almanac(word, word_length, start, WN_a, t_0a, PI)
        isnothing(packet) && continue
        almanacs = _merge_keyed(almanacs, packet.PRN_a, packet)
    end
    BeiDouB1CData(raw; reduced_almanacs = almanacs)
end

"""
Subframe 3, page type 3 — SISAI + EOP + BGTO (ICD Figures 6-10, 6-19, 6-20).
"""
function parse_b1c_sf3_page3(raw::BeiDouB1CData, word::UInt288)
    word_length = B1C_SF3_INFO_BITS
    raw = BeiDouB1CData(
        raw;
        SISAI_oe = Int(get_bits(word, word_length, 16, 5)),
        # EOP (Figure 6-19, Table 7-18), bits 21-158:
        beidou_eop_block(word, word_length, 21)...,
    )
    # BGTO (Figure 6-20, Table 7-21), bits 159-226. GNSS ID 0 means "not
    # available" (§7.13.1) — nothing is stored then.
    GNSS_ID = Int(get_bits(word, word_length, 159, 3))
    GNSS_ID == 0 && return raw
    bgto = BeiDouB1CBGTO(; beidou_bgto_block(word, word_length, 159)...)
    BeiDouB1CData(raw; bgtos = _merge_keyed(raw.bgtos, GNSS_ID, bgto))
end

"""
Subframe 3, page type 4 — SISAI + one midi almanac (ICD Figures 6-11, 6-21, Table 7-13).
"""
function parse_b1c_sf3_page4(raw::BeiDouB1CData, word::UInt288, PI::Float64)
    word_length = B1C_SF3_INFO_BITS
    raw = _parse_b1c_sisai_oc(raw, word, 16)
    # Midi almanac (Figure 6-21), bits 38-193.
    alm = beidou_midi_almanac(word, word_length, 38, PI)
    isnothing(alm) && return raw  # empty almanac slot
    BeiDouB1CData(raw; midi_almanacs = _merge_keyed(raw.midi_almanacs, alm.PRN_a, alm))
end

"""
    validate_data(state::GNSSDecoderState{<:BeiDouB1CData})

Promote `raw_data` to `data` once the positioning set is complete (SOH,
subframe 2, and the subframe-3 health status) and the ephemeris and clock sets
are a matched pair. The per-frame BCH re-check and monotonic-SOH enforcement
are performed inline by [`decode_syncro_sequence`](@ref); this hook publishes
validated data and arms the streaming counter.

Sharing one CRC-protected block does not make IODE and IODC a matched pair —
see [`is_decoding_completed_for_positioning`](@ref), which carries that gate.

The counter is armed at the epoch the SOH stamps — the leading edge of the
*current* subframe 1 (ICD §7.3) — which at promotion time is
`syncro_sequence_length + preamble_length` (1872 symbols = 18.72 s) behind the
newest buffered symbol, as in `b2a.jl` and `b2b.jl`. GPS L1C-D's TOI stamps
the *next* frame instead, so `gps/l1c_d.jl` anchors to `preamble_length`
alone; using that anchor here would report the time 18 s early.
"""
function validate_data(state::GNSSDecoderState{<:BeiDouB1CData})
    if is_decoding_completed_for_positioning(state.raw_data)
        return GNSSDecoderState(
            state;
            data = state.raw_data,
            num_bits_after_valid_syncro_sequence = state.constants.syncro_sequence_length +
                                                   state.constants.preamble_length,
        )
    end
    return state
end

"""
$(TYPEDSIGNATURES)

Check if the BeiDou B1C satellite is healthy and usable for positioning.

Examines the 2-bit satellite health status (HS) broadcast in every subframe-3
page (BDS-SIS-ICD-B1C-1.0 §7.14, Table 7-22): a satellite is healthy iff
`HS == 0` (the satellite provides services). The B1C integrity status flags
(`DIF`/`SIF`/`AIF`, §7.15) are reported separately on
[`BeiDouB1CData`](@ref) and are deliberately not folded in here — they flag
message/signal integrity for precision users, not the satellite's service
state.

!!! warning

    Requires a subframe-3 page to have been decoded and the positioning set
    validated; returns `false` until then.

# Arguments

  - `state::GNSSDecoderState{<:BeiDouB1CData}`: BeiDou B1C decoder state.

# Returns

  - `Bool`: `true` iff the health status word indicates a healthy satellite.
"""
function is_sat_healthy(state::GNSSDecoderState{<:BeiDouB1CData})
    # `HS` is stored as Int64; compare against a typed zero so the egal check
    # cannot fail on a platform where the literal `0` is Int32.
    state.data.HS === Int64(0)
end
