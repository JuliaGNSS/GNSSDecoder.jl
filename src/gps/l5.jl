# GPS L5I (CNAV) decoder — IS-GPS-705J.
#
# The L5 in-phase data channel carries the CNAV message as 300-bit messages at
# 50 bps, FEC-encoded with the rate-1/2, constraint-length-7 (K=7)
# non-systematic convolutional (NSC) code (G1 = 0o171, G2 = 0o133,
# IS-GPS-705J §3.3.3.1.1) to 600 channel symbols at 100 sps. Every message
# starts with the 8-bit preamble `10001011` and ends with a 24-bit CRC-24Q
# over the whole message (IS-GPS-705J §20.3.5).
#
# Unlike Galileo E1B pages, the CNAV FEC is convolved *continuously* across
# message boundaries — there are no tail bits and no encoder reset — so the
# preamble is only visible after FEC decoding and AFF3CT's tail-terminated
# `ConvViterbiDecoder` does not apply. Message sync therefore follows the
# generic streaming framework in `src/gnss.jl` with a window-decoding
# `try_sync`: the soft-symbol `CircularDeque` is sized to one full message
# (600 symbols) plus the 16-symbol encoding of the *next* message's start
# (616 total); each sync attempt Viterbi-decodes the whole window on the raw
# `Float32` LLRs and requires the preamble at *both ends* of the decoded
# window (in either polarity, mirroring `find_preamble`) plus a clean
# CRC-24Q before locking. Symbol-pair phase ambiguity needs no special
# handling — a misaligned window simply fails to sync and slides one symbol.

"""
One CNAV message: 300 bits = 600 channel symbols (rate 1/2) at 100 sps.
"""
const L5I_MESSAGE_BITS = 300
const L5I_MESSAGE_SYMBOLS = 2 * L5I_MESSAGE_BITS
"""
Trailing window: the next message's first 8 bits (16 symbols) confirm sync.
"""
const L5I_PREAMBLE_BITS = 8
const L5I_PREAMBLE_SYMBOLS = 2 * L5I_PREAMBLE_BITS
const L5I_WINDOW_SYMBOLS = L5I_MESSAGE_SYMBOLS + L5I_PREAMBLE_SYMBOLS  # 616
const L5I_WINDOW_BITS = L5I_MESSAGE_BITS + L5I_PREAMBLE_BITS  # 308

"""
CNAV preamble `10001011` (IS-GPS-705J §20.3.3).
"""
const L5I_PREAMBLE = 0b10001011

# Semi-major axis reference (IS-GPS-705J Table 20-I, meters).
const L5I_A_REF = 26_559_710.0
# Rate-of-right-ascension reference (IS-GPS-705J Table 20-I, semi-circles/sec).
const L5I_OMEGA_DOT_REF = -2.6e-9

"""
    GPSL5IConstants

WGS 84 constants and CNAV message structure parameters for GPS L5I decoding.

The message is modelled through the generic streaming framework:
`syncro_sequence_length` is the 600-symbol message that is drained once
decoded, and `preamble_length` is the 16-symbol encoding of the *next*
message's preamble retained at the tail of the sync window.

# Fields

$(TYPEDFIELDS)

# Reference

IS-GPS-705J, Sections 20.3.3 and 20.3.4.3.
"""
Base.@kwdef struct GPSL5IConstants <: AbstractGNSSConstants
    """
    Message length drained after each decoded message (600 symbols)
    """
    syncro_sequence_length::Int = L5I_MESSAGE_SYMBOLS
    """
    Trailing next-message preamble segment retained for sync (16 symbols)
    """
    preamble_length::Int = L5I_PREAMBLE_SYMBOLS
    """
    Mathematical constant π (IS-GPS-705J Table 20-II)
    """
    PI::Float64 = 3.1415926535898
    """
    WGS 84 Earth rotation rate (rad/s)
    """
    Ω_dot_e::Float64 = 7.2921151467e-5
    """
    Speed of light (m/s)
    """
    c::Float64 = 2.99792458e8
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
    GPSL5IReducedAlmanac

One satellite's reduced-almanac packet from CNAV message types 12 or 31
(IS-GPS-705J Figure 20-16, Table 20-VI).

The reduced almanac gives a very coarse ephemeris for satellite selection.
Message type 12 carries seven 31-bit packets, message type 31 four; each
packet is complete in itself, so `GPSL5IData.reduced_almanacs` entries are
inserted whole, keyed by `PRN_a` (mirrors [`GPSL1C_DReducedAlmanac`](@ref)).

Reference values to apply (Table 20-VI footnotes): `e = 0`,
`δi = +0.0056 semi-circles` (so `i = 55°`),
`Ω̇ = -2.6e-9 semi-circles/s`, `A = A_ref + δA` with `A_ref = 26 559 710 m`,
`Φ₀ = M₀ + ω`. Semi-circle fields are converted to radians on decode.

# Fields

  - `PRN_a::Int`: Almanac satellite PRN (1-63; 0 marks an empty packet).
  - `WN_a::Int`: Almanac reference week number (mod 8192).
  - `t_oa::Int`: Almanac reference time of week (seconds).
  - `δA::Float64`: Semi-major-axis delta from `A_ref` (meters).
  - `Ω_0::Float64`: Longitude of ascending node at weekly epoch (rad).
  - `Φ_0::Float64`: Argument of latitude at reference time, `M₀+ω` (rad).
  - `l1_health::Bool`, `l2_health::Bool`, `l5_health::Bool`: per-band health
    (false = OK, true = some/all signals bad).

# Reference

IS-GPS-705J, Figures 20-4 / 20-11 / 20-16, Table 20-VI.
"""
Base.@kwdef struct GPSL5IReducedAlmanac
    PRN_a::Int
    WN_a::Int
    t_oa::Int
    δA::Float64
    Ω_0::Float64
    Φ_0::Float64
    l1_health::Bool
    l2_health::Bool
    l5_health::Bool
end

"""
    GPSL5IMidiAlmanac

One satellite's Midi almanac from CNAV message type 37 (IS-GPS-705J
Figure 20-10, Table 20-V).

The Midi almanac is a medium-precision single-SV almanac, complete in a
single message, so `GPSL5IData.midi_almanacs` entries are inserted whole,
keyed by `PRN_a` (mirrors [`GPSL1C_DMidiAlmanac`](@ref)). Inclination is `δi`
relative to `i₀ = 0.30 semi-circles` (54°); semi-circle fields are converted
to radians on decode.

# Fields

  - `PRN_a::Int`: Almanac satellite PRN.
  - `WN_a::Int`: Almanac reference week number (mod 8192).
  - `t_oa::Int`: Almanac reference time of week (seconds).
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

IS-GPS-705J, Figure 20-10, Table 20-V.
"""
Base.@kwdef struct GPSL5IMidiAlmanac
    PRN_a::Int
    WN_a::Int
    t_oa::Int
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
    GPSL5IClockDifferentialCorrection

One satellite's clock differential-correction (CDC) packet from CNAV message
types 13 or 34, keyed by `PRN_a` (IS-GPS-705J Figure 20-17, Table 20-X).

`dc_data_type` selects the data the corrections apply to: `false` ⇒ CNAV
(message types 30-37), `true` ⇒ legacy NAV. A `PRN ID` of all-ones marks an
empty packet (not stored).

# Fields

  - `PRN_a::Int`: PRN the corrections apply to.
  - `t_op_D::Int`: DC data predict time of week (seconds).
  - `t_OD::Int`: Time of DC data (seconds).
  - `dc_data_type::Bool`: false ⇒ CNAV, true ⇒ legacy NAV.
  - `δa_f0::Float64`, `δa_f1::Float64`: Clock bias / drift corrections (s, s/s).
  - `UDRA_index::Int`: UDRA index (signed).

# Reference

IS-GPS-705J, Figures 20-7 / 20-12 / 20-17, Table 20-X.
"""
Base.@kwdef struct GPSL5IClockDifferentialCorrection
    PRN_a::Int
    t_op_D::Int
    t_OD::Int
    dc_data_type::Bool
    δa_f0::Float64
    δa_f1::Float64
    UDRA_index::Int
end

"""
    GPSL5IEphemerisDifferentialCorrection

One satellite's ephemeris differential-correction (EDC) packet from CNAV
message types 14 or 34, keyed by `PRN_a` (IS-GPS-705J Figure 20-17,
Table 20-X). Semi-circle fields are converted to radians on decode.

# Fields

  - `PRN_a::Int`: PRN the corrections apply to.
  - `t_op_D::Int`: DC data predict time of week (seconds).
  - `t_OD::Int`: Time of DC data (seconds).
  - `dc_data_type::Bool`: false ⇒ CNAV, true ⇒ legacy NAV.
  - `Δα::Float64`, `Δβ::Float64`: Ephemeris α/β corrections (dimensionless).
  - `Δγ::Float64`: Ephemeris γ correction (rad).
  - `Δi::Float64`, `ΔΩ::Float64`: Inclination / right-ascension corrections (rad).
  - `ΔA::Float64`: Semi-major-axis correction (meters).
  - `UDRA_dot_index::Int`: Rate-of-UDRA index (signed).

# Reference

IS-GPS-705J, Figures 20-7 / 20-13 / 20-17, Table 20-X.
"""
Base.@kwdef struct GPSL5IEphemerisDifferentialCorrection
    PRN_a::Int
    t_op_D::Int
    t_OD::Int
    dc_data_type::Bool
    Δα::Float64
    Δβ::Float64
    Δγ::Float64
    Δi::Float64
    ΔΩ::Float64
    ΔA::Float64
    UDRA_dot_index::Int
end

"""
    GPSL5IIntegritySupportMessage

Integrity Support Message from CNAV message type 40 (ARAIM), complete in a
single message (IS-GPS-705J Figure 20-14a, Table 20-XIa).

# Fields

  - `GNSS_ID::Int`: GNSS identifier the ISM applies to.
  - `WN_ISM::Int`, `TOW_ISM::Int`: ISM reference week / time-of-week counts.
  - `t_correl::Int`, `b_nom::Int`, `γ_nom::Int`, `R_sat::Int`, `P_const::Int`,
    `MFD::Int`, `service_level::Int`: encoded ARAIM parameter indices.
  - `mask::UInt64`: 63-bit SV mask (MSB = PRN 1).

# Reference

IS-GPS-705J, Figure 20-14a, Table 20-XIa.
"""
Base.@kwdef struct GPSL5IIntegritySupportMessage
    GNSS_ID::Int
    WN_ISM::Int
    TOW_ISM::Int
    t_correl::Int
    b_nom::Int
    γ_nom::Int
    R_sat::Int
    P_const::Int
    MFD::Int
    service_level::Int
    mask::UInt64
end

"""
    GPSL5IData

Decoded GPS L5I CNAV navigation message data.

Holds the parameters decoded from CNAV message types 10, 11, 12, 13, 14, 15,
30-37, and 40 (IS-GPS-705J §20.3.3). The decoder fills fields incrementally
as the corresponding message types are received. Field-naming follows
[`GPSL1C_DData`](@ref) (CNAV-2 broadcasts nearly the same parameter set):
semi-circle quantities are converted to radians on decode (multiplied by π),
all `Union{Nothing,…}` until first decoded.

# Header (every message)

  - `last_message_id::Int`: Most recently decoded message type (0 until then).
  - `TOW::Int64`: SV time in seconds at the start of the *next* 6-second
    message (message TOW count × 6).
  - `alert_flag::Bool`: Raised when the signal URA may be worse than indicated.

# Health / accuracy (message types 10, 30-37)

  - `l1_health::Bool`, `l2_health::Bool`, `l5_health::Bool`: per-band signal
    health (false = OK).
  - `ura_ed_index::Int64`: Ephemeris URA index (signed).
  - `ura_ned0_index::Int64`, `ura_ned1_index::Int64`, `ura_ned2_index::Int64`:
    Clock URA indices.

# Ephemeris (message types 10 + 11, Table 20-I)

  - `WN::Int64`: Transmission week number, modulo-8192.
  - `t_op::Int64`: Data predict time of week (seconds).
  - `t_0e::Int64`: Ephemeris data reference time of week (seconds).
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
  - `i_0_dot::Float64`: Rate of inclination angle (rad/s).
  - `C_is::Float64`, `C_ic::Float64`: Sine/cosine inclination harmonic corrections (rad).
  - `C_rs::Float64`, `C_rc::Float64`: Sine/cosine orbit-radius harmonic corrections (m).
  - `C_us::Float64`, `C_uc::Float64`: Sine/cosine argument-of-latitude harmonic corrections (rad).
  - `integrity_status_flag::Bool`, `l2c_phasing::Bool`: message type 10 flags.

# Clock (message types 30-37, Table 20-III)

  - `t_0c::Int64`: Clock data reference time of week (seconds).
  - `a_f0::Float64`, `a_f1::Float64`, `a_f2::Float64`: Clock bias / drift / drift-rate.

# Group delay / ISC + ionosphere (message type 30, Tables 20-III / 20-IV)

  - `T_GD::Float64`: L1/L2 P(Y) inter-signal correction (seconds).
  - `ISC_L1CA,ISC_L2C,ISC_L5I5,ISC_L5Q5::Float64`: inter-signal corrections (s).
  - `α0,α1,α2,α3,β0,β1,β2,β3::Float64`: Klobuchar ionospheric coefficients.
  - `WN_op::Int64`: Data predict week number (mod 256).

# EOP (message type 32, Table 20-VII)

  - `t_EOP::Int64`: EOP reference time of week (s).
  - `PM_X,PM_X_dot,PM_Y,PM_Y_dot::Float64`: polar-motion values/rates (arcsec, arcsec/day).
  - `ΔUT_GPS,ΔUT_GPS_dot::Float64`: UT1-GPS difference (s) and rate (s/day).

# UTC (message type 33, Table 20-IX)

  - `A0_UTC,A1_UTC,A2_UTC::Float64`: UTC polynomial (s, s/s, s/s²).
  - `Δt_LS,Δt_LSF::Int64`: current/past and future leap-second counts (s).
  - `t_ot::Int64`: UTC reference time of week (s).
  - `WN_ot,WN_LSF::Int64`: UTC and leap-second reference week numbers.
  - `DN::Int64`: leap-second reference day number (1-7).

# GGTO (message type 35, Table 20-XI)

  - `A0_GGTO,A1_GGTO,A2_GGTO::Float64`: GPS/GNSS time-offset polynomial.
  - `t_GGTO::Int64`, `WN_GGTO::Int64`: GGTO reference time/week.
  - `GNSS_ID::Int64`: 0 none, 1 Galileo, 2 GLONASS, 3 BeiDou, 4-7 reserved.

# Almanacs / corrections / text — keyed dictionaries (`nothing` until first decoded)

  - `reduced_almanacs::Dictionary{Int,GPSL5IReducedAlmanac}` (message types 12, 31).
  - `midi_almanacs::Dictionary{Int,GPSL5IMidiAlmanac}` (message type 37).
  - `clock_corrections::Dictionary{Int,GPSL5IClockDifferentialCorrection}`
    (message types 13, 34).
  - `ephemeris_corrections::Dictionary{Int,GPSL5IEphemerisDifferentialCorrection}`
    (message types 14, 34).
  - `text_mt15::String`, `text_page_mt15::Int64`: message type 15 text page
    (29 ASCII characters, control chars stripped).
  - `text_mt36::String`, `text_page_mt36::Int64`: message type 36 text page
    (18 ASCII characters).
  - `ism::GPSL5IIntegritySupportMessage`: message type 40 Integrity Support Message.

# Reference

IS-GPS-705J, Figures 20-1 through 20-17 and Tables 20-I through 20-XIa.
"""
Base.@kwdef struct GPSL5IData <: AbstractGNSSData
    last_message_id::Int = 0
    TOW::Union{Nothing,Int64} = nothing
    alert_flag::Union{Nothing,Bool} = nothing

    l1_health::Union{Nothing,Bool} = nothing
    l2_health::Union{Nothing,Bool} = nothing
    l5_health::Union{Nothing,Bool} = nothing
    ura_ed_index::Union{Nothing,Int64} = nothing
    ura_ned0_index::Union{Nothing,Int64} = nothing
    ura_ned1_index::Union{Nothing,Int64} = nothing
    ura_ned2_index::Union{Nothing,Int64} = nothing

    WN::Union{Nothing,Int64} = nothing
    t_op::Union{Nothing,Int64} = nothing
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
    i_0_dot::Union{Nothing,Float64} = nothing
    C_is::Union{Nothing,Float64} = nothing
    C_ic::Union{Nothing,Float64} = nothing
    C_rs::Union{Nothing,Float64} = nothing
    C_rc::Union{Nothing,Float64} = nothing
    C_us::Union{Nothing,Float64} = nothing
    C_uc::Union{Nothing,Float64} = nothing
    integrity_status_flag::Union{Nothing,Bool} = nothing
    l2c_phasing::Union{Nothing,Bool} = nothing

    t_0c::Union{Nothing,Int64} = nothing
    a_f0::Union{Nothing,Float64} = nothing
    a_f1::Union{Nothing,Float64} = nothing
    a_f2::Union{Nothing,Float64} = nothing

    T_GD::Union{Nothing,Float64} = nothing
    ISC_L1CA::Union{Nothing,Float64} = nothing
    ISC_L2C::Union{Nothing,Float64} = nothing
    ISC_L5I5::Union{Nothing,Float64} = nothing
    ISC_L5Q5::Union{Nothing,Float64} = nothing
    α0::Union{Nothing,Float64} = nothing
    α1::Union{Nothing,Float64} = nothing
    α2::Union{Nothing,Float64} = nothing
    α3::Union{Nothing,Float64} = nothing
    β0::Union{Nothing,Float64} = nothing
    β1::Union{Nothing,Float64} = nothing
    β2::Union{Nothing,Float64} = nothing
    β3::Union{Nothing,Float64} = nothing
    WN_op::Union{Nothing,Int64} = nothing

    t_EOP::Union{Nothing,Int64} = nothing
    PM_X::Union{Nothing,Float64} = nothing
    PM_X_dot::Union{Nothing,Float64} = nothing
    PM_Y::Union{Nothing,Float64} = nothing
    PM_Y_dot::Union{Nothing,Float64} = nothing
    ΔUT_GPS::Union{Nothing,Float64} = nothing
    ΔUT_GPS_dot::Union{Nothing,Float64} = nothing

    A0_UTC::Union{Nothing,Float64} = nothing
    A1_UTC::Union{Nothing,Float64} = nothing
    A2_UTC::Union{Nothing,Float64} = nothing
    Δt_LS::Union{Nothing,Int64} = nothing
    t_ot::Union{Nothing,Int64} = nothing
    WN_ot::Union{Nothing,Int64} = nothing
    WN_LSF::Union{Nothing,Int64} = nothing
    DN::Union{Nothing,Int64} = nothing
    Δt_LSF::Union{Nothing,Int64} = nothing

    A0_GGTO::Union{Nothing,Float64} = nothing
    A1_GGTO::Union{Nothing,Float64} = nothing
    A2_GGTO::Union{Nothing,Float64} = nothing
    t_GGTO::Union{Nothing,Int64} = nothing
    WN_GGTO::Union{Nothing,Int64} = nothing
    GNSS_ID::Union{Nothing,Int64} = nothing

    reduced_almanacs::Union{Nothing,Dictionary{Int,GPSL5IReducedAlmanac}} = nothing
    midi_almanacs::Union{Nothing,Dictionary{Int,GPSL5IMidiAlmanac}} = nothing
    clock_corrections::Union{Nothing,Dictionary{Int,GPSL5IClockDifferentialCorrection}} =
        nothing
    ephemeris_corrections::Union{
        Nothing,
        Dictionary{Int,GPSL5IEphemerisDifferentialCorrection},
    } = nothing

    text_mt15::Union{Nothing,String} = nothing
    text_page_mt15::Union{Nothing,Int64} = nothing
    text_mt36::Union{Nothing,String} = nothing
    text_page_mt36::Union{Nothing,Int64} = nothing

    ism::Union{Nothing,GPSL5IIntegritySupportMessage} = nothing
end

function GPSL5IData(
    data::GPSL5IData;
    last_message_id = data.last_message_id,
    TOW = data.TOW,
    alert_flag = data.alert_flag,
    l1_health = data.l1_health,
    l2_health = data.l2_health,
    l5_health = data.l5_health,
    ura_ed_index = data.ura_ed_index,
    ura_ned0_index = data.ura_ned0_index,
    ura_ned1_index = data.ura_ned1_index,
    ura_ned2_index = data.ura_ned2_index,
    WN = data.WN,
    t_op = data.t_op,
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
    i_0_dot = data.i_0_dot,
    C_is = data.C_is,
    C_ic = data.C_ic,
    C_rs = data.C_rs,
    C_rc = data.C_rc,
    C_us = data.C_us,
    C_uc = data.C_uc,
    integrity_status_flag = data.integrity_status_flag,
    l2c_phasing = data.l2c_phasing,
    t_0c = data.t_0c,
    a_f0 = data.a_f0,
    a_f1 = data.a_f1,
    a_f2 = data.a_f2,
    T_GD = data.T_GD,
    ISC_L1CA = data.ISC_L1CA,
    ISC_L2C = data.ISC_L2C,
    ISC_L5I5 = data.ISC_L5I5,
    ISC_L5Q5 = data.ISC_L5Q5,
    α0 = data.α0,
    α1 = data.α1,
    α2 = data.α2,
    α3 = data.α3,
    β0 = data.β0,
    β1 = data.β1,
    β2 = data.β2,
    β3 = data.β3,
    WN_op = data.WN_op,
    t_EOP = data.t_EOP,
    PM_X = data.PM_X,
    PM_X_dot = data.PM_X_dot,
    PM_Y = data.PM_Y,
    PM_Y_dot = data.PM_Y_dot,
    ΔUT_GPS = data.ΔUT_GPS,
    ΔUT_GPS_dot = data.ΔUT_GPS_dot,
    A0_UTC = data.A0_UTC,
    A1_UTC = data.A1_UTC,
    A2_UTC = data.A2_UTC,
    Δt_LS = data.Δt_LS,
    t_ot = data.t_ot,
    WN_ot = data.WN_ot,
    WN_LSF = data.WN_LSF,
    DN = data.DN,
    Δt_LSF = data.Δt_LSF,
    A0_GGTO = data.A0_GGTO,
    A1_GGTO = data.A1_GGTO,
    A2_GGTO = data.A2_GGTO,
    t_GGTO = data.t_GGTO,
    WN_GGTO = data.WN_GGTO,
    GNSS_ID = data.GNSS_ID,
    reduced_almanacs = data.reduced_almanacs,
    midi_almanacs = data.midi_almanacs,
    clock_corrections = data.clock_corrections,
    ephemeris_corrections = data.ephemeris_corrections,
    text_mt15 = data.text_mt15,
    text_page_mt15 = data.text_page_mt15,
    text_mt36 = data.text_mt36,
    text_page_mt36 = data.text_page_mt36,
    ism = data.ism,
)
    GPSL5IData(
        last_message_id,
        TOW,
        alert_flag,
        l1_health,
        l2_health,
        l5_health,
        ura_ed_index,
        ura_ned0_index,
        ura_ned1_index,
        ura_ned2_index,
        WN,
        t_op,
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
        i_0_dot,
        C_is,
        C_ic,
        C_rs,
        C_rc,
        C_us,
        C_uc,
        integrity_status_flag,
        l2c_phasing,
        t_0c,
        a_f0,
        a_f1,
        a_f2,
        T_GD,
        ISC_L1CA,
        ISC_L2C,
        ISC_L5I5,
        ISC_L5Q5,
        α0,
        α1,
        α2,
        α3,
        β0,
        β1,
        β2,
        β3,
        WN_op,
        t_EOP,
        PM_X,
        PM_X_dot,
        PM_Y,
        PM_Y_dot,
        ΔUT_GPS,
        ΔUT_GPS_dot,
        A0_UTC,
        A1_UTC,
        A2_UTC,
        Δt_LS,
        t_ot,
        WN_ot,
        WN_LSF,
        DN,
        Δt_LSF,
        A0_GGTO,
        A1_GGTO,
        A2_GGTO,
        t_GGTO,
        WN_GGTO,
        GNSS_ID,
        reduced_almanacs,
        midi_almanacs,
        clock_corrections,
        ephemeris_corrections,
        text_mt15,
        text_page_mt15,
        text_mt36,
        text_page_mt36,
        ism,
    )
end

# The default struct `==` falls back to `===` (reference equality), which fails
# for the mutable `Dictionary` fields even when their contents match. Compare
# field-by-field (mirrors `GPSL1C_DData`).
function Base.:(==)(a::GPSL5IData, b::GPSL5IData)
    for f in fieldnames(GPSL5IData)
        getfield(a, f) == getfield(b, f) || return false
    end
    return true
end

"""
$(TYPEDEF)

Per-decoder cache for the GPS L5I signal.

Holds the soft-symbol `CircularDeque{Float32}` (capacity = 616 = 600 message
symbols + 16 next-message preamble symbols). The K=7 NSC FEC is undone on the
raw `Float32` LLRs by the window-decoding Viterbi in `try_sync`; no other
mutable state is needed.

# Fields

$(TYPEDFIELDS)
"""
struct GPSL5ICache <: AbstractGNSSCache
    """
    Soft-symbol buffer (616 = 600 message + 16 next-message preamble)
    """
    soft_buffer::CircularDeque{Float32}
end

GPSL5ICache() = GPSL5ICache(CircularDeque{Float32}(L5I_WINDOW_SYMBOLS))

function Base.:(==)(a::GPSL5ICache, b::GPSL5ICache)
    deques_equal(a.soft_buffer, b.soft_buffer)
end

function is_ephemeris_decoded(data::GPSL5IData)
    # Message type 10
    !isnothing(data.WN) &&
        !isnothing(data.ΔA) &&
        !isnothing(data.A_dot) &&
        !isnothing(data.Δn_0) &&
        !isnothing(data.M_0) &&
        !isnothing(data.e) &&
        !isnothing(data.ω) &&
        # Message type 11
        !isnothing(data.t_0e) &&
        !isnothing(data.Ω_0) &&
        !isnothing(data.i_0) &&
        !isnothing(data.ΔΩ_dot) &&
        !isnothing(data.i_0_dot) &&
        !isnothing(data.C_is) &&
        !isnothing(data.C_ic) &&
        !isnothing(data.C_rs) &&
        !isnothing(data.C_rc) &&
        !isnothing(data.C_us) &&
        !isnothing(data.C_uc)
end

function is_clock_correction_decoded(data::GPSL5IData)
    # Message types 30-37 shared clock block + message type 30 group delay
    !isnothing(data.t_0c) &&
        !isnothing(data.a_f0) &&
        !isnothing(data.a_f1) &&
        !isnothing(data.a_f2) &&
        !isnothing(data.T_GD)
end

function is_decoding_completed_for_positioning(data::GPSL5IData)
    !isnothing(data.TOW) && is_ephemeris_decoded(data) && is_clock_correction_decoded(data)
end

"""
$(TYPEDSIGNATURES)

Create a decoder state for GPS L5I CNAV navigation messages.

Initializes a [`GNSSDecoderState`](@ref) configured for decoding GPS L5I
civil navigation (CNAV) messages from FEC-encoded 100 sps soft symbols. Each
sync attempt Viterbi-decodes the buffered 616-symbol window, locates the
8-bit preamble (`0b10001011`) at both ends of the decoded bit window,
validates the 300-bit message with CRC-24Q, and dispatches it to per-type
parsers (message types 10-15, 30-37, and 40, IS-GPS-705J §20.3.3).

# Arguments

  - `prn::Int`: Pseudo-Random Noise code identifier (1-63 for GPS satellites)

# Returns

  - `GNSSDecoderState{GPSL5IData}`: Initialized decoder state for GPS L5I

# Example

```julia
state = GPSL5IDecoderState(1)            # PRN 1
state = decode(state, soft_symbols, num_symbols)
if is_sat_healthy(state)
    # Use state.data for positioning
end
```

# See Also

  - [`GNSSDecoderState`](@ref): The underlying state structure
  - [`decode`](@ref): Decode soft symbols using this state
  - [`reset_decoder_state`](@ref): Reset after signal loss
  - [`is_sat_healthy`](@ref): Check satellite health status
"""
function GPSL5IDecoderState(prn)
    GNSSDecoderState(
        prn,
        GPSL5IData(),
        GPSL5IData(),
        GPSL5IConstants(),
        GPSL5ICache(),
        nothing,
        false,
    )
end

function GNSSDecoderState(system::GPSL5I, prn)
    GPSL5IDecoderState(prn)
end

"""
$(TYPEDSIGNATURES)

Reset the GPS L5I decoder state after a signal loss or reacquisition.

Clears the soft-symbol buffer and the time-of-week (TOW) field while
preserving the remaining decoded ephemeris and clock data in `raw_data`, so a
`GNSSReceiver` can re-use the satellite after reacquisition without
re-decoding all message types. Mirrors the semantics of the GPS L1 C/A and
Galileo E1B implementations.

# Arguments

  - `state::GNSSDecoderState{<:GPSL5IData}`: Current GPS L5I decoder state

# Returns

  - `GNSSDecoderState{<:GPSL5IData}`: Reset decoder state with cleared buffers

# See Also

  - [`GPSL5IDecoderState`](@ref): Create a fresh decoder state
  - [`decode`](@ref): Continue decoding after reset
"""
function reset_decoder_state(state::GNSSDecoderState{<:GPSL5IData})
    empty!(state.cache.soft_buffer)
    GNSSDecoderState(
        state;
        raw_data = GPSL5IData(state.raw_data; TOW = nothing),
        data = GPSL5IData(),
        num_bits_after_valid_syncro_sequence = nothing,
        is_shifted_by_180_degrees = false,
    )
end

# ---- Viterbi ----------------------------------------------------------------
#
# The CNAV FEC is the rate-1/2, K=7 NSC code with G1 = 0o171, G2 = 0o133
# (IS-GPS-705J §3.3.3.1.1, Figure 3-7) — the same code Galileo E1B uses, but
# convolved continuously across message boundaries with no tail bits. The
# window decoder below therefore starts from an unknown trellis state (all
# path metrics equal), runs a full add-compare-select forward pass over the
# 616-symbol window, and chains back from the best final state. Noise-free
# the maximum-likelihood path is exact at both window edges; under noise the
# preamble + CRC gate rejects any residual edge errors.

const L5I_VITERBI_NUM_STATES = 64  # 2^(K-1) for K=7

# Encoder output (G1, G2) for every (state, input_bit) pair. State numbering:
# bits [s1 s2 s3 s4 s5 s6] where s1 (MSB) is the oldest register stage. Input
# bit u shifts in as the new s6:
#   y1 = u ⊕ s3 ⊕ s4 ⊕ s5 ⊕ s6   (G1 = 0o171 = 0b1111001)
#   y2 = u ⊕ s1 ⊕ s3 ⊕ s4 ⊕ s6   (G2 = 0o133 = 0b1011011)
# Reference: IS-GPS-705J Figure 3-7.
function _l5i_encoder_output(state::UInt8, u::UInt8)
    s1 = (state >> 5) & 0x01
    s3 = (state >> 3) & 0x01
    s4 = (state >> 2) & 0x01
    s5 = (state >> 1) & 0x01
    s6 = state & 0x01
    y1 = u ⊻ s3 ⊻ s4 ⊻ s5 ⊻ s6
    y2 = u ⊻ s1 ⊻ s3 ⊻ s4 ⊻ s6
    return (y1, y2)
end

const _L5I_VITERBI_G1 = [
    Bool(_l5i_encoder_output(UInt8(s), UInt8(u))[1]) for
    s = 0:(L5I_VITERBI_NUM_STATES-1), u = 0:1
]
const _L5I_VITERBI_G2 = [
    Bool(_l5i_encoder_output(UInt8(s), UInt8(u))[2]) for
    s = 0:(L5I_VITERBI_NUM_STATES-1), u = 0:1
]

# Soft branch penalty for one expected encoder output bit `y` given the
# received LLR (positive ⇒ bit 0): charge `abs(llr)` iff the hard slice
# disagrees with `y`. Summed over a path this equals the maximum-likelihood
# correlation metric up to a path-independent constant, so the decoder is
# exact ML on soft inputs and degrades gracefully to hard-decision Hamming
# distance for ±1 inputs.
@inline _l5i_branch_penalty(llr::Float32, y::Bool) = hard_slice(llr) != y ? abs(llr) : 0.0f0

"""
    gps_l5i_viterbi(soft_window) -> Vector{Bool}

Viterbi-decode a window of FEC soft symbols (LLR convention, length must be
even) into `length(soft_window) ÷ 2` bits. Starts from an unknown encoder
state (the CNAV FEC runs continuously across message boundaries, so the
window's initial state is the tail of the previous message) and traces back
from the best final state.
"""
function gps_l5i_viterbi(soft_window::AbstractVector{Float32})
    num_steps = length(soft_window) ÷ 2
    metrics = zeros(Float32, L5I_VITERBI_NUM_STATES)
    next_metrics = Vector{Float32}(undef, L5I_VITERBI_NUM_STATES)
    # decisions[s+1, t]: whether the survivor into state s at step t came from
    # the odd predecessor (predecessor LSB = 1).
    decisions = Matrix{Bool}(undef, L5I_VITERBI_NUM_STATES, num_steps)
    @inbounds for t = 1:num_steps
        llr1 = soft_window[2t-1]
        llr2 = soft_window[2t]
        for sp = 0:(L5I_VITERBI_NUM_STATES-1)
            # Predecessors of `sp` share their top 5 bits with sp's bottom 5;
            # both are reached via the same input bit u (top bit of sp).
            ps0 = UInt8((sp & 0x1f) << 1)
            ps1 = ps0 | 0x01
            u = ((sp >> 5) & 0x01) + 1
            m0 =
                metrics[ps0+1] +
                _l5i_branch_penalty(llr1, _L5I_VITERBI_G1[ps0+1, u]) +
                _l5i_branch_penalty(llr2, _L5I_VITERBI_G2[ps0+1, u])
            m1 =
                metrics[ps1+1] +
                _l5i_branch_penalty(llr1, _L5I_VITERBI_G1[ps1+1, u]) +
                _l5i_branch_penalty(llr2, _L5I_VITERBI_G2[ps1+1, u])
            from_odd = m1 < m0
            decisions[sp+1, t] = from_odd
            next_metrics[sp+1] = from_odd ? m1 : m0
        end
        # Renormalize so the accumulated Float32 penalties stay small and
        # precision is preserved regardless of the caller's LLR magnitudes.
        min_metric = minimum(next_metrics)
        for i = 1:L5I_VITERBI_NUM_STATES
            metrics[i] = next_metrics[i] - min_metric
        end
    end
    # Chain back from the best final state; the decoded bit at step t is the
    # input bit that produced that step's state (its top bit).
    best_state = argmin(metrics) - 1
    bits = Vector{Bool}(undef, num_steps)
    @inbounds for t = num_steps:-1:1
        bits[t] = (best_state >> 5) & 0x01 == 0x01
        best_state = ((best_state & 0x1f) << 1) | (decisions[best_state+1, t] ? 1 : 0)
    end
    return bits
end

# ---- Sync ------------------------------------------------------------------
#
# Override the generic packed-buffer sync: the CNAV preamble only exists in
# the *decoded* bit domain, so each sync attempt Viterbi-decodes the whole
# 616-symbol window into 308 bits and mirrors `find_preamble` on the result —
# the preamble must be visible at both the start of this message (bits 1-8)
# and the start of the next (bits 301-308), either both upright or both
# inverted (the NSC generators have odd weight, so a 180°-inverted symbol
# stream decodes to the complemented bit stream). The 300-bit message must
# then pass CRC-24Q; preamble matches with a failing CRC are treated as no
# sync, which lets the window slide on by one symbol.

"""
Sync result for GPS L5I: the polarity-resolved 300-bit message and the detected polarity.
"""
struct GPSL5ISync
    """
    Polarity-resolved message, packed MSB-first (bit 1 of the message at the MSB)
    """
    word::UInt320
    """
    Whether the symbol stream is 180-degrees phase shifted
    """
    polarity_flipped::Bool
end

"""
Pack `bits[start:start+7]` MSB-first into a `UInt8`.
"""
function _pack_preamble(bits::Vector{Bool}, start::Int)
    candidate = 0x00
    @inbounds for i = 0:(L5I_PREAMBLE_BITS-1)
        candidate = UInt8(candidate << 1) | UInt8(bits[start+i])
    end
    return candidate
end

"""
    try_sync(state::GNSSDecoderState{<:GPSL5IData}) -> Union{Nothing,GPSL5ISync}

Per-signal sync hook for GPS L5I. Viterbi-decodes the buffered 616-symbol
window, requires the CNAV preamble at both ends of the decoded 308-bit window
(in either polarity) and a clean CRC-24Q over the polarity-resolved 300-bit
message. Returns the [`GPSL5ISync`](@ref) on a match (carrying the message
bits and the detected polarity flip) or `nothing`.
"""
function try_sync(state::GNSSDecoderState{<:GPSL5IData})
    deque = soft_buffer(state)
    window = _deque_slice(deque, 1, L5I_WINDOW_SYMBOLS)
    bits = gps_l5i_viterbi(window)
    leading = _pack_preamble(bits, 1)
    trailing = _pack_preamble(bits, L5I_MESSAGE_BITS + 1)
    if leading == L5I_PREAMBLE && trailing == L5I_PREAMBLE
        polarity_flipped = false
    elseif leading == ~L5I_PREAMBLE && trailing == ~L5I_PREAMBLE
        polarity_flipped = true
        # Inverted symbols decode to complemented bits; un-complement the
        # message portion before the CRC check and field extraction.
        @inbounds for i = 1:L5I_MESSAGE_BITS
            bits[i] = !bits[i]
        end
    else
        return nothing
    end
    crc24q(view(bits, 1:L5I_MESSAGE_BITS)) == 0 || return nothing
    word = UInt320(0)
    @inbounds for i = 1:L5I_MESSAGE_BITS
        word = (word << 1) | UInt320(bits[i])
    end
    return GPSL5ISync(word, polarity_flipped)
end

"""
    complement_buffer_if_necessary(state::GNSSDecoderState{<:GPSL5IData}, sync)

Resolve the 180° polarity ambiguity for GPS L5I. `sync::GPSL5ISync` already
carries the polarity-resolved message bits from `try_sync`; record the
polarity on the state and pass the sync object through unchanged for
`decode_syncro_sequence`.
"""
function complement_buffer_if_necessary(
    state::GNSSDecoderState{<:GPSL5IData},
    sync::GPSL5ISync,
)
    GNSSDecoderState(state; is_shifted_by_180_degrees = sync.polarity_flipped), sync
end

# ---- Message pipeline -------------------------------------------------------

"""
    decode_syncro_sequence(state::GNSSDecoderState{<:GPSL5IData}, sync::GPSL5ISync)

Process one CRC-validated 300-bit CNAV message: parse the common header
(message type ID at bits 15-20, message TOW count at 21-37, alert flag at
38 — IS-GPS-705J §20.3.3), then dispatch to the per-message-type parser.
Unknown or reserved message types keep the decoded header but no further
fields.
"""
function decode_syncro_sequence(state::GNSSDecoderState{<:GPSL5IData}, sync::GPSL5ISync)
    word = sync.word
    word_length = L5I_MESSAGE_BITS
    PI = state.constants.PI

    message_id = Int(get_bits(word, word_length, 15, 6))
    # The message TOW count is the SV time at the start of the *next*
    # 6-second message, in units of 6 seconds.
    TOW = Int64(get_bits(word, word_length, 21, 17)) * 6
    alert_flag = get_bit(word, word_length, 38)
    raw = GPSL5IData(state.raw_data; last_message_id = message_id, TOW, alert_flag)

    raw = if message_id == 10
        parse_mt10(raw, word, PI)
    elseif message_id == 11
        parse_mt11(raw, word, PI)
    elseif message_id == 12
        parse_mt12(raw, word, PI)
    elseif message_id == 13
        parse_mt13(raw, word)
    elseif message_id == 14
        parse_mt14(raw, word, PI)
    elseif message_id == 15
        parse_mt15(raw, word)
    elseif message_id == 30
        parse_mt30(raw, word)
    elseif message_id == 31
        parse_mt31(raw, word, PI)
    elseif message_id == 32
        parse_mt32(raw, word)
    elseif message_id == 33
        parse_mt33(raw, word)
    elseif message_id == 34
        parse_mt34(raw, word, PI)
    elseif message_id == 35
        parse_mt35(raw, word)
    elseif message_id == 36
        parse_mt36(raw, word)
    elseif message_id == 37
        parse_mt37(raw, word, PI)
    elseif message_id == 40
        parse_mt40(raw, word)
    else
        raw  # unknown/reserved message type: header only
    end
    GNSSDecoderState(state; raw_data = raw)
end

# ---- Per-message-type bit-field extraction (IS-GPS-705J §20.3.3) ------------
#
# `word` is the 300-bit message packed MSB-first into a `UInt320` (bit 1 = the
# leftmost preamble bit). Fields are read by 1-based start bit and length
# through the shared `get_bits` / `get_twos_complement_num` / `get_bit`
# helpers, matching the other decoders (#48).

"""
Message type 10 — ephemeris 1 + health (IS-GPS-705J Fig 20-1, Table 20-I).
"""
function parse_mt10(raw::GPSL5IData, word::UInt320, PI::Float64)
    word_length = L5I_MESSAGE_BITS
    GPSL5IData(
        raw;
        WN = Int64(get_bits(word, word_length, 39, 13)),
        l1_health = get_bit(word, word_length, 52),
        l2_health = get_bit(word, word_length, 53),
        l5_health = get_bit(word, word_length, 54),
        t_op = Int64(get_bits(word, word_length, 55, 11)) * 300,
        ura_ed_index = get_twos_complement_num(word, word_length, 66, 5),
        t_0e = Int64(get_bits(word, word_length, 71, 11)) * 300,
        ΔA = get_twos_complement_num(word, word_length, 82, 26) * 2.0^-9,
        A_dot = get_twos_complement_num(word, word_length, 108, 25) * 2.0^-21,
        Δn_0 = get_twos_complement_num(word, word_length, 133, 17) * 2.0^-44 * PI,
        Δn_0_dot = get_twos_complement_num(word, word_length, 150, 23) * 2.0^-57 * PI,
        M_0 = get_twos_complement_num(word, word_length, 173, 33) * 2.0^-32 * PI,
        e = Int(get_bits(word, word_length, 206, 33)) * 2.0^-34,
        ω = get_twos_complement_num(word, word_length, 239, 33) * 2.0^-32 * PI,
        integrity_status_flag = get_bit(word, word_length, 272),
        l2c_phasing = get_bit(word, word_length, 273),
    )
end

"""
Message type 11 — ephemeris 2 (IS-GPS-705J Fig 20-2, Table 20-I).
"""
function parse_mt11(raw::GPSL5IData, word::UInt320, PI::Float64)
    word_length = L5I_MESSAGE_BITS
    GPSL5IData(
        raw;
        t_0e = Int64(get_bits(word, word_length, 39, 11)) * 300,
        Ω_0 = get_twos_complement_num(word, word_length, 50, 33) * 2.0^-32 * PI,
        i_0 = get_twos_complement_num(word, word_length, 83, 33) * 2.0^-32 * PI,
        ΔΩ_dot = get_twos_complement_num(word, word_length, 116, 17) * 2.0^-44 * PI,
        i_0_dot = get_twos_complement_num(word, word_length, 133, 15) * 2.0^-44 * PI,
        C_is = get_twos_complement_num(word, word_length, 148, 16) * 2.0^-30,
        C_ic = get_twos_complement_num(word, word_length, 164, 16) * 2.0^-30,
        C_rs = get_twos_complement_num(word, word_length, 180, 24) * 2.0^-8,
        C_rc = get_twos_complement_num(word, word_length, 204, 24) * 2.0^-8,
        C_us = get_twos_complement_num(word, word_length, 228, 21) * 2.0^-30,
        C_uc = get_twos_complement_num(word, word_length, 249, 21) * 2.0^-30,
    )
end

"""
Shared clock block of message types 30-37, bits 39-127 (IS-GPS-705J Fig 20-3, Table 20-III).
"""
function parse_clock_block(raw::GPSL5IData, word::UInt320)
    word_length = L5I_MESSAGE_BITS
    GPSL5IData(
        raw;
        t_op = Int64(get_bits(word, word_length, 39, 11)) * 300,
        ura_ned0_index = get_twos_complement_num(word, word_length, 50, 5),
        ura_ned1_index = Int64(get_bits(word, word_length, 55, 3)),
        ura_ned2_index = Int64(get_bits(word, word_length, 58, 3)),
        t_0c = Int64(get_bits(word, word_length, 61, 11)) * 300,
        a_f0 = get_twos_complement_num(word, word_length, 72, 26) * 2.0^-35,
        a_f1 = get_twos_complement_num(word, word_length, 98, 20) * 2.0^-48,
        a_f2 = get_twos_complement_num(word, word_length, 118, 10) * 2.0^-60,
    )
end

"""
Message type 30 — clock, iono & group delay (IS-GPS-705J Fig 20-3, Tables 20-III/20-IV).
"""
function parse_mt30(raw::GPSL5IData, word::UInt320)
    word_length = L5I_MESSAGE_BITS
    raw = parse_clock_block(raw, word)
    GPSL5IData(
        raw;
        T_GD = get_twos_complement_num(word, word_length, 128, 13) * 2.0^-35,
        ISC_L1CA = get_twos_complement_num(word, word_length, 141, 13) * 2.0^-35,
        ISC_L2C = get_twos_complement_num(word, word_length, 154, 13) * 2.0^-35,
        ISC_L5I5 = get_twos_complement_num(word, word_length, 167, 13) * 2.0^-35,
        ISC_L5Q5 = get_twos_complement_num(word, word_length, 180, 13) * 2.0^-35,
        α0 = get_twos_complement_num(word, word_length, 193, 8) * 2.0^-30,
        α1 = get_twos_complement_num(word, word_length, 201, 8) * 2.0^-27,
        α2 = get_twos_complement_num(word, word_length, 209, 8) * 2.0^-24,
        α3 = get_twos_complement_num(word, word_length, 217, 8) * 2.0^-24,
        β0 = get_twos_complement_num(word, word_length, 225, 8) * 2.0^11,
        β1 = get_twos_complement_num(word, word_length, 233, 8) * 2.0^14,
        β2 = get_twos_complement_num(word, word_length, 241, 8) * 2.0^16,
        β3 = get_twos_complement_num(word, word_length, 249, 8) * 2.0^16,
        WN_op = Int64(get_bits(word, word_length, 257, 8)),
    )
end

"""
Decode one 31-bit reduced-almanac packet starting at 1-based bit `start` (IS-GPS-705J Fig 20-16).
"""
function _l5i_reduced_almanac_packet(
    word::UInt320,
    start::Int,
    WN_a::Int,
    t_oa::Int,
    PI::Float64,
)
    word_length = L5I_MESSAGE_BITS
    PRN_a = Int(get_bits(word, word_length, start, 6))
    PRN_a == 0 && return nothing  # empty packet ⇒ no further packets follow
    GPSL5IReducedAlmanac(;
        PRN_a,
        WN_a,
        t_oa,
        δA = get_twos_complement_num(word, word_length, start + 6, 8) * 2.0^9,
        Ω_0 = get_twos_complement_num(word, word_length, start + 14, 7) * 2.0^-6 * PI,
        Φ_0 = get_twos_complement_num(word, word_length, start + 21, 7) * 2.0^-6 * PI,
        l1_health = get_bit(word, word_length, start + 28),
        l2_health = get_bit(word, word_length, start + 29),
        l5_health = get_bit(word, word_length, start + 30),
    )
end

"""
Message type 12 — seven reduced-almanac packets (IS-GPS-705J Fig 20-11).
"""
function parse_mt12(raw::GPSL5IData, word::UInt320, PI::Float64)
    word_length = L5I_MESSAGE_BITS
    WN_a = Int(get_bits(word, word_length, 39, 13))
    t_oa = Int(get_bits(word, word_length, 52, 8)) * 2^12
    almanacs = raw.reduced_almanacs
    # Seven 31-bit packets at bits 60, 91, 122, 153, 184, 215, 246.
    for start in (60, 91, 122, 153, 184, 215, 246)
        packet = _l5i_reduced_almanac_packet(word, start, WN_a, t_oa, PI)
        isnothing(packet) && break  # PRN_a==0 ⇒ remaining packets are filler
        almanacs = _merge_keyed(almanacs, packet.PRN_a, packet)
    end
    GPSL5IData(raw; reduced_almanacs = almanacs)
end

"""
Message type 31 — clock & four reduced-almanac packets (IS-GPS-705J Fig 20-4).
"""
function parse_mt31(raw::GPSL5IData, word::UInt320, PI::Float64)
    word_length = L5I_MESSAGE_BITS
    raw = parse_clock_block(raw, word)
    WN_a = Int(get_bits(word, word_length, 128, 13))
    t_oa = Int(get_bits(word, word_length, 141, 8)) * 2^12
    almanacs = raw.reduced_almanacs
    # Four 31-bit packets at bits 149, 180, 211, 242.
    for start in (149, 180, 211, 242)
        packet = _l5i_reduced_almanac_packet(word, start, WN_a, t_oa, PI)
        isnothing(packet) && break  # PRN_a==0 ⇒ remaining packets are filler
        almanacs = _merge_keyed(almanacs, packet.PRN_a, packet)
    end
    GPSL5IData(raw; reduced_almanacs = almanacs)
end

"""
Message type 32 — clock & EOP (IS-GPS-705J Fig 20-5, Table 20-VII).
"""
function parse_mt32(raw::GPSL5IData, word::UInt320)
    word_length = L5I_MESSAGE_BITS
    raw = parse_clock_block(raw, word)
    GPSL5IData(
        raw;
        t_EOP = Int64(get_bits(word, word_length, 128, 16)) * 16,
        PM_X = get_twos_complement_num(word, word_length, 144, 21) * 2.0^-20,
        PM_X_dot = get_twos_complement_num(word, word_length, 165, 15) * 2.0^-21,
        PM_Y = get_twos_complement_num(word, word_length, 180, 21) * 2.0^-20,
        PM_Y_dot = get_twos_complement_num(word, word_length, 201, 15) * 2.0^-21,
        ΔUT_GPS = get_twos_complement_num(word, word_length, 216, 31) * 2.0^-23,
        ΔUT_GPS_dot = get_twos_complement_num(word, word_length, 247, 19) * 2.0^-25,
    )
end

"""
Message type 33 — clock & UTC (IS-GPS-705J Fig 20-6, Table 20-IX).
"""
function parse_mt33(raw::GPSL5IData, word::UInt320)
    word_length = L5I_MESSAGE_BITS
    raw = parse_clock_block(raw, word)
    GPSL5IData(
        raw;
        A0_UTC = get_twos_complement_num(word, word_length, 128, 16) * 2.0^-35,
        A1_UTC = get_twos_complement_num(word, word_length, 144, 13) * 2.0^-51,
        A2_UTC = get_twos_complement_num(word, word_length, 157, 7) * 2.0^-68,
        Δt_LS = get_twos_complement_num(word, word_length, 164, 8),
        t_ot = Int64(get_bits(word, word_length, 172, 16)) * 16,
        WN_ot = Int64(get_bits(word, word_length, 188, 13)),
        WN_LSF = Int64(get_bits(word, word_length, 201, 13)),
        DN = Int64(get_bits(word, word_length, 214, 4)),
        Δt_LSF = get_twos_complement_num(word, word_length, 218, 8),
    )
end

"""
Message type 35 — clock & GGTO (IS-GPS-705J Fig 20-8, Table 20-XI).
"""
function parse_mt35(raw::GPSL5IData, word::UInt320)
    word_length = L5I_MESSAGE_BITS
    raw = parse_clock_block(raw, word)
    GPSL5IData(
        raw;
        t_GGTO = Int64(get_bits(word, word_length, 128, 16)) * 16,
        WN_GGTO = Int64(get_bits(word, word_length, 144, 13)),
        GNSS_ID = Int64(get_bits(word, word_length, 157, 3)),
        A0_GGTO = get_twos_complement_num(word, word_length, 160, 16) * 2.0^-35,
        A1_GGTO = get_twos_complement_num(word, word_length, 176, 13) * 2.0^-51,
        A2_GGTO = get_twos_complement_num(word, word_length, 189, 7) * 2.0^-68,
    )
end

"""
Message type 37 — clock & one Midi almanac (IS-GPS-705J Fig 20-10, Table 20-V).
"""
function parse_mt37(raw::GPSL5IData, word::UInt320, PI::Float64)
    word_length = L5I_MESSAGE_BITS
    raw = parse_clock_block(raw, word)
    PRN_a = Int(get_bits(word, word_length, 149, 6))
    PRN_a == 0 && return raw  # empty almanac
    alm = GPSL5IMidiAlmanac(;
        PRN_a,
        WN_a = Int(get_bits(word, word_length, 128, 13)),
        t_oa = Int(get_bits(word, word_length, 141, 8)) * 2^12,
        l1_health = get_bit(word, word_length, 155),
        l2_health = get_bit(word, word_length, 156),
        l5_health = get_bit(word, word_length, 157),
        e = Int(get_bits(word, word_length, 158, 11)) * 2.0^-16,
        δi = get_twos_complement_num(word, word_length, 169, 11) * 2.0^-14 * PI,
        Ω_dot = get_twos_complement_num(word, word_length, 180, 11) * 2.0^-33 * PI,
        sqrt_A = Int(get_bits(word, word_length, 191, 17)) * 2.0^-4,
        Ω_0 = get_twos_complement_num(word, word_length, 208, 16) * 2.0^-15 * PI,
        ω = get_twos_complement_num(word, word_length, 224, 16) * 2.0^-15 * PI,
        M_0 = get_twos_complement_num(word, word_length, 240, 16) * 2.0^-15 * PI,
        a_f0 = get_twos_complement_num(word, word_length, 256, 11) * 2.0^-20,
        a_f1 = get_twos_complement_num(word, word_length, 267, 10) * 2.0^-37,
    )
    GPSL5IData(raw; midi_almanacs = _merge_keyed(raw.midi_almanacs, PRN_a, alm))
end

# All-ones PRN ID in a CDC/EDC packet ⇒ no DC data in the remainder of the
# data block (IS-GPS-705J §20.3.3.7.2.3).
const L5I_DC_EMPTY_PRN = 0xff

"""
Decode one 34-bit CDC packet starting at 1-based bit `start` (IS-GPS-705J Fig 20-17).
"""
function _l5i_cdc_packet(
    word::UInt320,
    start::Int,
    t_op_D::Int,
    t_OD::Int,
    dc_data_type::Bool,
)
    word_length = L5I_MESSAGE_BITS
    PRN_a = Int(get_bits(word, word_length, start, 8))
    PRN_a == L5I_DC_EMPTY_PRN && return nothing
    GPSL5IClockDifferentialCorrection(;
        PRN_a,
        t_op_D,
        t_OD,
        dc_data_type,
        δa_f0 = get_twos_complement_num(word, word_length, start + 8, 13) * 2.0^-35,
        δa_f1 = get_twos_complement_num(word, word_length, start + 21, 8) * 2.0^-51,
        UDRA_index = get_twos_complement_num(word, word_length, start + 29, 5),
    )
end

"""
Decode one 92-bit EDC packet starting at 1-based bit `start` (IS-GPS-705J Fig 20-17).
"""
function _l5i_edc_packet(
    word::UInt320,
    start::Int,
    t_op_D::Int,
    t_OD::Int,
    dc_data_type::Bool,
    PI::Float64,
)
    word_length = L5I_MESSAGE_BITS
    PRN_a = Int(get_bits(word, word_length, start, 8))
    PRN_a == L5I_DC_EMPTY_PRN && return nothing
    GPSL5IEphemerisDifferentialCorrection(;
        PRN_a,
        t_op_D,
        t_OD,
        dc_data_type,
        Δα = get_twos_complement_num(word, word_length, start + 8, 14) * 2.0^-34,
        Δβ = get_twos_complement_num(word, word_length, start + 22, 14) * 2.0^-34,
        Δγ = get_twos_complement_num(word, word_length, start + 36, 15) * 2.0^-32 * PI,
        Δi = get_twos_complement_num(word, word_length, start + 51, 12) * 2.0^-32 * PI,
        ΔΩ = get_twos_complement_num(word, word_length, start + 63, 12) * 2.0^-32 * PI,
        ΔA = get_twos_complement_num(word, word_length, start + 75, 12) * 2.0^-9,
        UDRA_dot_index = get_twos_complement_num(word, word_length, start + 87, 5),
    )
end

"""
Message type 13 — clock differential correction (IS-GPS-705J Fig 20-12).
"""
function parse_mt13(raw::GPSL5IData, word::UInt320)
    word_length = L5I_MESSAGE_BITS
    t_op_D = Int(get_bits(word, word_length, 39, 11)) * 300
    t_OD = Int(get_bits(word, word_length, 50, 11)) * 300
    corrections = raw.clock_corrections
    # Six 35-bit packets (1 DC data type bit + 34-bit CDC) at bits 61, 96,
    # 131, 166, 201, 236.
    for start in (61, 96, 131, 166, 201, 236)
        dc_data_type = get_bit(word, word_length, start)
        packet = _l5i_cdc_packet(word, start + 1, t_op_D, t_OD, dc_data_type)
        isnothing(packet) && continue
        corrections = _merge_keyed(corrections, packet.PRN_a, packet)
    end
    GPSL5IData(raw; clock_corrections = corrections)
end

"""
Message type 14 — ephemeris differential correction (IS-GPS-705J Fig 20-13).
"""
function parse_mt14(raw::GPSL5IData, word::UInt320, PI::Float64)
    word_length = L5I_MESSAGE_BITS
    t_op_D = Int(get_bits(word, word_length, 39, 11)) * 300
    t_OD = Int(get_bits(word, word_length, 50, 11)) * 300
    corrections = raw.ephemeris_corrections
    # Two 93-bit packets (1 DC data type bit + 92-bit EDC) at bits 61, 154.
    for start in (61, 154)
        dc_data_type = get_bit(word, word_length, start)
        packet = _l5i_edc_packet(word, start + 1, t_op_D, t_OD, dc_data_type, PI)
        isnothing(packet) && continue
        corrections = _merge_keyed(corrections, packet.PRN_a, packet)
    end
    GPSL5IData(raw; ephemeris_corrections = corrections)
end

"""
Message type 34 — clock & one CDC+EDC pair (IS-GPS-705J Fig 20-7).
"""
function parse_mt34(raw::GPSL5IData, word::UInt320, PI::Float64)
    word_length = L5I_MESSAGE_BITS
    raw = parse_clock_block(raw, word)
    t_op_D = Int(get_bits(word, word_length, 128, 11)) * 300
    t_OD = Int(get_bits(word, word_length, 139, 11)) * 300
    dc_data_type = get_bit(word, word_length, 150)
    cdc = _l5i_cdc_packet(word, 151, t_op_D, t_OD, dc_data_type)
    edc = _l5i_edc_packet(word, 185, t_op_D, t_OD, dc_data_type, PI)
    clock_corrections =
        isnothing(cdc) ? raw.clock_corrections :
        _merge_keyed(raw.clock_corrections, cdc.PRN_a, cdc)
    ephemeris_corrections =
        isnothing(edc) ? raw.ephemeris_corrections :
        _merge_keyed(raw.ephemeris_corrections, edc.PRN_a, edc)
    GPSL5IData(raw; clock_corrections, ephemeris_corrections)
end

"""
Decode `num_chars` 8-bit ASCII characters starting at 1-based bit `start`, stripping control chars.
"""
function _l5i_text(word::UInt320, start::Int, num_chars::Int)
    word_length = L5I_MESSAGE_BITS
    chars = Char[]
    for k = 0:(num_chars-1)
        code = Int(get_bits(word, word_length, start + 8k, 8))
        # Keep printable ASCII; skip NUL/control padding so the message is clean.
        (code >= 0x20 && code < 0x7f) && push!(chars, Char(code))
    end
    return String(chars)
end

"""
Message type 15 — 29 ASCII characters at bits 39-270 (IS-GPS-705J Fig 20-14).
"""
function parse_mt15(raw::GPSL5IData, word::UInt320)
    word_length = L5I_MESSAGE_BITS
    GPSL5IData(
        raw;
        text_mt15 = _l5i_text(word, 39, 29),
        text_page_mt15 = Int64(get_bits(word, word_length, 271, 4)),
    )
end

"""
Message type 36 — clock & 18 ASCII characters at bits 128-271 (IS-GPS-705J Fig 20-9).
"""
function parse_mt36(raw::GPSL5IData, word::UInt320)
    word_length = L5I_MESSAGE_BITS
    raw = parse_clock_block(raw, word)
    GPSL5IData(
        raw;
        text_mt36 = _l5i_text(word, 128, 18),
        text_page_mt36 = Int64(get_bits(word, word_length, 272, 4)),
    )
end

"""
Message type 40 — Integrity Support Message (IS-GPS-705J Fig 20-14a).
"""
function parse_mt40(raw::GPSL5IData, word::UInt320)
    word_length = L5I_MESSAGE_BITS
    # 63-bit SV mask split across the bit-89 boundary (12 + 51 bits).
    mask = UInt64(get_bits(word, word_length, 89, 63))
    ism = GPSL5IIntegritySupportMessage(;
        GNSS_ID = Int(get_bits(word, word_length, 39, 4)),
        WN_ISM = Int(get_bits(word, word_length, 43, 13)),
        TOW_ISM = Int(get_bits(word, word_length, 56, 6)),
        t_correl = Int(get_bits(word, word_length, 62, 4)),
        b_nom = Int(get_bits(word, word_length, 66, 4)),
        γ_nom = Int(get_bits(word, word_length, 70, 4)),
        R_sat = Int(get_bits(word, word_length, 74, 4)),
        P_const = Int(get_bits(word, word_length, 78, 4)),
        MFD = Int(get_bits(word, word_length, 82, 4)),
        service_level = Int(get_bits(word, word_length, 86, 3)),
        mask,
    )
    GPSL5IData(raw; ism)
end

"""
    validate_data(state::GNSSDecoderState{<:GPSL5IData})

Promote `raw_data` to `data` once the minimum positioning set (message types
10 + 11 ephemeris and the message type 30 clock + group delay block) has been
decoded. Every CNAV message carries its own TOW, so each validated message
re-arms the streaming counter: the TOW refers to the start of the *next*
message, whose first `preamble_length` symbols are already buffered.
"""
function validate_data(state::GNSSDecoderState{<:GPSL5IData})
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

Check if the GPS L5 satellite is healthy and usable for positioning.

Examines the L5 signal health bit decoded from the most recent message
type 10 (IS-GPS-705J §20.3.3.1.1.2): a satellite is healthy iff the health
bit is 0 (all navigation data on the L5 signal are OK).

!!! warning

    Requires message type 10 to have been decoded and the positioning set to
    have been validated; returns `false` until then.

# Arguments

  - `state::GNSSDecoderState{<:GPSL5IData}`: GPS L5I decoder state.

# Returns

  - `Bool`: `true` iff the L5 signal-health bit indicates OK.
"""
function is_sat_healthy(state::GNSSDecoderState{<:GPSL5IData})
    state.data.l5_health === false
end
