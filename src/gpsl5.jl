# UInt320 buffer for GPS L5 CNAV decoded-bit buffer.
# Holds at least one full 300-bit CNAV message plus 8 bits of preamble lookahead.
# UInt320 is already defined by gpsl1.jl.

# UInt256 holds the largest CNAV text payload (MT 15 has 236 bits of text + flags).
BitIntegers.@define_integers 256

"""
    GPSL5Constants

WGS 84 constants and CNAV message structure parameters for GPS L5 signal decoding.

The physical constants and message structure are defined in IS-GPS-705 (Interface
Specification) and are used for computing satellite positions and clock corrections
from broadcast ephemeris data.

# Fields
- `message_length::Int`: Length of one CNAV message in bits (300)
- `crc_length::Int`: Length of CRC-24Q parity in bits (24)
- `preamble::UInt8`: CNAV preamble pattern (10001011 binary, 0x8B)
- `preamble_length::Int`: Length of preamble in bits (8)
- `viterbi_constraint::Int`: Convolutional code constraint length (7)
- `viterbi_poly_g1::UInt8`: Generator polynomial G1, octal 171 = 0b1111001
- `viterbi_poly_g2::UInt8`: Generator polynomial G2, octal 133 = 0b1011011
- `viterbi_traceback::Int`: Viterbi traceback depth in bits (32, ~5 × constraint length)
- `A_REF::Float64`: Reference semi-major axis (26 559 710 m, IS-GPS-705 Table 20-I)
- `Ω_dot_REF::Float64`: Reference rate of right ascension (-2.6×10⁻⁹ semi-circles/s)
- `PI::Float64`: Mathematical constant π = 3.1415926535898 (IS-GPS-705 Table 20-II)
- `Ω_dot_e::Float64`: WGS 84 Earth rotation rate = 7.2921151467×10⁻⁵ rad/s
- `c::Float64`: Speed of light = 2.99792458×10⁸ m/s
- `μ::Float64`: WGS 84 Earth gravitational parameter = 3.986005×10¹⁴ m³/s²
- `F::Float64`: Relativistic correction constant = -4.442807633×10⁻¹⁰ s/√m

# Reference
IS-GPS-705J, Section 20.3
"""
Base.@kwdef struct GPSL5Constants <: AbstractGNSSConstants
    syncro_sequence_length::Int = 300  # Compatibility with generic field name; equals message_length
    message_length::Int = 300
    crc_length::Int = 24
    preamble::UInt8 = 0b10001011
    preamble_length::Int = 8
    viterbi_constraint::Int = 7
    viterbi_poly_g1::UInt8 = 0b1111001  # 171 octal
    viterbi_poly_g2::UInt8 = 0b1011011  # 133 octal
    viterbi_traceback::Int = 32
    A_REF::Float64 = 26_559_710.0
    Ω_dot_REF::Float64 = -2.6e-9
    PI::Float64 = 3.1415926535898
    Ω_dot_e::Float64 = 7.2921151467e-5
    c::Float64 = 2.99792458e8
    μ::Float64 = 3.986005e14
    F::Float64 = -4.442807633e-10
end

"""
    GPSL5ReducedAlmanacPacket

Reduced almanac packet decoded from CNAV message types 12 or 31.

# Fields
- `PRN_a::UInt`: Satellite PRN (1-63; 0 indicates an unused packet)
- `δA::Float64`: Difference of semi-major axis from `A_REF` (meters)
- `Ω_0::Float64`: Longitude of ascending node (semi-circles)
- `Φ_0::Float64`: Argument of latitude at reference time = M0 + ω (semi-circles)
- `signal_health_l1::Bool`: L1 signal health (false = OK, true = bad)
- `signal_health_l2::Bool`: L2 signal health
- `signal_health_l5::Bool`: L5 signal health

# Reference
IS-GPS-705J, Figure 20-16, Table 20-VI
"""
Base.@kwdef struct GPSL5ReducedAlmanacPacket
    PRN_a::Union{Nothing,UInt} = nothing
    δA::Union{Nothing,Float64} = nothing
    Ω_0::Union{Nothing,Float64} = nothing
    Φ_0::Union{Nothing,Float64} = nothing
    signal_health_l1::Union{Nothing,Bool} = nothing
    signal_health_l2::Union{Nothing,Bool} = nothing
    signal_health_l5::Union{Nothing,Bool} = nothing
end

"""
    GPSL5MidiAlmanac

Midi almanac data decoded from CNAV message type 37.

# Reference
IS-GPS-705J, Figure 20-10, Table 20-V
"""
Base.@kwdef struct GPSL5MidiAlmanac
    PRN_a::Union{Nothing,UInt} = nothing
    signal_health_l1::Union{Nothing,Bool} = nothing
    signal_health_l2::Union{Nothing,Bool} = nothing
    signal_health_l5::Union{Nothing,Bool} = nothing
    e::Union{Nothing,Float64} = nothing
    δi::Union{Nothing,Float64} = nothing
    Ω_dot::Union{Nothing,Float64} = nothing
    sqrt_A::Union{Nothing,Float64} = nothing
    Ω_0::Union{Nothing,Float64} = nothing
    ω::Union{Nothing,Float64} = nothing
    M_0::Union{Nothing,Float64} = nothing
    a_f0::Union{Nothing,Float64} = nothing
    a_f1::Union{Nothing,Float64} = nothing
end

"""
    GPSL5ClockDifferentialCorrection

Clock differential correction packet decoded from CNAV message types 13 or 34.

# Reference
IS-GPS-705J, Figure 20-17 (CDC), Table 20-X
"""
Base.@kwdef struct GPSL5ClockDifferentialCorrection
    PRN_id::Union{Nothing,UInt} = nothing
    δa_f0::Union{Nothing,Float64} = nothing
    δa_f1::Union{Nothing,Float64} = nothing
    UDRA::Union{Nothing,Int} = nothing
end

"""
    GPSL5EphemerisDifferentialCorrection

Ephemeris differential correction packet decoded from CNAV message types 14 or 34.

# Reference
IS-GPS-705J, Figure 20-17 (EDC), Table 20-X
"""
Base.@kwdef struct GPSL5EphemerisDifferentialCorrection
    PRN_id::Union{Nothing,UInt} = nothing
    Δα::Union{Nothing,Float64} = nothing
    Δβ::Union{Nothing,Float64} = nothing
    Δγ::Union{Nothing,Float64} = nothing
    Δi::Union{Nothing,Float64} = nothing
    ΔΩ::Union{Nothing,Float64} = nothing
    ΔA::Union{Nothing,Float64} = nothing
    UDRA_dot::Union{Nothing,Int} = nothing
end

"""
    GPSL5IntegritySupportMessage

Integrity Support Message decoded from CNAV message type 40 (ARAIM).

# Reference
IS-GPS-705J, Figure 20-14a, Table 20-XIa
"""
Base.@kwdef struct GPSL5IntegritySupportMessage
    GNSS_id::Union{Nothing,UInt} = nothing
    WN_ISM::Union{Nothing,UInt} = nothing
    TOW_ISM::Union{Nothing,UInt} = nothing
    t_correl::Union{Nothing,UInt} = nothing
    b_nom::Union{Nothing,UInt} = nothing
    γ_nom::Union{Nothing,UInt} = nothing
    R_sat::Union{Nothing,UInt} = nothing
    P_const::Union{Nothing,UInt} = nothing
    MFD::Union{Nothing,UInt} = nothing
    service_level::Union{Nothing,UInt} = nothing
    mask::Union{Nothing,UInt128} = nothing
end

"""
    GPSL5Data

Decoded GPS L5 CNAV navigation message data.

Contains all parameters that may be decoded from L5 CNAV message types 10, 11,
12, 13, 14, 15, 30-37, and 40 as defined in IS-GPS-705J. The decoder fills
fields incrementally as the corresponding message types are received.

# CNAV Common Header
- `last_message_id::Int`: Most recently decoded message type (10-15, 30-37, 40)
- `TOW::Int64`: Time of Week at start of next 6-second message (seconds)
- `alert_flag::Bool`: Alert flag (1 indicates URA may be worse than indicated)

# Ephemeris (MT 10 + 11)
See IS-GPS-705J Table 20-I.

# Clock + Ionospheric + Group Delay (MT 30)
See IS-GPS-705J Tables 20-III and 20-IV.

# UTC Parameters (MT 33)
See IS-GPS-705J Table 20-IX.

# GPS/GNSS Time Offset (MT 35)
See IS-GPS-705J Table 20-XI.

# Earth Orientation Parameters (MT 32)
See IS-GPS-705J Table 20-VII.

# Almanacs (MT 12, 31, 37)
See IS-GPS-705J Tables 20-V and 20-VI.

# Differential Correction (MT 13, 14, 34)
See IS-GPS-705J Table 20-X.

# Text (MT 15, 36)
ASCII text bits stored verbatim. Caller must convert to characters.

# Integrity Support Message (MT 40)
See IS-GPS-705J Table 20-XIa.

# Reference
IS-GPS-705J, Section 20.3.3
"""
Base.@kwdef struct GPSL5Data <: AbstractGNSSData
    last_message_id::Int = 0
    TOW::Union{Nothing,Int64} = nothing
    alert_flag::Union{Nothing,Bool} = nothing

    # MT 10 - Ephemeris 1 + Health
    WN::Union{Nothing,Int64} = nothing
    signal_health_l1::Union{Nothing,Bool} = nothing
    signal_health_l2::Union{Nothing,Bool} = nothing
    signal_health_l5::Union{Nothing,Bool} = nothing
    t_op::Union{Nothing,Int64} = nothing
    URA_ED::Union{Nothing,Int} = nothing
    ΔA::Union{Nothing,Float64} = nothing
    A_dot::Union{Nothing,Float64} = nothing
    Δn_0::Union{Nothing,Float64} = nothing
    Δn_0_dot::Union{Nothing,Float64} = nothing
    M_0::Union{Nothing,Float64} = nothing
    e::Union{Nothing,Float64} = nothing
    ω::Union{Nothing,Float64} = nothing
    integrity_status_flag::Union{Nothing,Bool} = nothing
    l2c_phasing::Union{Nothing,Bool} = nothing

    # MT 11 - Ephemeris 2
    t_oe::Union{Nothing,Int64} = nothing
    Ω_0::Union{Nothing,Float64} = nothing
    i_0::Union{Nothing,Float64} = nothing
    Ω_dot::Union{Nothing,Float64} = nothing  # Δ(Ω̇) before adding the reference
    i_dot::Union{Nothing,Float64} = nothing
    C_is::Union{Nothing,Float64} = nothing
    C_ic::Union{Nothing,Float64} = nothing
    C_rs::Union{Nothing,Float64} = nothing
    C_rc::Union{Nothing,Float64} = nothing
    C_us::Union{Nothing,Float64} = nothing
    C_uc::Union{Nothing,Float64} = nothing

    # MT 30 - Clock + ISC + IONO
    t_oc::Union{Nothing,Int64} = nothing
    a_f0::Union{Nothing,Float64} = nothing
    a_f1::Union{Nothing,Float64} = nothing
    a_f2::Union{Nothing,Float64} = nothing
    URA_NED0::Union{Nothing,Int} = nothing
    URA_NED1::Union{Nothing,Int} = nothing
    URA_NED2::Union{Nothing,Int} = nothing
    T_GD::Union{Nothing,Float64} = nothing
    ISC_L1CA::Union{Nothing,Float64} = nothing
    ISC_L2C::Union{Nothing,Float64} = nothing
    ISC_L5I5::Union{Nothing,Float64} = nothing
    ISC_L5Q5::Union{Nothing,Float64} = nothing
    α_0::Union{Nothing,Float64} = nothing
    α_1::Union{Nothing,Float64} = nothing
    α_2::Union{Nothing,Float64} = nothing
    α_3::Union{Nothing,Float64} = nothing
    β_0::Union{Nothing,Float64} = nothing
    β_1::Union{Nothing,Float64} = nothing
    β_2::Union{Nothing,Float64} = nothing
    β_3::Union{Nothing,Float64} = nothing
    WN_op::Union{Nothing,UInt} = nothing

    # MT 32 - EOP
    t_EOP::Union{Nothing,Int64} = nothing
    PM_X::Union{Nothing,Float64} = nothing
    PM_X_dot::Union{Nothing,Float64} = nothing
    PM_Y::Union{Nothing,Float64} = nothing
    PM_Y_dot::Union{Nothing,Float64} = nothing
    ΔUT_GPS::Union{Nothing,Float64} = nothing
    ΔUT_GPS_dot::Union{Nothing,Float64} = nothing

    # MT 33 - UTC
    A_0_utc::Union{Nothing,Float64} = nothing
    A_1_utc::Union{Nothing,Float64} = nothing
    A_2_utc::Union{Nothing,Float64} = nothing
    Δt_LS::Union{Nothing,Int} = nothing
    t_ot::Union{Nothing,Int} = nothing
    WN_ot::Union{Nothing,UInt} = nothing
    WN_LSF::Union{Nothing,UInt} = nothing
    DN::Union{Nothing,UInt} = nothing
    Δt_LSF::Union{Nothing,Int} = nothing

    # MT 35 - GGTO
    t_GGTO::Union{Nothing,Int} = nothing
    WN_GGTO::Union{Nothing,UInt} = nothing
    GNSS_id::Union{Nothing,UInt} = nothing
    A_0_GGTO::Union{Nothing,Float64} = nothing
    A_1_GGTO::Union{Nothing,Float64} = nothing
    A_2_GGTO::Union{Nothing,Float64} = nothing

    # MT 12 / 31 / 37 - Almanacs
    WN_a::Union{Nothing,UInt} = nothing
    t_oa::Union{Nothing,Int} = nothing
    reduced_almanacs::Vector{GPSL5ReducedAlmanacPacket} = [GPSL5ReducedAlmanacPacket() for _ = 1:63]
    midi_almanacs::Vector{GPSL5MidiAlmanac} = [GPSL5MidiAlmanac() for _ = 1:63]

    # MT 13 / 14 / 34 - Differential Correction
    t_op_D::Union{Nothing,Int64} = nothing
    t_OD::Union{Nothing,Int64} = nothing
    cdc_packets::Vector{GPSL5ClockDifferentialCorrection} = GPSL5ClockDifferentialCorrection[]
    edc_packets::Vector{GPSL5EphemerisDifferentialCorrection} = GPSL5EphemerisDifferentialCorrection[]

    # MT 15 / 36 - Text (raw bits, callers can decode the ASCII themselves)
    text15_bits::Union{Nothing,UInt256} = nothing
    text15_page::Union{Nothing,UInt} = nothing
    text36_bits::Union{Nothing,UInt256} = nothing
    text36_page::Union{Nothing,UInt} = nothing

    # MT 40 - Integrity Support Message
    ism::GPSL5IntegritySupportMessage = GPSL5IntegritySupportMessage()
end

function GPSL5Data(
    data::GPSL5Data;
    last_message_id = data.last_message_id,
    TOW = data.TOW,
    alert_flag = data.alert_flag,
    WN = data.WN,
    signal_health_l1 = data.signal_health_l1,
    signal_health_l2 = data.signal_health_l2,
    signal_health_l5 = data.signal_health_l5,
    t_op = data.t_op,
    URA_ED = data.URA_ED,
    ΔA = data.ΔA,
    A_dot = data.A_dot,
    Δn_0 = data.Δn_0,
    Δn_0_dot = data.Δn_0_dot,
    M_0 = data.M_0,
    e = data.e,
    ω = data.ω,
    integrity_status_flag = data.integrity_status_flag,
    l2c_phasing = data.l2c_phasing,
    t_oe = data.t_oe,
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
    t_oc = data.t_oc,
    a_f0 = data.a_f0,
    a_f1 = data.a_f1,
    a_f2 = data.a_f2,
    URA_NED0 = data.URA_NED0,
    URA_NED1 = data.URA_NED1,
    URA_NED2 = data.URA_NED2,
    T_GD = data.T_GD,
    ISC_L1CA = data.ISC_L1CA,
    ISC_L2C = data.ISC_L2C,
    ISC_L5I5 = data.ISC_L5I5,
    ISC_L5Q5 = data.ISC_L5Q5,
    α_0 = data.α_0,
    α_1 = data.α_1,
    α_2 = data.α_2,
    α_3 = data.α_3,
    β_0 = data.β_0,
    β_1 = data.β_1,
    β_2 = data.β_2,
    β_3 = data.β_3,
    WN_op = data.WN_op,
    t_EOP = data.t_EOP,
    PM_X = data.PM_X,
    PM_X_dot = data.PM_X_dot,
    PM_Y = data.PM_Y,
    PM_Y_dot = data.PM_Y_dot,
    ΔUT_GPS = data.ΔUT_GPS,
    ΔUT_GPS_dot = data.ΔUT_GPS_dot,
    A_0_utc = data.A_0_utc,
    A_1_utc = data.A_1_utc,
    A_2_utc = data.A_2_utc,
    Δt_LS = data.Δt_LS,
    t_ot = data.t_ot,
    WN_ot = data.WN_ot,
    WN_LSF = data.WN_LSF,
    DN = data.DN,
    Δt_LSF = data.Δt_LSF,
    t_GGTO = data.t_GGTO,
    WN_GGTO = data.WN_GGTO,
    GNSS_id = data.GNSS_id,
    A_0_GGTO = data.A_0_GGTO,
    A_1_GGTO = data.A_1_GGTO,
    A_2_GGTO = data.A_2_GGTO,
    WN_a = data.WN_a,
    t_oa = data.t_oa,
    reduced_almanacs = data.reduced_almanacs,
    midi_almanacs = data.midi_almanacs,
    t_op_D = data.t_op_D,
    t_OD = data.t_OD,
    cdc_packets = data.cdc_packets,
    edc_packets = data.edc_packets,
    text15_bits = data.text15_bits,
    text15_page = data.text15_page,
    text36_bits = data.text36_bits,
    text36_page = data.text36_page,
    ism = data.ism,
)
    GPSL5Data(
        last_message_id, TOW, alert_flag,
        WN, signal_health_l1, signal_health_l2, signal_health_l5,
        t_op, URA_ED, ΔA, A_dot, Δn_0, Δn_0_dot, M_0, e, ω,
        integrity_status_flag, l2c_phasing,
        t_oe, Ω_0, i_0, Ω_dot, i_dot, C_is, C_ic, C_rs, C_rc, C_us, C_uc,
        t_oc, a_f0, a_f1, a_f2, URA_NED0, URA_NED1, URA_NED2,
        T_GD, ISC_L1CA, ISC_L2C, ISC_L5I5, ISC_L5Q5,
        α_0, α_1, α_2, α_3, β_0, β_1, β_2, β_3, WN_op,
        t_EOP, PM_X, PM_X_dot, PM_Y, PM_Y_dot, ΔUT_GPS, ΔUT_GPS_dot,
        A_0_utc, A_1_utc, A_2_utc, Δt_LS, t_ot, WN_ot, WN_LSF, DN, Δt_LSF,
        t_GGTO, WN_GGTO, GNSS_id, A_0_GGTO, A_1_GGTO, A_2_GGTO,
        WN_a, t_oa, reduced_almanacs, midi_almanacs,
        t_op_D, t_OD, cdc_packets, edc_packets,
        text15_bits, text15_page, text36_bits, text36_page,
        ism,
    )
end

# Field-by-field equality so Vector and struct fields with mutable element types compare correctly.
function Base.:(==)(a::GPSL5Data, b::GPSL5Data)
    for f in fieldnames(GPSL5Data)
        getfield(a, f) == getfield(b, f) || return false
    end
    return true
end

# =============================================================================
# Streaming Viterbi decoder for L5 CNAV (K=7, rate 1/2, G1=171, G2=133 octal)
# =============================================================================
#
# CNAV FEC is convolved continuously across message boundaries (IS-GPS-705J
# §3.3.3.1.1) — there are no tail bits and no encoder reset between messages.
# A receiver that block-decodes individual messages would lose ~5 bits of
# accuracy at each message boundary because the trellis cannot anchor to a
# known final state.
#
# This decoder maintains 64 path metrics and a 32-bit decision history per
# state. Each symbol pair triggers an add-compare-select (ACS) update over all
# 64 states, after which one decoded bit is emitted by chainback from the best
# state through 32 prior decisions. Persistent state lets us decode a
# continuous symbol stream with no boundary loss.
#
# Symbol-pair phase (which received symbol is the "first" of a pair) is not
# known at acquisition time, so the receiver runs two Viterbi instances
# offset by one symbol; whichever instance achieves a successful preamble +
# CRC match acquires "message lock" and the other is discarded. This mirrors
# the GNSS-SDR libswiftcnav design.

const GPSL5_VITERBI_TRACEBACK = 32  # bits; ~5 × constraint length
const GPSL5_VITERBI_NUM_STATES = 64  # 2^(K-1) for K=7

# Precompute the encoded output (G1, G2) for every (state, input_bit) pair.
# State numbering: bits [b1 b2 b3 b4 b5 b6] where b1 is the most-recently-shifted-in.
# Input bit u feeds into stage 0; the encoder output bits are
#   y1 = u ⊕ s3 ⊕ s4 ⊕ s5 ⊕ s6   (G1 = 171 octal = 0b1111001)
#   y2 = u ⊕ s1 ⊕ s3 ⊕ s4 ⊕ s6   (G2 = 133 octal = 0b1011011)
# where s1..s6 are the six register stages (state bits, MSB first).
# Reference: IS-GPS-705J Figure 3-7.
function _gpsl5_viterbi_output(state::UInt8, u::UInt8)
    # state is a 6-bit value: bit 5 = s1 (oldest), bit 0 = s6 (newest)
    s1 = (state >> 5) & 0x01
    s3 = (state >> 3) & 0x01
    s4 = (state >> 2) & 0x01
    s5 = (state >> 1) & 0x01
    s6 = state & 0x01
    # G1 = 1 + x³ + x⁴ + x⁵ + x⁶  (octal 171 = 0b1111001) → taps u, s3, s4, s5, s6
    y1 = u ⊻ s3 ⊻ s4 ⊻ s5 ⊻ s6
    # G2 = 1 + x + x³ + x⁴ + x⁶   (octal 133 = 0b1011011) → taps u, s1, s3, s4, s6
    y2 = u ⊻ s1 ⊻ s3 ⊻ s4 ⊻ s6
    return (y1, y2)
end

# Lookup tables: for each (state, u), the next state and the (y1, y2) pair.
# Size: 64 × 2. Stored as flat tuples for type stability.
const _GPSL5_VITERBI_NEXT_STATE = let
    arr = Array{UInt8}(undef, GPSL5_VITERBI_NUM_STATES, 2)
    for s = 0:GPSL5_VITERBI_NUM_STATES-1, u = 0:1
        # New state shifts u into s6 (LSB) and drops the oldest (s1, MSB).
        ns = (UInt8(u) << 5) | (UInt8(s) >> 1)
        arr[s+1, u+1] = ns & 0x3f
    end
    arr
end

const _GPSL5_VITERBI_OUTPUT_Y1 = let
    arr = Array{UInt8}(undef, GPSL5_VITERBI_NUM_STATES, 2)
    for s = 0:GPSL5_VITERBI_NUM_STATES-1, u = 0:1
        y1, _ = _gpsl5_viterbi_output(UInt8(s), UInt8(u))
        arr[s+1, u+1] = y1
    end
    arr
end

const _GPSL5_VITERBI_OUTPUT_Y2 = let
    arr = Array{UInt8}(undef, GPSL5_VITERBI_NUM_STATES, 2)
    for s = 0:GPSL5_VITERBI_NUM_STATES-1, u = 0:1
        _, y2 = _gpsl5_viterbi_output(UInt8(s), UInt8(u))
        arr[s+1, u+1] = y2
    end
    arr
end

# Per-state predecessor lookup — for ACS we need, for each *new* state s', the
# two predecessor (state, u) pairs whose transitions lead into s'.
# A new state s' = (u << 5) | (s >> 1) means the predecessor state s has its
# top 5 bits equal to s' & 0x1F shifted left by one (with low bit either 0 or 1),
# and the input u is the top bit of s'.
const _GPSL5_VITERBI_PRED = let
    # arr[s'+1, k] gives (predecessor_state, input_bit) for k = 1, 2.
    arr_state = Array{UInt8}(undef, GPSL5_VITERBI_NUM_STATES, 2)
    arr_input = Array{UInt8}(undef, GPSL5_VITERBI_NUM_STATES, 2)
    for sp = 0:GPSL5_VITERBI_NUM_STATES-1
        u = UInt8((sp >> 5) & 0x01)
        # Predecessor's top 5 bits are sp's bottom 5 bits shifted up.
        base = UInt8((sp & 0x1f) << 1)
        arr_state[sp+1, 1] = base | 0x00
        arr_state[sp+1, 2] = base | 0x01
        arr_input[sp+1, 1] = u
        arr_input[sp+1, 2] = u
    end
    (arr_state, arr_input)
end

# Streaming Viterbi state. Designed to be functionally updated (no in-place
# mutation) so that GNSSDecoderState semantics carry over. Path metrics are
# stored as `Vector{UInt32}` (length 64) and decision histories as
# `Vector{UInt32}` (length 64), each holding the most recent 32 decisions per
# state.
struct GPSL5ViterbiState
    metrics::Vector{UInt32}
    decisions::Vector{UInt32}
    n_steps::Int  # how many ACS steps have happened so far (capped at traceback)
end

GPSL5ViterbiState() = GPSL5ViterbiState(
    # Start from an unknown initial state: bias state 0 to be "best" (metric 0),
    # all others penalized. After ~5×K = 35 steps the trellis converges
    # regardless of which start we pick.
    UInt32[i == 1 ? 0x00000000 : 0x000000ff for i = 1:GPSL5_VITERBI_NUM_STATES],
    zeros(UInt32, GPSL5_VITERBI_NUM_STATES),
    0,
)

# Single ACS step. Inputs are two received symbols (each 0 or 1 for hard
# decision). Returns updated state and an Int8: -1 if no decoded bit is
# available yet, otherwise the decoded bit (0 or 1).
#
# Implementation note: this is "register-exchange" Viterbi — each state carries
# the recent input-bit history that produced its current path metric. On each
# ACS, the surviving predecessor's input-bit history is shifted left by one
# and the *input bit* `u` (the top bit of the new state) is appended in the
# LSB. Once `n_steps >= TRACEBACK`, the bit at position (TRACEBACK-1) of the
# best state's history is the decoded bit at depth TRACEBACK.
function _gpsl5_viterbi_step(v::GPSL5ViterbiState, sym1::UInt8, sym2::UInt8)
    new_metrics = Vector{UInt32}(undef, GPSL5_VITERBI_NUM_STATES)
    new_decisions = Vector{UInt32}(undef, GPSL5_VITERBI_NUM_STATES)
    pred_state, _ = _GPSL5_VITERBI_PRED
    @inbounds for sp = 0:GPSL5_VITERBI_NUM_STATES-1
        ps0 = pred_state[sp+1, 1]
        ps1 = pred_state[sp+1, 2]
        # Both predecessors of `sp` arrive via the same input bit u (top bit of sp).
        u = UInt8((sp >> 5) & 0x01)
        y1_0 = _GPSL5_VITERBI_OUTPUT_Y1[ps0+1, u+1]
        y2_0 = _GPSL5_VITERBI_OUTPUT_Y2[ps0+1, u+1]
        y1_1 = _GPSL5_VITERBI_OUTPUT_Y1[ps1+1, u+1]
        y2_1 = _GPSL5_VITERBI_OUTPUT_Y2[ps1+1, u+1]
        # Hamming distance branch metric (hard decision).
        bm0 = UInt32((y1_0 ⊻ sym1) + (y2_0 ⊻ sym2))
        bm1 = UInt32((y1_1 ⊻ sym1) + (y2_1 ⊻ sym2))
        m0 = v.metrics[ps0+1] + bm0
        m1 = v.metrics[ps1+1] + bm1
        if m0 <= m1
            new_metrics[sp+1] = m0
            # Append the input bit `u` to the surviving predecessor's history.
            new_decisions[sp+1] = (v.decisions[ps0+1] << 1) | UInt32(u)
        else
            new_metrics[sp+1] = m1
            new_decisions[sp+1] = (v.decisions[ps1+1] << 1) | UInt32(u)
        end
    end

    # Renormalize metrics: subtract the minimum to keep numbers small. With
    # 32-bit metrics this is just hygiene since branch metrics are 0..2 per step.
    min_m = minimum(new_metrics)
    if min_m > 0
        @inbounds for i = 1:GPSL5_VITERBI_NUM_STATES
            new_metrics[i] -= min_m
        end
    end

    n_steps = min(v.n_steps + 1, GPSL5_VITERBI_TRACEBACK)
    decoded_bit::Int8 = -1
    if v.n_steps + 1 >= GPSL5_VITERBI_TRACEBACK
        # Find the state with smallest path metric.
        best_state = 0
        best_metric = new_metrics[1]
        @inbounds for i = 2:GPSL5_VITERBI_NUM_STATES
            if new_metrics[i] < best_metric
                best_metric = new_metrics[i]
                best_state = i - 1
            end
        end
        # Bit at depth (TRACEBACK-1) — counted from the most recently appended
        # bit (LSB) — is the decoded bit at time t-(TRACEBACK-1).
        decoded_bit =
            Int8((new_decisions[best_state+1] >> (GPSL5_VITERBI_TRACEBACK - 1)) & UInt32(1))
    end

    return GPSL5ViterbiState(new_metrics, new_decisions, n_steps), decoded_bit
end

# =============================================================================
# Per-Viterbi-instance framing state (preamble sync, message buffer, lock)
# =============================================================================

mutable struct _GPSL5UnusedSentinel end  # placeholder for type discrimination

"""
Per-Viterbi-instance framing state. Two of these run in parallel
(symbol-pair-phase 0 and 1) until one acquires message lock by matching the
preamble and CRC.
"""
struct GPSL5FramingPart
    viterbi::GPSL5ViterbiState
    pending_symbol::Int8     # -1 if no pending symbol; otherwise 0 or 1
    decoded::UInt320         # bit buffer (MSB = oldest decoded bit)
    n_decoded::Int           # number of bits currently in `decoded`
    preamble_seen::Bool
    invert::Bool             # true if the preamble matched its inverse
    message_lock::Bool
    n_crc_fail::Int
end

GPSL5FramingPart() =
    GPSL5FramingPart(GPSL5ViterbiState(), Int8(-1), UInt320(0), 0, false, false, false, 0)

function GPSL5FramingPart(
    p::GPSL5FramingPart;
    viterbi = p.viterbi,
    pending_symbol = p.pending_symbol,
    decoded = p.decoded,
    n_decoded = p.n_decoded,
    preamble_seen = p.preamble_seen,
    invert = p.invert,
    message_lock = p.message_lock,
    n_crc_fail = p.n_crc_fail,
)
    GPSL5FramingPart(
        viterbi, pending_symbol, decoded, n_decoded,
        preamble_seen, invert, message_lock, n_crc_fail,
    )
end

"""
    GPSL5Cache

Decoder cache state for GPS L5 CNAV. Holds two parallel Viterbi/framing
instances to resolve symbol-pair phase ambiguity, and a counter for symbols
ingested before either instance acquires message lock.
"""
struct GPSL5Cache <: AbstractGNSSCache
    part1::GPSL5FramingPart  # symbol-pair-phase 0 (default)
    part2::GPSL5FramingPart  # symbol-pair-phase 1 (offset by one symbol)
    locked_part::Int         # 0 = neither, 1 = part1, 2 = part2
end

GPSL5Cache() = GPSL5Cache(GPSL5FramingPart(), GPSL5FramingPart(), 0)

function GPSL5Cache(
    c::GPSL5Cache;
    part1 = c.part1,
    part2 = c.part2,
    locked_part = c.locked_part,
)
    GPSL5Cache(part1, part2, locked_part)
end

# `Vector` fields (path metrics, decisions) compare by reference under the
# default `==`, so we override these with field-by-field equality. This
# matches the pattern used by `GalileoE1BData` for its `almanacs::Vector` field.
function Base.:(==)(a::GPSL5ViterbiState, b::GPSL5ViterbiState)
    a.metrics == b.metrics && a.decisions == b.decisions && a.n_steps == b.n_steps
end

function Base.:(==)(a::GPSL5FramingPart, b::GPSL5FramingPart)
    for f in fieldnames(GPSL5FramingPart)
        getfield(a, f) == getfield(b, f) || return false
    end
    return true
end

function Base.:(==)(a::GPSL5Cache, b::GPSL5Cache)
    a.part1 == b.part1 && a.part2 == b.part2 && a.locked_part == b.locked_part
end

# =============================================================================
# Constructor / reset
# =============================================================================

"""
$(TYPEDSIGNATURES)

Create a decoder state for GPS L5 CNAV navigation messages.

Initializes a [`GNSSDecoderState`](@ref) configured for decoding GPS L5
(L5-I) civil navigation (CNAV) messages. The decoder accepts FEC-encoded
100 sps symbols, runs an internal continuous Viterbi decoder, locates the
preamble, validates CRC-24Q, and parses message types 10, 11, 12, 13, 14,
15, 30-37, and 40 as defined in IS-GPS-705J.

# Arguments
- `prn::Int`: Pseudo-Random Noise code identifier (1-63 for GPS satellites)

# Returns
- `GNSSDecoderState{GPSL5Data}`: Initialized decoder state for GPS L5

# Example
```julia
state = GPSL5DecoderState(1)  # Create decoder for PRN 1
state = decode(state, symbols, num_symbols)  # symbols are 100 sps FEC symbols
if is_sat_healthy(state)
    # Use state.data for positioning
end
```

# See Also
- [`GNSSDecoderState`](@ref): The underlying state structure
- [`decode`](@ref): Decode symbols using this state
- [`reset_decoder_state`](@ref): Reset after signal loss
- [`is_sat_healthy`](@ref): Check satellite health status
"""
function GPSL5DecoderState(prn)
    GNSSDecoderState(
        prn,
        UInt320(0),
        UInt320(0),
        GPSL5Data(),
        GPSL5Data(),
        GPSL5Constants(),
        GPSL5Cache(),
        0,
        nothing,
        false,
    )
end

function GNSSDecoderState(::GPSL5, prn)
    GPSL5DecoderState(prn)
end

"""
$(TYPEDSIGNATURES)

Reset the GPS L5 decoder state after a signal loss or reacquisition.

Clears Viterbi state, framing buffers, and the time-of-week (TOW) field
while preserving the remaining decoded ephemeris and clock data in `raw_data`.
This allows the decoder to be re-used after a brief signal outage without
discarding all previously collected navigation data.

# Arguments
- `state::GNSSDecoderState{<:GPSL5Data}`: Current GPS L5 decoder state

# Returns
- `GNSSDecoderState{<:GPSL5Data}`: Reset decoder state with cleared Viterbi/framing state

# See Also
- [`GPSL5DecoderState`](@ref): Create a fresh decoder state
- [`decode`](@ref): Continue decoding after reset
"""
function reset_decoder_state(state::GNSSDecoderState{<:GPSL5Data})
    GNSSDecoderState(
        state;
        raw_buffer = UInt320(0),
        buffer = UInt320(0),
        raw_data = GPSL5Data(state.raw_data; TOW = nothing),
        data = GPSL5Data(),
        cache = GPSL5Cache(),
        num_bits_buffered = 0,
        num_bits_after_valid_syncro_sequence = nothing,
    )
end

# =============================================================================
# Preamble + CRC framing on the per-part decoded bit buffer
# =============================================================================

# CRC-24Q for CNAV — same generator polynomial as Galileo (0x864CFB), seed 0,
# applied to the first 276 information bits (preamble + payload, excluding
# the 24 CRC bits). Reuses the existing `galCRC24` instance defined in
# GNSSDecoder.jl.
#
# To check a 300-bit message we feed all 300 bits and verify the output is 0
# (a property of CRC-24Q: appending the checksum makes the codeword
# CRC-clean). The `galCRC24` instance was constructed for byte-aligned input,
# so we right-align the 300 bits in 38 bytes (304 bits, 4 leading zeros) and
# evaluate the CRC over those bytes.
function _gpsl5_check_crc(message_bits::UInt320)
    # message_bits holds 300 bits in its low 300 positions. Pad with 4 leading
    # zero bits to get 304 bits = 38 bytes, then check CRC over those.
    padded = UInt320(message_bits) << 4
    bytes = reverse(digits(UInt8, padded; base = 256, pad = 38))
    return galCRC24(bytes) == 0
end

# Append a decoded bit to the part's buffer, pushing oldest bits up (MSB).
# Buffer holds at most 304 bits (300 message + extra for preamble lookahead).
function _gpsl5_push_bit(part::GPSL5FramingPart, bit::Int8)
    # bit is 0 or 1.
    new_decoded = (part.decoded << 1) | UInt320(bit)
    new_n = min(part.n_decoded + 1, 308)  # cap so we don't overflow UInt320 logically
    return GPSL5FramingPart(part; decoded = new_decoded, n_decoded = new_n)
end

# Try to find a preamble (0x8B normal or 0x74 inverted) in the buffered bits.
# When found, shift the buffer so the preamble starts at bit position
# (n_decoded - GPSL5_PREAMBLE_LENGTH).
# Mirrors GNSS-SDR's cnav_rescan_preamble_.
function _gpsl5_rescan_preamble(part::GPSL5FramingPart)
    n = part.n_decoded
    if n < 8
        return GPSL5FramingPart(part; preamble_seen = false, invert = false)
    end
    # Search positions: the preamble must lie at the *oldest* bits if the
    # buffer is exactly 8 bits wide. As more bits accumulate, the oldest 8
    # bits represent the candidate preamble. We try every alignment
    # (1..n-7) and keep the buffer trimmed to where the preamble starts.
    for offset = 0:(n - 8)
        # Extract 8 bits starting `offset` positions from the oldest bit
        # (i.e. shift right by (n - 8 - offset) and mask).
        candidate = UInt8((part.decoded >> UInt(n - 8 - offset)) & UInt320(0xff))
        if candidate == 0x8b
            # Drop the `offset` oldest bits (mask them out).
            mask = (UInt320(1) << UInt(n - offset)) - UInt320(1)
            new_decoded = part.decoded & mask
            return GPSL5FramingPart(
                part;
                decoded = new_decoded,
                n_decoded = n - offset,
                preamble_seen = true,
                invert = false,
            )
        elseif candidate == 0x74
            mask = (UInt320(1) << UInt(n - offset)) - UInt320(1)
            new_decoded = part.decoded & mask
            return GPSL5FramingPart(
                part;
                decoded = new_decoded,
                n_decoded = n - offset,
                preamble_seen = true,
                invert = true,
            )
        end
    end
    # Not found — keep only the most recent 7 bits so that the next pushed
    # bit can complete a candidate preamble at the oldest end.
    keep = min(n, 7)
    mask = (UInt320(1) << UInt(keep)) - UInt320(1)
    new_decoded = part.decoded & mask
    return GPSL5FramingPart(
        part;
        decoded = new_decoded,
        n_decoded = keep,
        preamble_seen = false,
        invert = false,
    )
end

# Extract the 300 oldest bits as a complete message. Caller must check that
# n_decoded >= 300.
function _gpsl5_extract_message(part::GPSL5FramingPart)
    # The 300 oldest bits live at positions (n_decoded-1 .. n_decoded-300).
    msg_bits = (part.decoded >> UInt(part.n_decoded - 300)) &
        ((UInt320(1) << 300) - UInt320(1))
    if part.invert
        msg_bits = (~msg_bits) & ((UInt320(1) << 300) - UInt320(1))
    end
    # Drop the consumed 300 oldest bits from the buffer.
    keep = part.n_decoded - 300
    if keep > 0
        mask = (UInt320(1) << UInt(keep)) - UInt320(1)
        new_decoded = part.decoded & mask
    else
        new_decoded = UInt320(0)
    end
    new_part = GPSL5FramingPart(part; decoded = new_decoded, n_decoded = keep)
    return new_part, msg_bits
end

const GPSL5_LOCK_MAX_CRC_FAILS = 10

# Drop the oldest bit from a part's buffer (used after a spurious-preamble
# CRC failure so the next rescan considers a different alignment).
function _gpsl5_drop_oldest_bit(part::GPSL5FramingPart)
    n = part.n_decoded
    if n == 0
        return part
    end
    keep = n - 1
    mask = keep > 0 ? (UInt320(1) << UInt(keep)) - UInt320(1) : UInt320(0)
    return GPSL5FramingPart(part; decoded = part.decoded & mask, n_decoded = keep)
end

# After a CRC failure on a locked stream, advance one bit and try to re-sync
# (drop the lock).
function _gpsl5_drop_lock(part::GPSL5FramingPart)
    GPSL5FramingPart(
        part;
        preamble_seen = false,
        invert = false,
        message_lock = false,
        n_crc_fail = 0,
    )
end

# =============================================================================
# Bit-field access on a 300-bit CNAV message
# =============================================================================

# Extract `length` bits from a 300-bit message starting at 1-indexed position `start`.
# Bit 1 is the MSB of the broadcast (the leftmost bit of the preamble).
function _gpsl5_get_bits(msg::UInt320, start::Int, length::Int)
    return UInt(get_bits(msg, 300, start, length))
end

function _gpsl5_get_bit(msg::UInt320, position::Int)
    return Bool(_gpsl5_get_bits(msg, position, 1))
end

function _gpsl5_get_signed(msg::UInt320, start::Int, length::Int)
    return get_twos_complement_num(msg, 300, start, length)
end

# =============================================================================
# Per-message-type decoders
# =============================================================================

# Decode the common 38-bit header (preamble, PRN, MT ID, TOW count, alert flag).
function _gpsl5_decode_common_header(state::GNSSDecoderState{<:GPSL5Data}, msg::UInt320)
    msg_id = Int(_gpsl5_get_bits(msg, 15, 6))
    tow_count = Int64(_gpsl5_get_bits(msg, 21, 17))
    TOW = tow_count * 6  # message TOW is 17 MSBs of actual TOW at start of *next* 6-second message
    alert_flag = _gpsl5_get_bit(msg, 38)
    return GNSSDecoderState(
        state;
        raw_data = GPSL5Data(
            state.raw_data;
            last_message_id = msg_id,
            TOW = TOW,
            alert_flag = alert_flag,
        ),
    )
end

function _gpsl5_decode_mt10(state::GNSSDecoderState{<:GPSL5Data}, msg::UInt320)
    constants = state.constants
    WN = Int64(_gpsl5_get_bits(msg, 39, 13))
    signal_health_l1 = _gpsl5_get_bit(msg, 52)
    signal_health_l2 = _gpsl5_get_bit(msg, 53)
    signal_health_l5 = _gpsl5_get_bit(msg, 54)
    t_op = Int64(_gpsl5_get_bits(msg, 55, 11)) * 300
    URA_ED = Int(_gpsl5_get_signed(msg, 66, 5))
    t_oe = Int64(_gpsl5_get_bits(msg, 71, 11)) * 300
    ΔA = _gpsl5_get_signed(msg, 82, 26) / Float64(1 << 9)
    A_dot = _gpsl5_get_signed(msg, 108, 25) / Float64(1 << 21)
    Δn_0 = _gpsl5_get_signed(msg, 133, 17) * constants.PI / Float64(1 << 44)
    Δn_0_dot = _gpsl5_get_signed(msg, 150, 23) * constants.PI / Float64(1 << 57)
    M_0 = _gpsl5_get_signed_60(msg, 173, 33) * constants.PI / Float64(1 << 32)
    e = _gpsl5_get_bits_60(msg, 206, 33) / Float64(1 << 34)
    ω = _gpsl5_get_signed_60(msg, 239, 33) * constants.PI / Float64(1 << 32)
    integrity_status_flag = _gpsl5_get_bit(msg, 272)
    l2c_phasing = _gpsl5_get_bit(msg, 273)
    return GNSSDecoderState(
        state;
        raw_data = GPSL5Data(
            state.raw_data;
            WN, signal_health_l1, signal_health_l2, signal_health_l5,
            t_op, URA_ED, t_oe, ΔA, A_dot, Δn_0, Δn_0_dot, M_0, e, ω,
            integrity_status_flag, l2c_phasing,
        ),
    )
end

# Wide-bit accessors for fields longer than 32 bits (33 bits for M_0, e, ω).
# get_bits returns the field zero-extended into UInt; for ≥33-bit fields we
# need to widen first.
function _gpsl5_get_bits_60(msg::UInt320, start::Int, length::Int)
    return UInt64((msg >> UInt(300 - start - length + 1)) &
        ((UInt320(1) << length) - UInt320(1)))
end

function _gpsl5_get_signed_60(msg::UInt320, start::Int, length::Int)
    value = _gpsl5_get_bits_60(msg, start, length)
    shift = 64 - length
    return signed(value << shift) >> shift
end

function _gpsl5_decode_mt11(state::GNSSDecoderState{<:GPSL5Data}, msg::UInt320)
    constants = state.constants
    t_oe = Int64(_gpsl5_get_bits(msg, 39, 11)) * 300
    Ω_0 = _gpsl5_get_signed_60(msg, 50, 33) * constants.PI / Float64(1 << 32)
    i_0 = _gpsl5_get_signed_60(msg, 83, 33) * constants.PI / Float64(1 << 32)
    Ω_dot = _gpsl5_get_signed(msg, 116, 17) * constants.PI / Float64(1 << 44)
    i_dot = _gpsl5_get_signed(msg, 133, 15) * constants.PI / Float64(1 << 44)
    C_is = _gpsl5_get_signed(msg, 148, 16) / Float64(1 << 30)
    C_ic = _gpsl5_get_signed(msg, 164, 16) / Float64(1 << 30)
    C_rs = _gpsl5_get_signed(msg, 180, 24) / Float64(1 << 8)
    C_rc = _gpsl5_get_signed(msg, 204, 24) / Float64(1 << 8)
    C_us = _gpsl5_get_signed(msg, 228, 21) / Float64(1 << 30)
    C_uc = _gpsl5_get_signed(msg, 249, 21) / Float64(1 << 30)
    return GNSSDecoderState(
        state;
        raw_data = GPSL5Data(
            state.raw_data;
            t_oe, Ω_0, i_0, Ω_dot, i_dot,
            C_is, C_ic, C_rs, C_rc, C_us, C_uc,
        ),
    )
end

# MT 30-37 share a common clock-correction header at bits 39..127.
function _gpsl5_decode_mt30_37_clock(state::GNSSDecoderState{<:GPSL5Data}, msg::UInt320)
    t_op = Int64(_gpsl5_get_bits(msg, 39, 11)) * 300
    URA_NED0 = Int(_gpsl5_get_signed(msg, 50, 5))
    URA_NED1 = Int(_gpsl5_get_bits(msg, 55, 3))
    URA_NED2 = Int(_gpsl5_get_bits(msg, 58, 3))
    t_oc = Int64(_gpsl5_get_bits(msg, 61, 11)) * 300
    a_f0 = _gpsl5_get_signed(msg, 72, 26) / Float64(1 << 35)
    a_f1 = _gpsl5_get_signed(msg, 98, 20) / Float64(1 << 48)
    a_f2 = _gpsl5_get_signed(msg, 118, 10) / Float64(1 << 60)
    return GNSSDecoderState(
        state;
        raw_data = GPSL5Data(
            state.raw_data;
            t_op, URA_NED0, URA_NED1, URA_NED2, t_oc, a_f0, a_f1, a_f2,
        ),
    )
end

function _gpsl5_decode_mt30(state::GNSSDecoderState{<:GPSL5Data}, msg::UInt320)
    state = _gpsl5_decode_mt30_37_clock(state, msg)
    T_GD = _gpsl5_get_signed(msg, 128, 13) / Float64(1 << 35)
    ISC_L1CA = _gpsl5_get_signed(msg, 141, 13) / Float64(1 << 35)
    ISC_L2C = _gpsl5_get_signed(msg, 154, 13) / Float64(1 << 35)
    ISC_L5I5 = _gpsl5_get_signed(msg, 167, 13) / Float64(1 << 35)
    ISC_L5Q5 = _gpsl5_get_signed(msg, 180, 13) / Float64(1 << 35)
    α_0 = _gpsl5_get_signed(msg, 193, 8) / Float64(1 << 30)
    α_1 = _gpsl5_get_signed(msg, 201, 8) / Float64(1 << 27)
    α_2 = _gpsl5_get_signed(msg, 209, 8) / Float64(1 << 24)
    α_3 = _gpsl5_get_signed(msg, 217, 8) / Float64(1 << 24)
    β_0 = _gpsl5_get_signed(msg, 225, 8) * Float64(1 << 11)
    β_1 = _gpsl5_get_signed(msg, 233, 8) * Float64(1 << 14)
    β_2 = _gpsl5_get_signed(msg, 241, 8) * Float64(1 << 16)
    β_3 = _gpsl5_get_signed(msg, 249, 8) * Float64(1 << 16)
    WN_op = UInt(_gpsl5_get_bits(msg, 257, 8))
    return GNSSDecoderState(
        state;
        raw_data = GPSL5Data(
            state.raw_data;
            T_GD, ISC_L1CA, ISC_L2C, ISC_L5I5, ISC_L5Q5,
            α_0, α_1, α_2, α_3, β_0, β_1, β_2, β_3, WN_op,
        ),
    )
end

function _gpsl5_decode_mt31(state::GNSSDecoderState{<:GPSL5Data}, msg::UInt320)
    state = _gpsl5_decode_mt30_37_clock(state, msg)
    WN_a = UInt(_gpsl5_get_bits(msg, 128, 13))
    t_oa = Int(_gpsl5_get_bits(msg, 141, 8)) * 4096
    # Four reduced-almanac packets at bits 149, 180, 211, 242 (each 31 bits).
    new_almanacs = copy(state.raw_data.reduced_almanacs)
    for start in (149, 180, 211, 242)
        pkt = _gpsl5_decode_reduced_almanac_packet(msg, start)
        if !isnothing(pkt.PRN_a) && pkt.PRN_a > 0 && pkt.PRN_a <= length(new_almanacs)
            new_almanacs[pkt.PRN_a] = pkt
        end
    end
    return GNSSDecoderState(
        state;
        raw_data = GPSL5Data(
            state.raw_data;
            WN_a, t_oa, reduced_almanacs = new_almanacs,
        ),
    )
end

function _gpsl5_decode_reduced_almanac_packet(msg::UInt320, start::Int)
    PRN_a = UInt(_gpsl5_get_bits(msg, start, 6))
    δA = _gpsl5_get_signed(msg, start + 6, 8) * 512.0
    Ω_0 = _gpsl5_get_signed(msg, start + 14, 7) / Float64(1 << 6)
    Φ_0 = _gpsl5_get_signed(msg, start + 21, 7) / Float64(1 << 6)
    sh_l1 = _gpsl5_get_bit(msg, start + 28)
    sh_l2 = _gpsl5_get_bit(msg, start + 29)
    sh_l5 = _gpsl5_get_bit(msg, start + 30)
    return GPSL5ReducedAlmanacPacket(;
        PRN_a, δA, Ω_0, Φ_0,
        signal_health_l1 = sh_l1,
        signal_health_l2 = sh_l2,
        signal_health_l5 = sh_l5,
    )
end

function _gpsl5_decode_mt32(state::GNSSDecoderState{<:GPSL5Data}, msg::UInt320)
    state = _gpsl5_decode_mt30_37_clock(state, msg)
    t_EOP = Int64(_gpsl5_get_bits(msg, 128, 16)) * 16
    PM_X = _gpsl5_get_signed(msg, 144, 21) / Float64(1 << 20)
    PM_X_dot = _gpsl5_get_signed(msg, 165, 15) / Float64(1 << 21)
    PM_Y = _gpsl5_get_signed(msg, 180, 21) / Float64(1 << 20)
    PM_Y_dot = _gpsl5_get_signed(msg, 201, 15) / Float64(1 << 21)
    ΔUT_GPS = _gpsl5_get_signed_60(msg, 216, 31) / Float64(1 << 23)
    ΔUT_GPS_dot = _gpsl5_get_signed(msg, 247, 19) / Float64(1 << 25)
    return GNSSDecoderState(
        state;
        raw_data = GPSL5Data(
            state.raw_data;
            t_EOP, PM_X, PM_X_dot, PM_Y, PM_Y_dot, ΔUT_GPS, ΔUT_GPS_dot,
        ),
    )
end

function _gpsl5_decode_mt33(state::GNSSDecoderState{<:GPSL5Data}, msg::UInt320)
    state = _gpsl5_decode_mt30_37_clock(state, msg)
    A_0_utc = _gpsl5_get_signed(msg, 128, 16) / Float64(1 << 35)
    A_1_utc = _gpsl5_get_signed(msg, 144, 13) / Float64(1 << 51)
    A_2_utc = _gpsl5_get_signed(msg, 157, 7) / Float64(1 << 68)
    Δt_LS = Int(_gpsl5_get_signed(msg, 164, 8))
    t_ot = Int(_gpsl5_get_bits(msg, 172, 16)) * 16
    WN_ot = UInt(_gpsl5_get_bits(msg, 188, 13))
    WN_LSF = UInt(_gpsl5_get_bits(msg, 201, 13))
    DN = UInt(_gpsl5_get_bits(msg, 214, 4))
    Δt_LSF = Int(_gpsl5_get_signed(msg, 218, 8))
    return GNSSDecoderState(
        state;
        raw_data = GPSL5Data(
            state.raw_data;
            A_0_utc, A_1_utc, A_2_utc, Δt_LS, t_ot, WN_ot, WN_LSF, DN, Δt_LSF,
        ),
    )
end

function _gpsl5_decode_mt35(state::GNSSDecoderState{<:GPSL5Data}, msg::UInt320)
    state = _gpsl5_decode_mt30_37_clock(state, msg)
    t_GGTO = Int(_gpsl5_get_bits(msg, 128, 16)) * 16
    WN_GGTO = UInt(_gpsl5_get_bits(msg, 144, 13))
    GNSS_id = UInt(_gpsl5_get_bits(msg, 157, 3))
    A_0_GGTO = _gpsl5_get_signed(msg, 160, 16) / Float64(1 << 35)
    A_1_GGTO = _gpsl5_get_signed(msg, 176, 13) / Float64(1 << 51)
    A_2_GGTO = _gpsl5_get_signed(msg, 189, 7) / Float64(1 << 68)
    return GNSSDecoderState(
        state;
        raw_data = GPSL5Data(
            state.raw_data;
            t_GGTO, WN_GGTO, GNSS_id, A_0_GGTO, A_1_GGTO, A_2_GGTO,
        ),
    )
end

function _gpsl5_decode_mt37(state::GNSSDecoderState{<:GPSL5Data}, msg::UInt320)
    state = _gpsl5_decode_mt30_37_clock(state, msg)
    WN_a = UInt(_gpsl5_get_bits(msg, 128, 13))
    t_oa = Int(_gpsl5_get_bits(msg, 141, 8)) * 4096
    PRN_a = UInt(_gpsl5_get_bits(msg, 149, 6))
    sh_l1 = _gpsl5_get_bit(msg, 155)
    sh_l2 = _gpsl5_get_bit(msg, 156)
    sh_l5 = _gpsl5_get_bit(msg, 157)
    e = _gpsl5_get_bits(msg, 158, 11) / Float64(1 << 16)
    δi = _gpsl5_get_signed(msg, 169, 11) / Float64(1 << 14)
    Ω_dot = _gpsl5_get_signed(msg, 180, 11) / Float64(1 << 33)
    sqrt_A = _gpsl5_get_bits(msg, 191, 17) / Float64(1 << 4)
    Ω_0 = _gpsl5_get_signed(msg, 208, 16) / Float64(1 << 15)
    ω = _gpsl5_get_signed(msg, 224, 16) / Float64(1 << 15)
    M_0 = _gpsl5_get_signed(msg, 240, 16) / Float64(1 << 15)
    a_f0 = _gpsl5_get_signed(msg, 256, 11) / Float64(1 << 20)
    a_f1 = _gpsl5_get_signed(msg, 267, 10) / Float64(1 << 37)
    almanac = GPSL5MidiAlmanac(;
        PRN_a,
        signal_health_l1 = sh_l1,
        signal_health_l2 = sh_l2,
        signal_health_l5 = sh_l5,
        e, δi, Ω_dot, sqrt_A, Ω_0, ω, M_0, a_f0, a_f1,
    )
    new_midi = copy(state.raw_data.midi_almanacs)
    if PRN_a > 0 && PRN_a <= length(new_midi)
        new_midi[PRN_a] = almanac
    end
    return GNSSDecoderState(
        state;
        raw_data = GPSL5Data(
            state.raw_data;
            WN_a, t_oa, midi_almanacs = new_midi,
        ),
    )
end

function _gpsl5_decode_mt12(state::GNSSDecoderState{<:GPSL5Data}, msg::UInt320)
    WN_a = UInt(_gpsl5_get_bits(msg, 39, 13))
    t_oa = Int(_gpsl5_get_bits(msg, 52, 8)) * 4096
    new_almanacs = copy(state.raw_data.reduced_almanacs)
    # Seven reduced almanac packets at bits 60, 91, 122, 153, 184, 215, 246
    for start in (60, 91, 122, 153, 184, 215, 246)
        pkt = _gpsl5_decode_reduced_almanac_packet(msg, start)
        if !isnothing(pkt.PRN_a) && pkt.PRN_a > 0 && pkt.PRN_a <= length(new_almanacs)
            new_almanacs[pkt.PRN_a] = pkt
        end
    end
    return GNSSDecoderState(
        state;
        raw_data = GPSL5Data(
            state.raw_data;
            WN_a, t_oa, reduced_almanacs = new_almanacs,
        ),
    )
end

function _gpsl5_decode_cdc_packet(msg::UInt320, start::Int)
    PRN_id = UInt(_gpsl5_get_bits(msg, start, 8))
    if PRN_id == 0xff
        return nothing
    end
    δa_f0 = _gpsl5_get_signed(msg, start + 8, 13) / Float64(1 << 35)
    δa_f1 = _gpsl5_get_signed(msg, start + 21, 8) / Float64(1 << 51)
    UDRA = Int(_gpsl5_get_signed(msg, start + 29, 5))
    return GPSL5ClockDifferentialCorrection(; PRN_id, δa_f0, δa_f1, UDRA)
end

function _gpsl5_decode_edc_packet(msg::UInt320, start::Int)
    PRN_id = UInt(_gpsl5_get_bits(msg, start, 8))
    if PRN_id == 0xff
        return nothing
    end
    Δα = _gpsl5_get_signed(msg, start + 8, 14) / Float64(1 << 34)
    Δβ = _gpsl5_get_signed(msg, start + 22, 14) / Float64(1 << 34)
    Δγ = _gpsl5_get_signed(msg, start + 36, 15) / Float64(1 << 32)
    Δi = _gpsl5_get_signed(msg, start + 51, 12) / Float64(1 << 32)
    ΔΩ = _gpsl5_get_signed(msg, start + 63, 12) / Float64(1 << 32)
    ΔA = _gpsl5_get_signed(msg, start + 75, 12) / Float64(1 << 9)
    UDRA_dot = Int(_gpsl5_get_signed(msg, start + 87, 5))
    return GPSL5EphemerisDifferentialCorrection(;
        PRN_id, Δα, Δβ, Δγ, Δi, ΔΩ, ΔA, UDRA_dot,
    )
end

function _gpsl5_decode_mt13(state::GNSSDecoderState{<:GPSL5Data}, msg::UInt320)
    t_op_D = Int64(_gpsl5_get_bits(msg, 39, 11)) * 300
    t_OD = Int64(_gpsl5_get_bits(msg, 50, 11)) * 300
    # Six CDC packets at bits 61, 96, 131, 166, 201, 236 (each 35 bits = 1 DC type bit + 34 CDC).
    new_cdcs = GPSL5ClockDifferentialCorrection[]
    for start in (62, 97, 132, 167, 202, 237)  # skip the leading DC Data Type bit
        pkt = _gpsl5_decode_cdc_packet(msg, start)
        if pkt !== nothing
            push!(new_cdcs, pkt)
        end
    end
    return GNSSDecoderState(
        state;
        raw_data = GPSL5Data(
            state.raw_data;
            t_op_D, t_OD, cdc_packets = new_cdcs,
        ),
    )
end

function _gpsl5_decode_mt14(state::GNSSDecoderState{<:GPSL5Data}, msg::UInt320)
    t_op_D = Int64(_gpsl5_get_bits(msg, 39, 11)) * 300
    t_OD = Int64(_gpsl5_get_bits(msg, 50, 11)) * 300
    # Two EDC packets at bits 62 (DC type at 61) and 155 (DC type at 154).
    new_edcs = GPSL5EphemerisDifferentialCorrection[]
    for start in (62, 155)
        pkt = _gpsl5_decode_edc_packet(msg, start)
        if pkt !== nothing
            push!(new_edcs, pkt)
        end
    end
    return GNSSDecoderState(
        state;
        raw_data = GPSL5Data(
            state.raw_data;
            t_op_D, t_OD, edc_packets = new_edcs,
        ),
    )
end

function _gpsl5_decode_mt34(state::GNSSDecoderState{<:GPSL5Data}, msg::UInt320)
    state = _gpsl5_decode_mt30_37_clock(state, msg)
    # CDC packet at bits 151 (after DC Data Type at bit 150), then EDC packet at bit 185+ (after another DC Data Type bit somewhere).
    cdc = _gpsl5_decode_cdc_packet(msg, 151)
    edc = _gpsl5_decode_edc_packet(msg, 185)
    new_cdcs = cdc === nothing ? GPSL5ClockDifferentialCorrection[] : [cdc]
    new_edcs = edc === nothing ? GPSL5EphemerisDifferentialCorrection[] : [edc]
    return GNSSDecoderState(
        state;
        raw_data = GPSL5Data(
            state.raw_data;
            cdc_packets = new_cdcs, edc_packets = new_edcs,
        ),
    )
end

function _gpsl5_decode_mt15(state::GNSSDecoderState{<:GPSL5Data}, msg::UInt320)
    # 232 text bits at positions 39..270, then 4-bit text page at 271..274.
    text_bits = UInt256(0)
    for i = 0:231
        text_bits = (text_bits << 1) | UInt256(_gpsl5_get_bit(msg, 39 + i))
    end
    text15_page = UInt(_gpsl5_get_bits(msg, 271, 4))
    return GNSSDecoderState(
        state;
        raw_data = GPSL5Data(
            state.raw_data;
            text15_bits = text_bits,
            text15_page,
        ),
    )
end

function _gpsl5_decode_mt36(state::GNSSDecoderState{<:GPSL5Data}, msg::UInt320)
    # 144 text bits at positions 128..271, then 4-bit text page at 272..275.
    text_bits = UInt256(0)
    for i = 0:143
        text_bits = (text_bits << 1) | UInt256(_gpsl5_get_bit(msg, 128 + i))
    end
    text36_page = UInt(_gpsl5_get_bits(msg, 272, 4))
    state = _gpsl5_decode_mt30_37_clock(state, msg)
    return GNSSDecoderState(
        state;
        raw_data = GPSL5Data(
            state.raw_data;
            text36_bits = text_bits,
            text36_page,
        ),
    )
end

function _gpsl5_decode_mt40(state::GNSSDecoderState{<:GPSL5Data}, msg::UInt320)
    GNSS_id = UInt(_gpsl5_get_bits(msg, 39, 4))
    WN_ISM = UInt(_gpsl5_get_bits(msg, 43, 13))
    TOW_ISM = UInt(_gpsl5_get_bits(msg, 56, 6))
    t_correl = UInt(_gpsl5_get_bits(msg, 62, 4))
    b_nom = UInt(_gpsl5_get_bits(msg, 66, 4))
    γ_nom = UInt(_gpsl5_get_bits(msg, 70, 4))
    R_sat = UInt(_gpsl5_get_bits(msg, 74, 4))
    P_const = UInt(_gpsl5_get_bits(msg, 78, 4))
    MFD = UInt(_gpsl5_get_bits(msg, 82, 4))
    service_level = UInt(_gpsl5_get_bits(msg, 86, 3))
    # 63-bit mask split across bits 89..151.
    mask = UInt128(0)
    for i = 0:62
        mask = (mask << 1) | UInt128(_gpsl5_get_bit(msg, 89 + i))
    end
    ism = GPSL5IntegritySupportMessage(;
        GNSS_id, WN_ISM, TOW_ISM, t_correl, b_nom, γ_nom, R_sat, P_const,
        MFD, service_level, mask,
    )
    return GNSSDecoderState(state; raw_data = GPSL5Data(state.raw_data; ism))
end

# Dispatch a complete validated 300-bit message to its type-specific decoder.
function _gpsl5_dispatch_message(state::GNSSDecoderState{<:GPSL5Data}, msg::UInt320)
    state = _gpsl5_decode_common_header(state, msg)
    msg_id = state.raw_data.last_message_id
    if msg_id == 10
        return _gpsl5_decode_mt10(state, msg)
    elseif msg_id == 11
        return _gpsl5_decode_mt11(state, msg)
    elseif msg_id == 12
        return _gpsl5_decode_mt12(state, msg)
    elseif msg_id == 13
        return _gpsl5_decode_mt13(state, msg)
    elseif msg_id == 14
        return _gpsl5_decode_mt14(state, msg)
    elseif msg_id == 15
        return _gpsl5_decode_mt15(state, msg)
    elseif msg_id == 30
        return _gpsl5_decode_mt30(state, msg)
    elseif msg_id == 31
        return _gpsl5_decode_mt31(state, msg)
    elseif msg_id == 32
        return _gpsl5_decode_mt32(state, msg)
    elseif msg_id == 33
        return _gpsl5_decode_mt33(state, msg)
    elseif msg_id == 34
        return _gpsl5_decode_mt34(state, msg)
    elseif msg_id == 35
        return _gpsl5_decode_mt35(state, msg)
    elseif msg_id == 36
        return _gpsl5_decode_mt36(state, msg)
    elseif msg_id == 37
        return _gpsl5_decode_mt37(state, msg)
    elseif msg_id == 40
        return _gpsl5_decode_mt40(state, msg)
    else
        # Unknown / reserved message type — keep the decoded common header but
        # no further fields are extracted.
        return state
    end
end

# =============================================================================
# Decode entry point: feed a stream of 100 sps symbols
# =============================================================================

# Process one received symbol on a single Viterbi/framing instance.
# Returns the updated part and (optionally) a 300-bit message that just passed
# CRC, or `nothing` if no message is ready.
function _gpsl5_step_part(part::GPSL5FramingPart, symbol::UInt8)
    if part.pending_symbol < 0
        # First symbol of a pair — buffer it.
        return GPSL5FramingPart(part; pending_symbol = Int8(symbol)), nothing
    end
    # Second symbol — run a Viterbi step.
    sym1 = UInt8(part.pending_symbol)
    sym2 = symbol
    new_v, decoded_bit = _gpsl5_viterbi_step(part.viterbi, sym1, sym2)
    part = GPSL5FramingPart(part; viterbi = new_v, pending_symbol = Int8(-1))
    if decoded_bit < 0
        return part, nothing
    end
    part = _gpsl5_push_bit(part, decoded_bit)

    # Drive the framing state machine.
    while true
        if !part.preamble_seen
            part = _gpsl5_rescan_preamble(part)
        end
        if part.preamble_seen && part.n_decoded >= 300
            new_part, msg_bits = _gpsl5_extract_message(part)
            if _gpsl5_check_crc(msg_bits)
                # CRC passed.
                if part.message_lock
                    part = GPSL5FramingPart(new_part; n_crc_fail = 0, preamble_seen = false)
                    return part, msg_bits
                else
                    # First successful decode — acquire lock, return message.
                    part = GPSL5FramingPart(
                        new_part;
                        message_lock = true,
                        n_crc_fail = 0,
                        preamble_seen = false,
                    )
                    return part, msg_bits
                end
            else
                # CRC failed.
                if part.message_lock
                    nf = part.n_crc_fail + 1
                    if nf > GPSL5_LOCK_MAX_CRC_FAILS
                        # Drop the lock and re-search for preamble in the buffered bits.
                        part = _gpsl5_drop_lock(GPSL5FramingPart(part; n_crc_fail = 0))
                        # Loop to retry preamble search on whatever's still in the buffer.
                        continue
                    else
                        part = GPSL5FramingPart(new_part; n_crc_fail = nf, preamble_seen = false)
                        return part, nothing
                    end
                else
                    # Spurious preamble — drop the oldest bit so the next
                    # rescan considers a different alignment, then retry.
                    part = _gpsl5_drop_oldest_bit(part)
                    part = GPSL5FramingPart(part; preamble_seen = false)
                    continue
                end
            end
        else
            return part, nothing
        end
    end
end

"""
$(TYPEDSIGNATURES)

Decode a stream of 100 sps GPS L5 CNAV FEC symbols and update the decoder state.

The decoder accepts hard-decision symbols (0 or 1), feeds them through two
parallel persistent Viterbi decoders to resolve symbol-pair-phase ambiguity,
locates the 8-bit preamble (`0b10001011`) in the post-Viterbi 50 bps bit
stream, validates each 300-bit CNAV message with CRC-24Q, and dispatches
validated messages to per-type decoders.

# Arguments
- `state::GNSSDecoderState{<:GPSL5Data}`: Current decoder state
- `symbols::T`: Unsigned integer containing FEC symbols, MSB first
- `num_symbols::Int`: Number of valid symbols in `symbols` to process

# Keywords
- `decode_once::Bool=false`: If `true`, stops decoding after the minimum
  positioning data set (MT 10 + MT 11 + MT 30) has been received.

# Returns
- `GNSSDecoderState`: Updated decoder state with newly decoded message data

# Example
```julia
state = GPSL5DecoderState(1)
state = decode(state, UInt8(0b10110010), 8)  # 8 hard-decision symbols
```

# See Also
- [`GPSL5DecoderState`](@ref): Create a decoder state
- [`reset_decoder_state`](@ref): Reset after signal loss
- [`is_sat_healthy`](@ref): Check satellite health status
"""
function decode(
    state::GNSSDecoderState{<:GPSL5Data},
    symbols::T,
    num_symbols::Int;
    decode_once::Bool = false,
) where {T<:Unsigned}
    num_symbols > sizeof(symbols) * 8 &&
        throw(ArgumentError("Number of symbols is too large to fit type of symbols"))
    cache = state.cache
    for i = num_symbols-1:-1:0
        sym = UInt8((symbols >> i) & T(1))
        # Process on both parts.
        new_part1, msg1 = _gpsl5_step_part(cache.part1, sym)
        new_part2, msg2 = _gpsl5_step_part(cache.part2, sym)
        # Resolve which part wins lock if either succeeds for the first time.
        locked = cache.locked_part
        if locked == 0
            if new_part1.message_lock
                locked = 1
                # Reset the losing part to save memory churn (matches GNSS-SDR).
                new_part2 = GPSL5FramingPart()
            elseif new_part2.message_lock
                locked = 2
                new_part1 = GPSL5FramingPart()
            end
        elseif locked == 1 && !new_part1.message_lock
            # Lost lock on part1 — both parts re-enter free search.
            locked = 0
            new_part2 = GPSL5FramingPart()
        elseif locked == 2 && !new_part2.message_lock
            locked = 0
            new_part1 = GPSL5FramingPart()
        end
        cache = GPSL5Cache(cache; part1 = new_part1, part2 = new_part2, locked_part = locked)
        # Apply any successfully decoded message.
        msg_bits = nothing
        if locked == 1 && msg1 !== nothing
            msg_bits = msg1
        elseif locked == 2 && msg2 !== nothing
            msg_bits = msg2
        elseif locked == 0
            # Both parts may have happened to decode in the same call (e.g. on
            # the very first lock-acquiring message). Prefer part1's message
            # (matches the dispatch order; locked is now 0 since neither part
            # is currently locked here, but msg1/msg2 captured the lock-event).
            if msg1 !== nothing
                msg_bits = msg1
            elseif msg2 !== nothing
                msg_bits = msg2
            end
        end
        # Update num_bits_after_valid_syncro_sequence for caller bookkeeping
        # (counts symbols since last validated message).
        if !isnothing(state.num_bits_after_valid_syncro_sequence)
            state = GNSSDecoderState(
                state;
                num_bits_after_valid_syncro_sequence =
                    state.num_bits_after_valid_syncro_sequence + 1,
            )
        end
        state = GNSSDecoderState(state; cache)

        if msg_bits !== nothing
            state = _gpsl5_dispatch_message(state, msg_bits)
            if !decode_once || !is_decoding_completed_for_positioning(state.data)
                state = validate_data(state)
            end
            state = GNSSDecoderState(
                state;
                num_bits_after_valid_syncro_sequence = 0,
            )
        end
    end
    return state
end

# =============================================================================
# Validation: copy raw_data → data once enough has accumulated for positioning
# =============================================================================

function is_decoding_completed_for_positioning(data::GPSL5Data)
    !isnothing(data.TOW) &&
        !isnothing(data.WN) &&
        # MT10 ephemeris-1
        !isnothing(data.ΔA) &&
        !isnothing(data.A_dot) &&
        !isnothing(data.Δn_0) &&
        !isnothing(data.M_0) &&
        !isnothing(data.e) &&
        !isnothing(data.ω) &&
        # MT11 ephemeris-2
        !isnothing(data.t_oe) &&
        !isnothing(data.Ω_0) &&
        !isnothing(data.i_0) &&
        !isnothing(data.Ω_dot) &&
        !isnothing(data.i_dot) &&
        !isnothing(data.C_is) &&
        !isnothing(data.C_ic) &&
        !isnothing(data.C_rs) &&
        !isnothing(data.C_rc) &&
        !isnothing(data.C_us) &&
        !isnothing(data.C_uc) &&
        # MT30 clock + group delay
        !isnothing(data.t_oc) &&
        !isnothing(data.a_f0) &&
        !isnothing(data.a_f1) &&
        !isnothing(data.a_f2) &&
        !isnothing(data.T_GD)
end

function validate_data(state::GNSSDecoderState{<:GPSL5Data})
    if is_decoding_completed_for_positioning(state.raw_data)
        state = GNSSDecoderState(state; data = state.raw_data)
    end
    return state
end

"""
$(TYPEDSIGNATURES)

Check if the GPS L5 satellite is healthy and usable for positioning.

Examines the L5 health bit decoded from the most recent MT10 message. The
satellite is considered healthy when `signal_health_l5 == false`.

!!! warning
    This function requires that an MT10 message has been successfully decoded.
    Check that `state.data.signal_health_l5` is not `nothing` before relying
    on this result.

# Arguments
- `state::GNSSDecoderState{<:GPSL5Data}`: GPS L5 decoder state with decoded data

# Returns
- `Bool`: `true` if the L5 health bit indicates normal operation
"""
function is_sat_healthy(state::GNSSDecoderState{<:GPSL5Data})
    state.data.signal_health_l5 == false
end
