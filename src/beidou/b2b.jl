# BeiDou B2b (B-CNAV3) navigation message decoder — BDS-SIS-ICD-B2b-1.0
# (2020-07).
#
# The B-CNAV3 message is broadcast on the B2b_I component (1207.14 MHz) of the
# BDS-3 MEO/IGSO satellites at 1000 sps. Each frame is 1000 symbols = 1 s
# (ICD Figure 6-1):
#
#     Pre (16 sym, 0xEB90) | PRN (6 sym, unencoded) | Rev (6 sym) | 972 encoded
#
# The 972 encoded symbols are the 64-ary LDPC(162, 81) codeword of a 486-bit
# message `MesType(6) SOW(20) data(436) CRC(24)` (§6.2.1); MesType, SOW, and
# data participate in the CRC-24Q. The LDPC codeword is decoded through the
# binary image of the ICD's GF(2^6) parity-check matrix (`data/bcnv3.alist`,
# see `scripts/generate_beidou_alist.jl`) with the shared Aff3ct BP helper.
#
# Sync mirrors the shared preamble machinery (16-bit preamble at both ends of
# the 1016-symbol window, either polarity) and is hardened twice before a
# frame is accepted: the 6 unencoded PRN symbols following the preamble must
# hard-decode to the decoder's own PRN (they are known a priori, ICD §6.2.1 —
# a free 6-symbol preamble extension that rejects cross-correlation locks
# onto another satellite), and the LDPC-decoded message must pass CRC-24Q.
# Both gates run inside `try_sync` (mirroring the GPS CNAV window decoder),
# so a failed frame never reaches `decode_syncro_sequence` and can never
# re-arm the SOW counter with stale data.
#
# Message types (§6.2.3, Table 7-1): 10 (ephemeris + integrity flags),
# 30 (clock, TGD, BDGIM ionosphere, BDT-UTC, EOP, SISAI, health), and
# 40 (BGTO, one midi almanac, five reduced almanacs). `000000` is Invalid and
# the remaining values are reserved; both are dropped the same way.
#
# B-CNAV3 is the MEO/IGSO message. The BDS-3 GEO satellites broadcast the
# *PPP-B2b* service on the same component with the same framing (preamble, PRN,
# LDPC block, CRC-24Q), but with its own message types 1-7 carrying precise
# orbit/clock/code-bias corrections — a separate ICD, and not decoded here.
# Nothing special is needed for that: such a satellite syncs and passes CRC
# normally, its message types fall through the reserved branch, and only
# `last_message_type` and `SOW` are published. This mirrors how the legacy D2
# wide-area differential subframes are framed but not parsed (`dnav.jl`).

# ---- Frame constants (BDS-SIS-ICD-B2b-1.0 §6.2.1) ---------------------------

"""
One B-CNAV3 frame: 1000 symbols at 1000 sps = 1 second.
"""
const B2B_FRAME_SYMBOLS = 1000
"""
Leading preamble 0xEB90 = `1110101110010000` (16 symbols, MSB first).
"""
const B2B_PREAMBLE = 0xEB90
const B2B_PREAMBLE_SYMBOLS = 16
"""
Unencoded 6-symbol PRN field following the preamble.
"""
const B2B_PRN_SYMBOLS = 6
"""
Reserved 6 symbols between the PRN field and the encoded message.
"""
const B2B_REV_SYMBOLS = 6
"""
Encoded message: 162 GF(2^6) symbols x 6 bits = 972 channel symbols.
"""
const B2B_ENCODED_SYMBOLS = 972
"""
Message length before encoding: MesType(6) + SOW(20) + data(436) + CRC(24).
"""
const B2B_MESSAGE_BITS = 486
"""
Sync window: one frame plus the next frame's preamble.
"""
const B2B_WINDOW_SYMBOLS = B2B_FRAME_SYMBOLS + B2B_PREAMBLE_SYMBOLS  # 1016
"""
Offset of the first encoded symbol inside the frame (after Pre + PRN + Rev).
"""
const B2B_ENCODED_OFFSET = B2B_PREAMBLE_SYMBOLS + B2B_PRN_SYMBOLS + B2B_REV_SYMBOLS  # 28

# Semi-major axis reference values for the ephemeris ΔA and the reduced
# almanac δA (ICD Table 7-6 / Table 7-14 footnotes, meters).
const B2B_A_REF_MEO = 27_906_100.0
const B2B_A_REF_IGSO_GEO = 42_162_200.0

"""
$(TYPEDEF)

Constants for the BeiDou B2b (B-CNAV3) decoder — BDS-SIS-ICD-B2b-1.0.

# Fields

$(TYPEDFIELDS)
"""
Base.@kwdef struct BeiDouB2bConstants <: AbstractGNSSConstants
    """
    Frame length drained after each decoded frame (1000 symbols)
    """
    syncro_sequence_length::Int = B2B_FRAME_SYMBOLS
    """
    Preamble 0xEB90 (ICD §6.2.1), MSB first
    """
    preamble::UInt16 = B2B_PREAMBLE
    """
    Trailing next-frame preamble segment retained for sync (16 symbols)
    """
    preamble_length::Int = B2B_PREAMBLE_SYMBOLS
    """
    Mathematical constant π (BDS-SIS-ICD-B2b-1.0 Table 7-13)
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
    Relativistic correction constant F = -2√μ/c² (s/√m, ICD §7.4.2)
    """
    F::Float64 = -4.442807309e-10
end

# ---- Almanac records (ICD §7.8 / §7.9) --------------------------------------

# ---- Decoded data container --------------------------------------------------

"""
    BeiDouB2bData <: AbstractBeiDouData

Decoded BeiDou B-CNAV3 navigation data (BDS-SIS-ICD-B2b-1.0 §6.2.3 / §7).

Every field is `nothing` until first decoded from a CRC-validated frame.
Message type 10 carries a complete ephemeris in one frame, message type 30 a
complete clock set; B-CNAV3 broadcasts no issue-of-data stamps, so no
cross-message consistency vote is needed (or possible) before promotion —
each CRC-gated message is atomic.

# Header (every message type)

  - `last_message_type::Int`: MesType of the most recently decoded frame (Table 7-1; 0 until the first decode).
  - `SOW::Int64`: Seconds of week of the current frame's leading edge (s). Broadcast as a 20-bit count with LSB 1 s — the scale factor of the Chinese-language Table 7-2; the English edition misprints it as 3. See `decode_syncro_sequence` for why the English table is not followed.

# Message type 10 — ephemeris (Figures 6-3, 6-6, 6-7) and integrity flags

  - `sat_type::Int64`: Satellite orbit type (binary 01 = GEO, 10 = IGSO, 11 = MEO; Table 7-6).
  - `t_0e::Int64`: Ephemeris reference time (s, LSB 300; the ICD writes `toe`).
  - `ΔA::Float64`: Semi-major axis difference at reference time (m) relative to `A_ref = 27 906 100 m` (MEO) / `42 162 200 m` (IGSO/GEO).
  - `A_dot::Float64`: Change rate of semi-major axis (m/s).
  - `Δn_0::Float64`: Mean motion difference at reference time (rad/s).
  - `Δn_0_dot::Float64`: Rate of mean motion difference (rad/s²).
  - `M_0::Float64`: Mean anomaly at reference time (rad).
  - `e::Float64`: Eccentricity.
  - `ω::Float64`: Argument of perigee (rad).
  - `Ω_0::Float64`: Longitude of ascending node at weekly epoch (rad).
  - `i_0::Float64`: Inclination at reference time (rad).
  - `Ω_dot::Float64`: Rate of right ascension (rad/s) — broadcast in full, not as a difference.
  - `i_dot::Float64`: Rate of inclination (rad/s).
  - `C_is,C_ic::Float64`: Sine/cosine harmonic correction to inclination (rad).
  - `C_rs,C_rc::Float64`: Sine/cosine harmonic correction to orbit radius (m).
  - `C_us,C_uc::Float64`: Sine/cosine harmonic correction to argument of latitude (rad).
  - `DIF,SIF,AIF::Bool`: Data / signal / accuracy integrity flags for the B2b_I signal (`false` = OK, Table 7-21).
  - `SISMAI::Int64`: Signal in space monitoring accuracy index (4 bits; semantics deferred to a future ICD update, §7.16). B2b-1.0 Figure 6-3 labels the broadcast field `SISMA`, reserving `SISMAI` for its index; B1C-1.0 and B2a-1.0 label the same field `SISMAI`, which is the spelling this package follows.

# Message type 30 — clock, group delay, ionosphere, UTC, EOP, accuracy, health

  - `WN::Int64`: BDT week number at the current frame's epoch (Table 7-2).
  - `t_0c::Int64`: Clock correction reference time (s, LSB 300).
  - `a_f0::Float64`: Clock bias (s). `a_f1::Float64`: drift (s/s). `a_f2::Float64`: drift rate (s/s²). (Table 7-3; the ICD names them a0/a1/a2, and writes `toc` for `t_0c`.)
  - `T_GD_B2bI::Float64`: Group delay differential of the B2b_I signal relative to B3I (s, Table 7-4). Required, not optional: `a_f0` is referenced to B3I (§7.6), so a B2b_I receiver must add this to reach its own component. It arrives in the same MT30 as the clock, so `is_decoding_completed_for_positioning` gates on it at no cost in time to first fix.
  - `α_bdgim_1 … α_bdgim_9::Float64`: BDGIM ionospheric model parameters (TECu, Table 7-8; α_bdgim_5 is broadcast unsigned and the decoder applies its negative scale factor, -2⁻³). The ICD names them `α_1` … `α_9`; see `b1c.jl` for why the `bdgim` qualifier is added.
  - `A_0UTC,A_1UTC,A_2UTC::Float64`: BDT-UTC polynomial (s, s/s, s/s²; Table 7-18).
  - `Δt_LS::Int64`: Current or past leap-second count (s). `Δt_LSF::Int64`: current or future leap-second count (s).
  - `t_0t::Int64`: UTC reference time of week (s, LSB 2⁴). `WN_0t::Int64`: UTC reference week. (The ICD writes `t_ot`/`WN_ot`.)
  - `WN_LSF::Int64`: Leap-second reference week. `DN::Int64`: leap-second reference day (0-6).
  - `t_EOP::Int64`: EOP reference time (s, LSB 2⁴; Table 7-16).
  - `PM_X,PM_Y::Float64`: Polar motion (arc-seconds). `PM_X_dot,PM_Y_dot::Float64`: drift (arc-seconds/day).
  - `ΔUT1::Float64`: UT1-UTC difference (s). `ΔUT1_dot::Float64`: its rate (s/day).
  - `t_op,SISAI_ocb,SISAI_oc1,SISAI_oc2,SISAI_oe::Int64`: Signal-in-space accuracy index fields, raw broadcast values (11/5/3/3/5 bits; semantics deferred to a future ICD update, §7.15).
  - `HS::Int64`: Satellite health status (0 = healthy, 1 = unhealthy or in test, 2-3 reserved; Table 7-20).

# Message type 40 — BGTO and almanacs (Figure 6-5)

  - `GNSS_ID::Int64`: BGTO GNSS identification (0 = not available, 1 = GPS, 2 = Galileo, 3 = GLONASS; §7.12).
  - `WN_0BGTO::Int64`, `t_0BGTO::Int64`: BGTO reference week / time of week (s, LSB 2⁴).
  - `A_0BGTO,A_1BGTO,A_2BGTO::Float64`: BDT-GNSS time offset polynomial (s, s/s, s/s²; Table 7-19).
  - `midi_almanacs::Dictionary{Int,BeiDouMidiAlmanac}`: Midi almanacs keyed by `PRN_a` (§7.8).
  - `reduced_almanacs::Dictionary{Int,BeiDouReducedAlmanac}`: Reduced almanacs keyed by `PRN_a` (§7.9).
  - `WN_a::Int64`, `t_0a::Int64`: Almanac reference week / time (s, LSB 2¹²) for the reduced almanacs (Table 7-15).

# Reference

BDS-SIS-ICD-B2b-1.0, Figures 6-3 through 6-15 and Tables 7-1 through 7-21.
"""
Base.@kwdef struct BeiDouB2bData <: AbstractBeiDouData
    last_message_type::Int = 0
    SOW::Union{Nothing,Int64} = nothing

    # Message type 10: ephemeris + integrity
    sat_type::Union{Nothing,Int64} = nothing
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
    Ω_dot::Union{Nothing,Float64} = nothing
    i_dot::Union{Nothing,Float64} = nothing
    C_is::Union{Nothing,Float64} = nothing
    C_ic::Union{Nothing,Float64} = nothing
    C_rs::Union{Nothing,Float64} = nothing
    C_rc::Union{Nothing,Float64} = nothing
    C_us::Union{Nothing,Float64} = nothing
    C_uc::Union{Nothing,Float64} = nothing
    DIF::Union{Nothing,Bool} = nothing
    SIF::Union{Nothing,Bool} = nothing
    AIF::Union{Nothing,Bool} = nothing
    SISMAI::Union{Nothing,Int64} = nothing

    # Message type 30: clock + TGD + ionosphere + UTC + EOP + SISAI + health
    WN::Union{Nothing,Int64} = nothing
    t_0c::Union{Nothing,Int64} = nothing
    a_f0::Union{Nothing,Float64} = nothing
    a_f1::Union{Nothing,Float64} = nothing
    a_f2::Union{Nothing,Float64} = nothing
    T_GD_B2bI::Union{Nothing,Float64} = nothing
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
    t_EOP::Union{Nothing,Int64} = nothing
    PM_X::Union{Nothing,Float64} = nothing
    PM_X_dot::Union{Nothing,Float64} = nothing
    PM_Y::Union{Nothing,Float64} = nothing
    PM_Y_dot::Union{Nothing,Float64} = nothing
    ΔUT1::Union{Nothing,Float64} = nothing
    ΔUT1_dot::Union{Nothing,Float64} = nothing
    t_op::Union{Nothing,Int64} = nothing
    SISAI_ocb::Union{Nothing,Int64} = nothing
    SISAI_oc1::Union{Nothing,Int64} = nothing
    SISAI_oc2::Union{Nothing,Int64} = nothing
    SISAI_oe::Union{Nothing,Int64} = nothing
    HS::Union{Nothing,Int64} = nothing

    # Message type 40: BGTO + almanacs
    GNSS_ID::Union{Nothing,Int64} = nothing
    WN_0BGTO::Union{Nothing,Int64} = nothing
    t_0BGTO::Union{Nothing,Int64} = nothing
    A_0BGTO::Union{Nothing,Float64} = nothing
    A_1BGTO::Union{Nothing,Float64} = nothing
    A_2BGTO::Union{Nothing,Float64} = nothing
    midi_almanacs::Union{Nothing,Dictionary{Int,BeiDouMidiAlmanac}} = nothing
    reduced_almanacs::Union{Nothing,Dictionary{Int,BeiDouReducedAlmanac}} = nothing
    WN_a::Union{Nothing,Int64} = nothing
    t_0a::Union{Nothing,Int64} = nothing
end

# Field-by-field equality: the almanac `Dictionary` fields otherwise compare
# by identity through the default struct `==`.
Base.:(==)(a::BeiDouB2bData, b::BeiDouB2bData) = fields_equal(a, b)

function BeiDouB2bData(
    data::BeiDouB2bData;
    last_message_type = data.last_message_type,
    SOW = data.SOW,
    sat_type = data.sat_type,
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
    Ω_dot = data.Ω_dot,
    i_dot = data.i_dot,
    C_is = data.C_is,
    C_ic = data.C_ic,
    C_rs = data.C_rs,
    C_rc = data.C_rc,
    C_us = data.C_us,
    C_uc = data.C_uc,
    DIF = data.DIF,
    SIF = data.SIF,
    AIF = data.AIF,
    SISMAI = data.SISMAI,
    WN = data.WN,
    t_0c = data.t_0c,
    a_f0 = data.a_f0,
    a_f1 = data.a_f1,
    a_f2 = data.a_f2,
    T_GD_B2bI = data.T_GD_B2bI,
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
    t_op = data.t_op,
    SISAI_ocb = data.SISAI_ocb,
    SISAI_oc1 = data.SISAI_oc1,
    SISAI_oc2 = data.SISAI_oc2,
    SISAI_oe = data.SISAI_oe,
    HS = data.HS,
    GNSS_ID = data.GNSS_ID,
    WN_0BGTO = data.WN_0BGTO,
    t_0BGTO = data.t_0BGTO,
    A_0BGTO = data.A_0BGTO,
    A_1BGTO = data.A_1BGTO,
    A_2BGTO = data.A_2BGTO,
    midi_almanacs = data.midi_almanacs,
    reduced_almanacs = data.reduced_almanacs,
    WN_a = data.WN_a,
    t_0a = data.t_0a,
)
    BeiDouB2bData(
        last_message_type,
        SOW,
        sat_type,
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
        Ω_dot,
        i_dot,
        C_is,
        C_ic,
        C_rs,
        C_rc,
        C_us,
        C_uc,
        DIF,
        SIF,
        AIF,
        SISMAI,
        WN,
        t_0c,
        a_f0,
        a_f1,
        a_f2,
        T_GD_B2bI,
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
        t_op,
        SISAI_ocb,
        SISAI_oc1,
        SISAI_oc2,
        SISAI_oe,
        HS,
        GNSS_ID,
        WN_0BGTO,
        t_0BGTO,
        A_0BGTO,
        A_1BGTO,
        A_2BGTO,
        midi_almanacs,
        reduced_almanacs,
        WN_a,
        t_0a,
    )
end

# ---- Cache -------------------------------------------------------------------

"""
$(TYPEDEF)

Mutable per-decoder scratch for BeiDou B2b: the streaming soft-symbol buffer,
the LDPC window copied out per sync attempt, and the long-lived Aff3ct LDPC
belief-propagation decoder handle for the B-CNAV3 binary image (972 → 486
bits), shared by reference through the otherwise-immutable state.

# Fields

$(TYPEDFIELDS)
"""
struct BeiDouB2bCache <: AbstractGNSSCache
    """
    Soft-symbol buffer (1016 = 1000 frame + 16 next-frame preamble)
    """
    soft_buffer::CircularDeque{Float32}
    """
    Polarity-resolved 972-symbol LDPC window copied out per sync attempt
    """
    ldpc_window::Vector{Float32}
    """
    LDPC decoder and its scratch buffers for the B-CNAV3 binary image (K=486, N=972)
    """
    ldpc::LDPCScratch
end

BeiDouB2bCache() = BeiDouB2bCache(
    CircularDeque{Float32}(B2B_WINDOW_SYMBOLS),
    Vector{Float32}(undef, B2B_ENCODED_SYMBOLS),
    LDPCScratch(alist_path("bcnv3.alist")),
)

# The LDPC decoder handle is stateless w.r.t. equality (it is a runtime Aff3ct
# object); two B2b caches are equal when their soft buffers match.
function Base.:(==)(a::BeiDouB2bCache, b::BeiDouB2bCache)
    deques_equal(a.soft_buffer, b.soft_buffer)
end

# ---- Sync --------------------------------------------------------------------

# No `packed_buffer_type` method: B2b overrides `try_sync` and reads the few
# bits it needs straight from the soft buffer, so the packed window this
# would name is never built.

"""
$(TYPEDEF)

Result of a successful B-CNAV3 frame sync: the CRC-validated 486-bit message
and the resolved polarity. Produced by `try_sync`, threaded through
`complement_buffer_if_necessary` to `decode_syncro_sequence`.

# Fields

$(TYPEDFIELDS)
"""
struct BeiDouB2bSync
    """
    CRC-validated message, packed MSB-first (bit 1 of the message at the MSB)
    """
    word::UInt512
    """
    Whether the symbol stream is 180-degrees phase shifted
    """
    polarity_flipped::Bool
end

"""
    try_sync(state::GNSSDecoderState{<:BeiDouB2bData}) -> Union{Nothing,BeiDouB2bSync}

B-CNAV3 frame sync (BDS-SIS-ICD-B2b-1.0 §6.2.1). Three gates, cheapest first:

 1. The 16-bit preamble `0xEB90` must appear at both ends of the 1016-symbol
    window (start of this frame and start of the next), both upright or both
    inverted (the shared `find_preamble` rule).
 2. The 6 unencoded PRN symbols following the preamble must hard-decode (in
    the resolved polarity) to the decoder's own PRN — they are known a
    priori, so they act as a free 6-symbol preamble extension.
 3. The 972 encoded symbols must LDPC-decode to a message passing CRC-24Q.

Only then is the frame handed to `decode_syncro_sequence`, so a corrupted or
foreign frame can never update `raw_data` (and, via `validate_data`, can
never re-arm the SOW symbol counter with stale time).
"""
function try_sync(state::GNSSDecoderState{<:BeiDouB2bData})
    deque = soft_buffer(state)
    # Gate 1, read straight from the deque via the shared `find_preamble_in_deque`.
    # The default `pack_buffer` path would shift a 1016-bit integer once per
    # symbol — ~16 M word-operations per second per satellite at B2b's 1000 sps —
    # to expose 38 bits that sit at three known offsets. Slice those instead; the
    # packed window is never needed, since the payload is consumed as soft
    # symbols below.
    polarity_flipped = find_preamble_in_deque(
        deque,
        state.constants.preamble,
        B2B_PREAMBLE_SYMBOLS,
        B2B_FRAME_SYMBOLS,
    )
    isnothing(polarity_flipped) && return nothing
    # Gate 2: the unencoded PRN field (window bits 17-22, ICD Figure 6-1).
    prn_bits = pack_bits_msb_first(deque, B2B_PREAMBLE_SYMBOLS + 1, B2B_PRN_SYMBOLS)
    prn_mask = (UInt64(1) << B2B_PRN_SYMBOLS) - UInt64(1)
    prn = Int(polarity_flipped ? prn_bits ⊻ prn_mask : prn_bits)
    prn == state.prn || return nothing
    # Gate 3: LDPC + CRC on the soft symbols (deque indices 29 .. 1000).
    window = copy_soft_window!(
        state.cache.ldpc_window,
        deque,
        B2B_ENCODED_OFFSET,
        B2B_ENCODED_SYMBOLS,
        polarity_flipped,
    )
    word = ldpc_decode_word(state.cache.ldpc, window, UInt512)
    isnothing(word) && return nothing
    return BeiDouB2bSync(word, polarity_flipped)
end

"""
    complement_buffer_if_necessary(state::GNSSDecoderState{<:BeiDouB2bData}, sync)

Record the polarity resolved by `try_sync` on the state and pass the
[`BeiDouB2bSync`](@ref) through unchanged (its message bits are already
polarity-resolved and CRC-validated).
"""
function complement_buffer_if_necessary(
    state::GNSSDecoderState{<:BeiDouB2bData},
    sync::BeiDouB2bSync,
)
    GNSSDecoderState(state; is_shifted_by_180_degrees = sync.polarity_flipped), sync
end

# ---- Message parsing (BDS-SIS-ICD-B2b-1.0 §6.2.3 / §7) ------------------------
#
# `word` is the 486-bit message packed MSB-first into a `UInt512` (bit 1 = the
# MesType MSB). Fields are read by 1-based start bit and length through the
# shared `get_bits` / `get_twos_complement_num` / `get_bit` helpers. Bit
# positions follow Figures 6-3 (MT10), 6-4 (MT30), and 6-5 (MT40) with the
# data-block layouts of Figures 6-6 through 6-15.

"""
    decode_syncro_sequence(state::GNSSDecoderState{<:BeiDouB2bData}, sync::BeiDouB2bSync)

Parse one CRC-validated B-CNAV3 message into `raw_data`. The 20-bit SOW
(LSB 1 s — the Chinese-language ICD Table 7-2 is authoritative here; see the
provenance note in the body before "correcting" it) stamps the rising edge of
the current frame's preamble. Message types other than 10, 30, and 40 are reserved (Table 7-1)
and contribute only the header.
"""
function decode_syncro_sequence(
    state::GNSSDecoderState{<:BeiDouB2bData},
    sync::BeiDouB2bSync,
)
    word = sync.word
    word_length = B2B_MESSAGE_BITS
    PI = state.constants.PI

    message_type = Int(get_bits(word, word_length, 1, 6))
    # 20-bit SOW count in seconds (LSB 1 s), stamping the current frame's
    # leading edge (Table 7-2: effective range 0..604799 s). One frame is
    # 1000 symbols = 1 s, so consecutive frames step the field by exactly 1;
    # 2^20 = 1048576 spans the week.
    #
    # The scale factor 1 is authoritative from the Chinese-language original:
    # BDS-SIS-ICD-B2b-1.0 中文版 (2020-07) 表 7-2 prints 比例因子 (scale
    # factor) 1. The English edition's Table 7-2 misprints it as 3 — a slip
    # carried over from B2a's B-CNAV2, whose 18-bit SOW genuinely needs the
    # 3 s LSB (2^18 < 604800) — and the misprint survives the silently
    # re-uploaded 2023 English PDF. The English table is also internally
    # inconsistent: its own effective range 0..604799 is not a multiple of 3
    # and is unreachable under a 3 s LSB. Do not "fix" this back to the
    # English table's literal value.
    SOW = Int64(get_bits(word, word_length, 7, 20))
    raw = BeiDouB2bData(state.raw_data; last_message_type = message_type, SOW)

    raw = if message_type == 10
        parse_b2b_mt10(raw, word, PI)
    elseif message_type == 30
        parse_b2b_mt30(raw, word, PI)
    elseif message_type == 40
        parse_b2b_mt40(raw, word, PI)
    else
        raw  # reserved/invalid message type (Table 7-1): header only
    end
    # SOW refers to the start of the frame just decoded; that epoch lies
    # The symbol counter is deliberately NOT re-armed here (the same split as
    # `b2a.jl`, for the same reason): `get_time_of_week` reads the *validated*
    # `data.SOW`, so the counter must anchor to the last frame `validate_data`
    # promoted, not to the last frame that passed CRC. Re-arming per frame while
    # the published SOW stands still makes `tow + num_bits/rate` read a whole
    # frame low whenever the completeness gate skips a frame — on B2b the gate
    # normally holds once MT10 and MT30 are in, so this is defensive here where
    # on B2a (whose MT10/MT11 adjacency gate skips one frame per message cycle)
    # it was a live once-per-cycle 3-second error.
    GNSSDecoderState(state; raw_data = raw)
end

"""
Message type 10 — complete ephemeris + integrity flags (Figures 6-3, 6-6, 6-7).
"""
function parse_b2b_mt10(raw::BeiDouB2bData, word::UInt512, PI::Float64)
    word_length = B2B_MESSAGE_BITS
    BeiDouB2bData(
        raw;
        # bits 27-30 reserved. Ephemeris I (bits 31-233, Figure 6-6):
        t_0e = Int64(get_bits(word, word_length, 31, 11)) * 300,
        sat_type = Int64(get_bits(word, word_length, 42, 2)),
        ΔA = get_twos_complement_num(word, word_length, 44, 26) * 2.0^-9,
        A_dot = get_twos_complement_num(word, word_length, 70, 25) * 2.0^-21,
        Δn_0 = get_twos_complement_num(word, word_length, 95, 17) * 2.0^-44 * PI,
        Δn_0_dot = get_twos_complement_num(word, word_length, 112, 23) * 2.0^-57 * PI,
        M_0 = get_twos_complement_num(word, word_length, 135, 33) * 2.0^-32 * PI,
        e = Float64(get_bits(word, word_length, 168, 33)) * 2.0^-34,
        ω = get_twos_complement_num(word, word_length, 201, 33) * 2.0^-32 * PI,
        # Ephemeris II (bits 234-455, Figure 6-7):
        Ω_0 = get_twos_complement_num(word, word_length, 234, 33) * 2.0^-32 * PI,
        i_0 = get_twos_complement_num(word, word_length, 267, 33) * 2.0^-32 * PI,
        Ω_dot = get_twos_complement_num(word, word_length, 300, 19) * 2.0^-44 * PI,
        i_dot = get_twos_complement_num(word, word_length, 319, 15) * 2.0^-44 * PI,
        C_is = get_twos_complement_num(word, word_length, 334, 16) * 2.0^-30,
        C_ic = get_twos_complement_num(word, word_length, 350, 16) * 2.0^-30,
        C_rs = get_twos_complement_num(word, word_length, 366, 24) * 2.0^-8,
        C_rc = get_twos_complement_num(word, word_length, 390, 24) * 2.0^-8,
        C_us = get_twos_complement_num(word, word_length, 414, 21) * 2.0^-30,
        C_uc = get_twos_complement_num(word, word_length, 435, 21) * 2.0^-30,
        # Integrity flags + SISMAI (Figure 6-3, Table 7-21):
        DIF = get_bit(word, word_length, 456),
        SIF = get_bit(word, word_length, 457),
        AIF = get_bit(word, word_length, 458),
        SISMAI = Int64(get_bits(word, word_length, 459, 4)),
    )
end

"""
Message type 30 — clock, T_GD, BDGIM ionosphere, BDT-UTC, EOP, SISAI, health
(Figures 6-4, 6-8 through 6-11, 6-13).
"""
function parse_b2b_mt30(raw::BeiDouB2bData, word::UInt512, PI::Float64)
    word_length = B2B_MESSAGE_BITS
    BeiDouB2bData(
        raw;
        WN = Int64(get_bits(word, word_length, 27, 13)),
        # bits 40-43 reserved. Clock correction (bits 44-112, Figure 6-8):
        beidou_clock_block(word, word_length, 44)...,
        T_GD_B2bI = get_twos_complement_num(word, word_length, 113, 12) * 2.0^-34,
        # BDGIM ionosphere (bits 125-198, Figure 6-10, Table 7-8):
        beidou_bdgim_block(word, word_length, 125)...,
        # BDT-UTC (bits 199-295, Figure 6-11, Table 7-18):
        beidou_bdt_utc_block(word, word_length, 199)...,
        # EOP (bits 296-433, Figure 6-13, Table 7-16):
        beidou_eop_block(word, word_length, 296)...,
        # `t_op` and the SISAIoc triple occupy bits 434-455 (Figure 6-9, "Bit
        # allocation for SISAIoc"); `SISAI_oe` at 456-460 is its own field of
        # Figure 6-4, not part of that block. Raw broadcast values throughout —
        # their semantics are deferred to a future ICD update (§7.15).
        beidou_sisai_oc_block(word, word_length, 434)...,
        SISAI_oe = Int64(get_bits(word, word_length, 456, 5)),
        HS = Int64(get_bits(word, word_length, 461, 2)),
    )
end

"""
Message type 40 — BGTO, one midi almanac, five reduced almanacs
(Figures 6-5, 6-12, 6-14, 6-15).
"""
function parse_b2b_mt40(raw::BeiDouB2bData, word::UInt512, PI::Float64)
    word_length = B2B_MESSAGE_BITS
    # Almanac reference week/time for the five reduced almanacs (bits 251-271,
    # Table 7-15) — kept on the data container as the raw broadcast and copied
    # into every reduced-almanac record so each record is complete in itself.
    WN_a_reduced = Int64(get_bits(word, word_length, 251, 13))
    t_0a_reduced = Int64(get_bits(word, word_length, 264, 8)) * 4096
    raw = BeiDouB2bData(
        raw;
        # BGTO (bits 27-94, Figure 6-14, Table 7-19):
        beidou_bgto_block(word, word_length, 27)...,
        # Almanac reference time for the five reduced almanacs (Table 7-15):
        WN_a = WN_a_reduced,
        t_0a = t_0a_reduced,
    )
    # Midi almanac (bits 95-250, Figure 6-15, Table 7-11). PRN_a = 0 marks an
    # empty block (no almanac broadcast in this frame).
    alm = beidou_midi_almanac(word, word_length, 95, PI)
    if !isnothing(alm)
        raw = BeiDouB2bData(
            raw;
            midi_almanacs = _merge_keyed(raw.midi_almanacs, alm.PRN_a, alm),
        )
    end
    # Five 38-bit reduced almanacs (bits 272-461, Figure 6-12, Table 7-14).
    reduced = raw.reduced_almanacs
    for j = 0:4
        packet = beidou_reduced_almanac(
            word,
            word_length,
            272 + 38j,
            WN_a_reduced,
            t_0a_reduced,
            PI,
        )
        isnothing(packet) && continue  # empty block
        reduced = _merge_keyed(reduced, packet.PRN_a, packet)
    end
    BeiDouB2bData(raw; reduced_almanacs = reduced)
end

# ---- Validation --------------------------------------------------------------

function is_ephemeris_decoded(data::BeiDouB2bData)
    # Message type 10 carries the complete ephemeris in one CRC-gated frame.
    !isnothing(data.t_0e) &&
        is_known_sat_type(data.sat_type) &&
        !isnothing(data.ΔA) &&
        !isnothing(data.A_dot) &&
        !isnothing(data.Δn_0) &&
        !isnothing(data.Δn_0_dot) &&
        !isnothing(data.M_0) &&
        !isnothing(data.e) &&
        !isnothing(data.ω) &&
        !isnothing(data.Ω_0) &&
        !isnothing(data.i_0) &&
        !isnothing(data.Ω_dot) &&
        !isnothing(data.i_dot) &&
        !isnothing(data.C_is) &&
        !isnothing(data.C_ic) &&
        !isnothing(data.C_rs) &&
        !isnothing(data.C_rc) &&
        !isnothing(data.C_us) &&
        !isnothing(data.C_uc)
end

function is_clock_correction_decoded(data::BeiDouB2bData)
    # Message type 30, and `T_GD_B2bI` with it. It is required rather than
    # optional for the B1C reason (see `b1c.jl`): ICD §7.6 references a0 to the
    # B3I signal, so a B2b_I receiver must add `T_GD_B2bI` to reach its own
    # component, and the field is decoded by the same constructor call as the
    # clock out of the same CRC-24Q-protected MT30. Requiring it costs no extra
    # wait — the alternative is a consumer treating `nothing` as zero and eating
    # a metre-scale bias.
    !isnothing(data.t_0c) &&
        !isnothing(data.a_f0) &&
        !isnothing(data.a_f1) &&
        !isnothing(data.a_f2) &&
        !isnothing(data.T_GD_B2bI)
end

# Orbit class from the satellite's own broadcast `sat_type` (Ephemeris I, MT10);
# `nothing` until MT10 is decoded, and for the reserved code.
get_orbit_class(state::GNSSDecoderState{<:BeiDouB2bData}) =
    beidou_orbit_class(state.data.sat_type)

# B-CNAV3 broadcasts the seconds of week as a 20-bit count with LSB 1 s (see
# `decode_syncro_sequence` on why the English Table 7-2 scale of 3 is a misprint).
get_time_of_week(data::BeiDouB2bData) = data.SOW

"""
$(TYPEDSIGNATURES)

The BDT-to-`target` offset from message type 40, or `nothing` when this
satellite is not currently broadcasting one for `target`. One set at a time,
tagged by `GNSS_ID`; see `beidou_bgto_offset`.
"""
get_time_offset(state::GNSSDecoderState{<:BeiDouB2bData}, target::TimeSystem) =
    beidou_bgto_offset(state, target)

function is_decoding_completed_for_positioning(data::BeiDouB2bData)
    !isnothing(data.SOW) &&
        !isnothing(data.WN) &&
        is_ephemeris_decoded(data) &&
        is_clock_correction_decoded(data)
end

"""
    validate_data(state::GNSSDecoderState{<:BeiDouB2bData})

Promote `raw_data` to `data` once the minimum positioning set (the message
type 10 ephemeris plus the message type 30 clock and week number) has been
decoded. B-CNAV3 broadcasts no issue-of-data stamps, so there is no
cross-message consistency vote: every message is atomic behind its CRC.

Every frame carries its own SOW, so each validated frame re-arms the
streaming symbol counter: the SOW stamps the *current* frame's leading edge
(ICD Table 7-2), and by validation time that whole frame (1000 symbols) plus
the buffered next-frame preamble (16 symbols) have elapsed since that epoch —
`SOW + num_bits_after_valid_syncro_sequence / 1000 Hz` is the current time.
"""
function validate_data(state::GNSSDecoderState{<:BeiDouB2bData})
    if is_decoding_completed_for_positioning(state.raw_data)
        # Promotion publishes this frame's SOW; the symbol counter re-anchors to
        # this frame's sync epoch with it, so the pair `data.SOW` + counter that
        # `get_time_of_week` consumers combine always describes one frame. Runs
        # right after `decode_syncro_sequence`, so `raw_data.SOW` is the
        # just-decoded frame's.
        return GNSSDecoderState(
            state;
            data = state.raw_data,
            num_bits_after_valid_syncro_sequence = state.constants.syncro_sequence_length +
                                                   state.constants.preamble_length,
        )
    end
    return state
end

function reset_decoder_state(state::GNSSDecoderState{<:BeiDouB2bData})
    empty!(state.cache.soft_buffer)
    GNSSDecoderState(
        state;
        raw_data = BeiDouB2bData(state.raw_data; SOW = nothing),
        data = BeiDouB2bData(),
        num_bits_after_valid_syncro_sequence = nothing,
        is_shifted_by_180_degrees = false,
    )
end

# ---- Signal layer -------------------------------------------------------------

"""
$(TYPEDSIGNATURES)

Create a decoder state for BeiDou B2b (B-CNAV3) navigation messages.

Initializes a [`GNSSDecoderState`](@ref) configured for decoding B-CNAV3
messages from the 1000 sps soft symbols of the B2b_I component
(BDS-SIS-ICD-B2b-1.0). Each sync attempt matches the 16-symbol preamble
`0xEB90` at both ends of the 1016-symbol window, checks the 6 unencoded PRN
symbols against `prn`, LDPC-decodes the 972 encoded symbols through the
binary image of the ICD's 64-ary LDPC(162, 81) code, and validates the
486-bit message with CRC-24Q before dispatching it to the per-type parsers
(message types 10, 30, and 40). Decoded fields land in a
[`BeiDouB2bData`](@ref).

# Arguments

  - `prn::Int`: Pseudo-Random Noise code identifier (6-58 — the B2b_I ranging
    codes of Table 5-1. The wider 1-63 range is the almanac's `PRN_a` key, not a
    trackable B2b channel.)

# Returns

  - `GNSSDecoderState{BeiDouB2bData}`: Initialized decoder state for BeiDou B2b

# Example

```julia
state = BeiDouB2bDecoderState(26)        # PRN 26
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
function BeiDouB2bDecoderState(prn)
    1 <= prn <= 63 || throw(ArgumentError("BeiDou PRN must be in 1..63"))
    GNSSDecoderState(
        prn,
        BeiDouB2bData(),
        BeiDouB2bData(),
        BeiDouB2bConstants(),
        BeiDouB2bCache(),
        nothing,
        false,
    )
end

# Dispatch from a GNSSSignals system type, mirroring `GNSSDecoderState(::GPSL5I, …)`.
# The B2b ICD v1.0 specifies the open-service I component only, so `BeiDouB2bI`
# is the data-bearing signal.
function GNSSDecoderState(system::BeiDouB2bI, prn)
    BeiDouB2bDecoderState(prn)
end

# The signal this decoder demodulates — the B2b_I open-service data component.
# Signal metadata is forwarded through it (see `src/gnss.jl`).
get_signal_type(::BeiDouB2bConstants) = BeiDouB2bI

"""
$(TYPEDSIGNATURES)

Check if the BeiDou B2b satellite is healthy and usable for positioning.

Examines the 2-bit satellite health status (HS) decoded from the most recent
message type 30 (BDS-SIS-ICD-B2b-1.0 §7.13, Table 7-20): a satellite is
healthy iff `HS = 0` ("the satellite is healthy / provides services"); 1
means unhealthy or in test, 2-3 are reserved.

!!! warning

    Requires message type 30 to have been decoded and the positioning set to
    have been validated; returns `false` until then.

# Arguments

  - `state::GNSSDecoderState{<:BeiDouB2bData}`: BeiDou B2b decoder state.

# Returns

  - `Bool`: `true` iff the health status word indicates a healthy satellite.
"""
function is_sat_healthy(state::GNSSDecoderState{<:BeiDouB2bData})
    state.data.HS === Int64(0)
end
