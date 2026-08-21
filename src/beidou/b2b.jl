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
# 40 (BGTO, one midi almanac, five reduced almanacs). Others are reserved.

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
  - `SOW::Int64`: Seconds of week of the current frame's leading edge (s). Broadcast as a 20-bit count with LSB 1 s (Table 7-2).

# Message type 10 — ephemeris (Figures 6-3, 6-6, 6-7) and integrity flags

  - `sat_type::Int64`: Satellite orbit type (binary 01 = GEO, 10 = IGSO, 11 = MEO; Table 7-6).
  - `t_0e::Int64`: Ephemeris reference time (s, LSB 300).
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
  - `dif,sif,aif::Bool`: Data / signal / accuracy integrity flags for the B2b_I signal (`false` = OK, Table 7-21).
  - `sismai::Int64`: Signal in space monitoring accuracy index (4 bits; semantics deferred to a future ICD update, §7.16).

# Message type 30 — clock, group delay, ionosphere, UTC, EOP, accuracy, health

  - `WN::Int64`: BDT week number at the current frame's epoch (Table 7-2).
  - `t_0c::Int64`: Clock correction reference time (s, LSB 300).
  - `a_f0::Float64`: Clock bias (s). `a_f1::Float64`: drift (s/s). `a_f2::Float64`: drift rate (s/s²). (Table 7-3; the ICD names them a0/a1/a2.)
  - `T_GD_B2bI::Float64`: Group delay differential of the B2b_I signal relative to B3I (s, Table 7-4).
  - `α_1 … α_9::Float64`: BDGIM ionospheric model parameters (TECu, Table 7-8; α_5 is broadcast with the negative scale factor -2⁻³ already applied).
  - `A0_UTC,A1_UTC,A2_UTC::Float64`: BDT-UTC polynomial (s, s/s, s/s²; Table 7-18).
  - `Δt_LS::Int64`: Current or past leap-second count (s). `Δt_LSF::Int64`: current or future leap-second count (s).
  - `t_ot::Int64`: UTC reference time of week (s, LSB 2⁴). `WN_ot::Int64`: UTC reference week.
  - `WN_LSF::Int64`: Leap-second reference week. `DN::Int64`: leap-second reference day (0-6).
  - `t_EOP::Int64`: EOP reference time (s, LSB 2⁴; Table 7-16).
  - `PM_X,PM_Y::Float64`: Polar motion (arc-seconds). `PM_X_dot,PM_Y_dot::Float64`: drift (arc-seconds/day).
  - `ΔUT1::Float64`: UT1-UTC difference (s). `ΔUT1_dot::Float64`: its rate (s/day).
  - `sisai_t_op,sisai_ocb,sisai_oc1,sisai_oc2,sisai_oe::Int64`: Signal-in-space accuracy index fields, raw broadcast values (11/5/3/3/5 bits; semantics deferred to a future ICD update, §7.15).
  - `hs::Int64`: Satellite health status (0 = healthy, 1 = unhealthy or in test, 2-3 reserved; Table 7-20).

# Message type 40 — BGTO and almanacs (Figure 6-5)

  - `GNSS_ID::Int64`: BGTO GNSS identification (0 = not available, 1 = GPS, 2 = Galileo, 3 = GLONASS; §7.12).
  - `WN_0BGTO::Int64`, `t_0BGTO::Int64`: BGTO reference week / time of week (s, LSB 2⁴).
  - `A0_BGTO,A1_BGTO,A2_BGTO::Float64`: BDT-GNSS time offset polynomial (s, s/s, s/s²; Table 7-19).
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
    dif::Union{Nothing,Bool} = nothing
    sif::Union{Nothing,Bool} = nothing
    aif::Union{Nothing,Bool} = nothing
    sismai::Union{Nothing,Int64} = nothing

    # Message type 30: clock + TGD + ionosphere + UTC + EOP + SISAI + health
    WN::Union{Nothing,Int64} = nothing
    t_0c::Union{Nothing,Int64} = nothing
    a_f0::Union{Nothing,Float64} = nothing
    a_f1::Union{Nothing,Float64} = nothing
    a_f2::Union{Nothing,Float64} = nothing
    T_GD_B2bI::Union{Nothing,Float64} = nothing
    α_1::Union{Nothing,Float64} = nothing
    α_2::Union{Nothing,Float64} = nothing
    α_3::Union{Nothing,Float64} = nothing
    α_4::Union{Nothing,Float64} = nothing
    α_5::Union{Nothing,Float64} = nothing
    α_6::Union{Nothing,Float64} = nothing
    α_7::Union{Nothing,Float64} = nothing
    α_8::Union{Nothing,Float64} = nothing
    α_9::Union{Nothing,Float64} = nothing
    A0_UTC::Union{Nothing,Float64} = nothing
    A1_UTC::Union{Nothing,Float64} = nothing
    A2_UTC::Union{Nothing,Float64} = nothing
    Δt_LS::Union{Nothing,Int64} = nothing
    t_ot::Union{Nothing,Int64} = nothing
    WN_ot::Union{Nothing,Int64} = nothing
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
    sisai_t_op::Union{Nothing,Int64} = nothing
    sisai_ocb::Union{Nothing,Int64} = nothing
    sisai_oc1::Union{Nothing,Int64} = nothing
    sisai_oc2::Union{Nothing,Int64} = nothing
    sisai_oe::Union{Nothing,Int64} = nothing
    hs::Union{Nothing,Int64} = nothing

    # Message type 40: BGTO + almanacs
    GNSS_ID::Union{Nothing,Int64} = nothing
    WN_0BGTO::Union{Nothing,Int64} = nothing
    t_0BGTO::Union{Nothing,Int64} = nothing
    A0_BGTO::Union{Nothing,Float64} = nothing
    A1_BGTO::Union{Nothing,Float64} = nothing
    A2_BGTO::Union{Nothing,Float64} = nothing
    midi_almanacs::Union{Nothing,Dictionary{Int,BeiDouMidiAlmanac}} = nothing
    reduced_almanacs::Union{Nothing,Dictionary{Int,BeiDouReducedAlmanac}} = nothing
    WN_a::Union{Nothing,Int64} = nothing
    t_0a::Union{Nothing,Int64} = nothing
end

# Field-by-field equality: the almanac `Dictionary` fields otherwise compare
# by identity through the default struct `==`.
function Base.:(==)(a::BeiDouB2bData, b::BeiDouB2bData)
    all(getfield(a, f) == getfield(b, f) for f in fieldnames(BeiDouB2bData))
end

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
    dif = data.dif,
    sif = data.sif,
    aif = data.aif,
    sismai = data.sismai,
    WN = data.WN,
    t_0c = data.t_0c,
    a_f0 = data.a_f0,
    a_f1 = data.a_f1,
    a_f2 = data.a_f2,
    T_GD_B2bI = data.T_GD_B2bI,
    α_1 = data.α_1,
    α_2 = data.α_2,
    α_3 = data.α_3,
    α_4 = data.α_4,
    α_5 = data.α_5,
    α_6 = data.α_6,
    α_7 = data.α_7,
    α_8 = data.α_8,
    α_9 = data.α_9,
    A0_UTC = data.A0_UTC,
    A1_UTC = data.A1_UTC,
    A2_UTC = data.A2_UTC,
    Δt_LS = data.Δt_LS,
    t_ot = data.t_ot,
    WN_ot = data.WN_ot,
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
    sisai_t_op = data.sisai_t_op,
    sisai_ocb = data.sisai_ocb,
    sisai_oc1 = data.sisai_oc1,
    sisai_oc2 = data.sisai_oc2,
    sisai_oe = data.sisai_oe,
    hs = data.hs,
    GNSS_ID = data.GNSS_ID,
    WN_0BGTO = data.WN_0BGTO,
    t_0BGTO = data.t_0BGTO,
    A0_BGTO = data.A0_BGTO,
    A1_BGTO = data.A1_BGTO,
    A2_BGTO = data.A2_BGTO,
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
        dif,
        sif,
        aif,
        sismai,
        WN,
        t_0c,
        a_f0,
        a_f1,
        a_f2,
        T_GD_B2bI,
        α_1,
        α_2,
        α_3,
        α_4,
        α_5,
        α_6,
        α_7,
        α_8,
        α_9,
        A0_UTC,
        A1_UTC,
        A2_UTC,
        Δt_LS,
        t_ot,
        WN_ot,
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
        sisai_t_op,
        sisai_ocb,
        sisai_oc1,
        sisai_oc2,
        sisai_oe,
        hs,
        GNSS_ID,
        WN_0BGTO,
        t_0BGTO,
        A0_BGTO,
        A1_BGTO,
        A2_BGTO,
        midi_almanacs,
        reduced_almanacs,
        WN_a,
        t_0a,
    )
end

# ---- Cache -------------------------------------------------------------------

# Path to the committed LDPC `.alist` parity matrix. `pkgdir`-free: walk up
# from this file (src/beidou/) to the package root, then into data/.
_b2b_data_path(name) = joinpath(@__DIR__, "..", "..", "data", name)

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
    Aff3ct LDPC BP decoder for the B-CNAV3 binary image (K=486, N=972)
    """
    ldpc_decoder::LDPCBPDecoder
end

BeiDouB2bCache() = BeiDouB2bCache(
    CircularDeque{Float32}(B2B_WINDOW_SYMBOLS),
    Vector{Float32}(undef, B2B_ENCODED_SYMBOLS),
    load_ldpc_decoder(_b2b_data_path("bcnv3.alist")),
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
    # Gate 1, read straight from the deque. The shared `pack_buffer` would
    # shift a 1016-bit integer once per symbol — ~16 M word-operations per
    # second per satellite at B2b's 1000 sps — to expose 38 bits that sit at
    # three known offsets. Slice those instead; the packed window is never
    # needed, since the payload is consumed as soft symbols below.
    pre = UInt64(state.constants.preamble)
    mask = (UInt64(1) << B2B_PREAMBLE_SYMBOLS) - UInt64(1)
    head = pack_soft_bits(deque, 1, B2B_PREAMBLE_SYMBOLS)
    tail = pack_soft_bits(deque, B2B_FRAME_SYMBOLS + 1, B2B_PREAMBLE_SYMBOLS)
    inverted = pre ⊻ mask
    # Both ends must carry the preamble, and in a common polarity.
    polarity_flipped = head == inverted && tail == inverted
    (head == pre && tail == pre) || polarity_flipped || return nothing
    # Gate 2: the unencoded PRN field (window bits 17-22, ICD Figure 6-1).
    prn_bits = pack_soft_bits(deque, B2B_PREAMBLE_SYMBOLS + 1, B2B_PRN_SYMBOLS)
    prn_mask = (UInt64(1) << B2B_PRN_SYMBOLS) - UInt64(1)
    prn = Int(polarity_flipped ? prn_bits ⊻ prn_mask : prn_bits)
    prn == state.prn || return nothing
    # Gate 3: LDPC + CRC on the soft symbols (deque indices 29 .. 1000).
    sign = polarity_flipped ? -1.0f0 : 1.0f0
    window = state.cache.ldpc_window
    @inbounds for i = 1:B2B_ENCODED_SYMBOLS
        window[i] = sign * deque[B2B_ENCODED_OFFSET+i]
    end
    word = ldpc_decode_word(state.cache.ldpc_decoder, window, B2B_MESSAGE_BITS, UInt512)
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
    # syncro_sequence_length + preamble_length symbols before the newest
    # buffered symbol. Re-arming here rather than in `validate_data` (the same
    # split `b2a.jl` uses) means a time-only consumer gets a usable epoch from
    # the very first frame, instead of waiting for MT10 *and* MT30 to complete
    # the positioning set — and keeps the counter honest when a later frame
    # fails CRC, since it then simply keeps counting from the last frame that
    # actually delivered a SOW.
    GNSSDecoderState(
        state;
        raw_data = raw,
        num_bits_after_valid_syncro_sequence = state.constants.syncro_sequence_length +
                                               state.constants.preamble_length,
    )
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
        dif = get_bit(word, word_length, 456),
        sif = get_bit(word, word_length, 457),
        aif = get_bit(word, word_length, 458),
        sismai = Int64(get_bits(word, word_length, 459, 4)),
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
        t_0c = Int64(get_bits(word, word_length, 44, 11)) * 300,
        a_f0 = get_twos_complement_num(word, word_length, 55, 25) * 2.0^-34,
        a_f1 = get_twos_complement_num(word, word_length, 80, 22) * 2.0^-50,
        a_f2 = get_twos_complement_num(word, word_length, 102, 11) * 2.0^-66,
        T_GD_B2bI = get_twos_complement_num(word, word_length, 113, 12) * 2.0^-34,
        # BDGIM ionosphere (bits 125-198, Figure 6-10, Table 7-8): α₁, α₃, α₄
        # unsigned; α₂, α₆..α₉ two's complement; α₅ unsigned with the negative
        # scale factor -2⁻³.
        α_1 = Float64(get_bits(word, word_length, 125, 10)) * 2.0^-3,
        α_2 = get_twos_complement_num(word, word_length, 135, 8) * 2.0^-3,
        α_3 = Float64(get_bits(word, word_length, 143, 8)) * 2.0^-3,
        α_4 = Float64(get_bits(word, word_length, 151, 8)) * 2.0^-3,
        α_5 = Float64(get_bits(word, word_length, 159, 8)) * -(2.0^-3),
        α_6 = get_twos_complement_num(word, word_length, 167, 8) * 2.0^-3,
        α_7 = get_twos_complement_num(word, word_length, 175, 8) * 2.0^-3,
        α_8 = get_twos_complement_num(word, word_length, 183, 8) * 2.0^-3,
        α_9 = get_twos_complement_num(word, word_length, 191, 8) * 2.0^-3,
        # BDT-UTC (bits 199-295, Figure 6-11, Table 7-18):
        A0_UTC = get_twos_complement_num(word, word_length, 199, 16) * 2.0^-35,
        A1_UTC = get_twos_complement_num(word, word_length, 215, 13) * 2.0^-51,
        A2_UTC = get_twos_complement_num(word, word_length, 228, 7) * 2.0^-68,
        Δt_LS = Int64(get_twos_complement_num(word, word_length, 235, 8)),
        t_ot = Int64(get_bits(word, word_length, 243, 16)) * 16,
        WN_ot = Int64(get_bits(word, word_length, 259, 13)),
        WN_LSF = Int64(get_bits(word, word_length, 272, 13)),
        DN = Int64(get_bits(word, word_length, 285, 3)),
        Δt_LSF = Int64(get_twos_complement_num(word, word_length, 288, 8)),
        # EOP (bits 296-433, Figure 6-13, Table 7-16):
        t_EOP = Int64(get_bits(word, word_length, 296, 16)) * 16,
        PM_X = get_twos_complement_num(word, word_length, 312, 21) * 2.0^-20,
        PM_X_dot = get_twos_complement_num(word, word_length, 333, 15) * 2.0^-21,
        PM_Y = get_twos_complement_num(word, word_length, 348, 21) * 2.0^-20,
        PM_Y_dot = get_twos_complement_num(word, word_length, 369, 15) * 2.0^-21,
        ΔUT1 = get_twos_complement_num(word, word_length, 384, 31) * 2.0^-24,
        ΔUT1_dot = get_twos_complement_num(word, word_length, 415, 19) * 2.0^-25,
        # SISAI (bits 434-460, Figure 6-9) — raw broadcast values; their
        # semantics are deferred to a future ICD update (§7.15).
        sisai_t_op = Int64(get_bits(word, word_length, 434, 11)),
        sisai_ocb = Int64(get_bits(word, word_length, 445, 5)),
        sisai_oc1 = Int64(get_bits(word, word_length, 450, 3)),
        sisai_oc2 = Int64(get_bits(word, word_length, 453, 3)),
        sisai_oe = Int64(get_bits(word, word_length, 456, 5)),
        hs = Int64(get_bits(word, word_length, 461, 2)),
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
        GNSS_ID = Int64(get_bits(word, word_length, 27, 3)),
        WN_0BGTO = Int64(get_bits(word, word_length, 30, 13)),
        t_0BGTO = Int64(get_bits(word, word_length, 43, 16)) * 16,
        A0_BGTO = get_twos_complement_num(word, word_length, 59, 16) * 2.0^-35,
        A1_BGTO = get_twos_complement_num(word, word_length, 75, 13) * 2.0^-51,
        A2_BGTO = get_twos_complement_num(word, word_length, 88, 7) * 2.0^-68,
        # Almanac reference time for the five reduced almanacs (Table 7-15):
        WN_a = WN_a_reduced,
        t_0a = t_0a_reduced,
    )
    # Midi almanac (bits 95-250, Figure 6-15, Table 7-11). PRN_a = 0 marks an
    # empty block (no almanac broadcast in this frame).
    PRN_a = Int(get_bits(word, word_length, 95, 6))
    if PRN_a != 0
        alm = BeiDouMidiAlmanac(;
            PRN_a,
            sat_type = Int(get_bits(word, word_length, 101, 2)),
            WN_a = Int(get_bits(word, word_length, 103, 13)),
            t_oa = Int(get_bits(word, word_length, 116, 8)) * 4096,
            e = Float64(get_bits(word, word_length, 124, 11)) * 2.0^-16,
            δi = get_twos_complement_num(word, word_length, 135, 11) * 2.0^-14 * PI,
            sqrt_A = Float64(get_bits(word, word_length, 146, 17)) * 2.0^-4,
            Ω_0 = get_twos_complement_num(word, word_length, 163, 16) * 2.0^-15 * PI,
            Ω_dot = get_twos_complement_num(word, word_length, 179, 11) * 2.0^-33 * PI,
            ω = get_twos_complement_num(word, word_length, 190, 16) * 2.0^-15 * PI,
            M_0 = get_twos_complement_num(word, word_length, 206, 16) * 2.0^-15 * PI,
            a_f0 = get_twos_complement_num(word, word_length, 222, 11) * 2.0^-20,
            a_f1 = get_twos_complement_num(word, word_length, 233, 10) * 2.0^-37,
            health = Int(get_bits(word, word_length, 243, 8)),
        )
        raw =
            BeiDouB2bData(raw; midi_almanacs = _merge_keyed(raw.midi_almanacs, PRN_a, alm))
    end
    # Five 38-bit reduced almanacs (bits 272-461, Figure 6-12, Table 7-14).
    reduced = raw.reduced_almanacs
    for j = 0:4
        base = 272 + 38j
        PRN_a = Int(get_bits(word, word_length, base, 6))
        PRN_a == 0 && continue  # empty block
        packet = BeiDouReducedAlmanac(;
            PRN_a,
            sat_type = Int(get_bits(word, word_length, base + 6, 2)),
            WN_a = Int(WN_a_reduced),
            t_oa = Int(t_0a_reduced),
            δA = Float64(get_twos_complement_num(word, word_length, base + 8, 8)) * 512.0,
            Ω_0 = get_twos_complement_num(word, word_length, base + 16, 7) * 2.0^-6 * PI,
            Φ_0 = get_twos_complement_num(word, word_length, base + 23, 7) * 2.0^-6 * PI,
            health = Int(get_bits(word, word_length, base + 30, 8)),
        )
        reduced = _merge_keyed(reduced, PRN_a, packet)
    end
    BeiDouB2bData(raw; reduced_almanacs = reduced)
end

# ---- Validation --------------------------------------------------------------

function is_ephemeris_decoded(data::BeiDouB2bData)
    # Message type 10 carries the complete ephemeris in one CRC-gated frame.
    !isnothing(data.t_0e) &&
        !isnothing(data.sat_type) &&
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
    # Message type 30. T_GD_B2bI arrives in the same frame as the clock but is
    # a metre-level inter-signal correction: apply it downstream when present,
    # never block a fix on it (mirrors the GPS CNAV rationale).
    !isnothing(data.t_0c) &&
        !isnothing(data.a_f0) &&
        !isnothing(data.a_f1) &&
        !isnothing(data.a_f2)
end

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
        # The counter was already armed by `decode_syncro_sequence` for this
        # very frame, so promotion only publishes the data.
        return GNSSDecoderState(state; data = state.raw_data)
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

  - `prn::Int`: Pseudo-Random Noise code identifier (1-63 for BeiDou satellites)

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
    state.data.hs === Int64(0)
end
