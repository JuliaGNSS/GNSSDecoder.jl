# BeiDou B2a (B-CNAV2) navigation message decoder — BDS-SIS-ICD-B2a-1.0 (2017-12).
#
# The B-CNAV2 message rides on the B2a *data* component (`BeiDouB2aI` in GNSSSignals; the
# quadrature B2aQ component is a dataless pilot). Framing (ICD §6.2.1, Figure 6-1):
#
#   frame = 600 symbols / 3 s at 200 sps
#         = 24-symbol preamble 0xE24DE8 (unencoded, MSB first)
#         + 576 symbols of 64-ary LDPC(96,48)-encoded message
#
# The 288-bit message before encoding is PRN(6) MesType(6) SOW(18) data(234) CRC-24Q(24);
# everything but the CRC participates in the CRC (ICD §6.2.1, §6.1.2). The LDPC code is
# decoded through its exact binary image (`data/bcnv2.alist`, see
# `scripts/generate_beidou_alist.jl`) with the shared `load_ldpc_decoder` /
# `ldpc_decode_word` pipeline (`src/coding/ldpc.jl`). Per-frame flow: preamble sync (both
# ends of the 600-symbol window, either polarity) → LDPC BP decode → CRC-24Q gate → PRN
# check → per-message-type field extraction (ICD §6.2.3 Figures 6-3..6-20, §7 tables).
#
# Eight message types are defined in ICD v1.0 (10, 11, 30, 31, 32, 33, 34, 40; §6.2.3,
# Table 7-1), all parsed here. MT10 (ephemeris I) and MT11 (ephemeris II) form one
# ephemeris set identified by the MT10 IODE; the ICD requires them to be broadcast
# continuously together (§6.2.3), which pairs an MT11 with the IODE of the neighbouring
# MT10 (see `is_ephemeris_decoded`). The clock set (MT30-34) is identified by IODC, and
# §7.4.3 defines a matched pair as IODE == the 8 LSBs of IODC.

# ---- Frame constants (ICD §6.2.1) -------------------------------------------

"""
One B-CNAV2 frame: 24 preamble symbols + 576 encoded symbols = 600 symbols (3 s at 200 sps).
"""
const B2A_FRAME_SYMBOLS = 600
const B2A_PREAMBLE_SYMBOLS = 24
const B2A_WINDOW_SYMBOLS = B2A_FRAME_SYMBOLS + B2A_PREAMBLE_SYMBOLS  # 624
const B2A_ENCODED_SYMBOLS = B2A_FRAME_SYMBOLS - B2A_PREAMBLE_SYMBOLS  # 576

"""
One B-CNAV2 message before encoding: 288 bits (PRN 6, MesType 6, SOW 18, data 234, CRC 24).
"""
const B2A_MESSAGE_BITS = 288

"""
B-CNAV2 preamble `0xE24DE8` = `111000100100110111101000`, MSB transmitted first (ICD §6.2.1).
"""
const B2A_PREAMBLE = UInt64(0xE24DE8)

"""
$(TYPEDEF)

Constants for the BeiDou B2a (B-CNAV2) decoder (BDS-SIS-ICD-B2a-1.0).

# Fields

$(TYPEDFIELDS)
"""
Base.@kwdef struct BeiDouB2aConstants <: AbstractGNSSConstants
    """
    Frame length drained after each decoded frame (600 symbols)
    """
    sync_unit_length::Int = B2A_FRAME_SYMBOLS
    """
    Preamble length (24 symbols, `0xE24DE8`, ICD §6.2.1)
    """
    preamble_length::Int = B2A_PREAMBLE_SYMBOLS
    """
    Preamble bit pattern (MSB-first packing of `111000100100110111101000`)
    """
    preamble::UInt64 = B2A_PREAMBLE
end

# ---- Data container -----------------------------------------------------------

"""
$(TYPEDEF)

Decoded BeiDou B2a B-CNAV2 navigation data (BDS-SIS-ICD-B2a-1.0 §6.2.3, §7).

Every field is `nothing` until its carrying message type has been decoded.
Angles are stored in radians (ICD semicircle values scaled by π =
`GNSS_PI`), times in seconds, week numbers in weeks of BDT (epoch
2006-01-01T00:00:00 UTC, §7.3).

Ephemeris pairing: ephemeris I is in message type 10 (with the set's IODE), ephemeris II
in message type 11, which carries *no* IOD; the ICD requires MT10 and MT11 to be broadcast
continuously together (§6.2.3). `SOW_mt10` / `SOW_mt11` record the SOW of the frames that
delivered each half so `is_ephemeris_decoded` can require adjacent frames (|ΔSOW| = 3 s).

# Fields

$(TYPEDFIELDS)
"""
Base.@kwdef struct BeiDouCNAV2Data <: AbstractBeiDouCNAVData
    """
    Message type of the most recently decoded frame (0 = none yet)
    """
    last_message_type::Int = 0
    """
    Seconds of week of the most recent frame (s; epoch = rising edge of that frame's first preamble chip, §7.3)
    """
    SOW::Union{Nothing,Int64} = nothing
    """
    BDT week number (MT10, §7.3)
    """
    WN::Union{Nothing,Int64} = nothing

    """
    Satellite health status HS (2 bits: 0 healthy, 1 unhealthy or in test, 2-3 reserved; MT11/30-34/40, §7.14)
    """
    HS::Union{Nothing,Int64} = nothing
    """
    B2a data integrity flag (0 = message error within predictive accuracy, §7.15)
    """
    DIF_B2a::Union{Nothing,Bool} = nothing
    """
    B2a signal integrity flag (0 = signal normal, §7.15)
    """
    SIF_B2a::Union{Nothing,Bool} = nothing
    """
    B2a accuracy integrity flag (0 = SISMAI value valid, §7.15)
    """
    AIF_B2a::Union{Nothing,Bool} = nothing
    """
    Signal in space monitoring accuracy index (4 bits; definition deferred to a future ICD update, §7.17)
    """
    SISMAI::Union{Nothing,Int64} = nothing
    """
    B1C data integrity flag (also broadcast on B-CNAV2, §7.15)
    """
    DIF_B1C::Union{Nothing,Bool} = nothing
    """
    B1C signal integrity flag (§7.15)
    """
    SIF_B1C::Union{Nothing,Bool} = nothing
    """
    B1C accuracy integrity flag (§7.15)
    """
    AIF_B1C::Union{Nothing,Bool} = nothing

    # ---- Ephemeris (message types 10 + 11, §7.7) ----
    """
    Issue of data, ephemeris (8 bits, MT10, §7.4.1)
    """
    IODE::Union{Nothing,Int64} = nothing
    """
    Ephemeris reference time (s, ×300; the ICD writes `toe`)
    """
    t_0e::Union{Nothing,Int64} = nothing
    """
    Satellite orbit type (2 bits: 1 = GEO, 2 = IGSO, 3 = MEO, 0 reserved; Table 7-8)
    """
    sat_type::Union{Nothing,Int64} = nothing
    """
    Semi-major axis difference at reference time (m; vs A_ref = 27906100 m MEO / 42162200 m IGSO-GEO)
    """
    ΔA::Union{Nothing,Float64} = nothing
    """
    Change rate of semi-major axis (m/s)
    """
    A_dot::Union{Nothing,Float64} = nothing
    """
    Mean motion difference at reference time (rad/s)
    """
    Δn_0::Union{Nothing,Float64} = nothing
    """
    Rate of mean motion difference (rad/s²)
    """
    Δn_0_dot::Union{Nothing,Float64} = nothing
    """
    Mean anomaly at reference time (rad)
    """
    M_0::Union{Nothing,Float64} = nothing
    """
    Eccentricity (dimensionless)
    """
    e::Union{Nothing,Float64} = nothing
    """
    Argument of perigee (rad)
    """
    ω::Union{Nothing,Float64} = nothing
    """
    SOW of the frame that delivered ephemeris I (pairing bookkeeping, see type docstring)
    """
    SOW_mt10::Union{Nothing,Int64} = nothing

    """
    Longitude of ascending node at weekly epoch (rad, MT11)
    """
    Ω_0::Union{Nothing,Float64} = nothing
    """
    Inclination angle at reference time (rad)
    """
    i_0::Union{Nothing,Float64} = nothing
    """
    Rate of right ascension (rad/s)
    """
    Ω_dot::Union{Nothing,Float64} = nothing
    """
    Rate of inclination angle (rad/s)
    """
    i_dot::Union{Nothing,Float64} = nothing
    """
    Amplitude of sine harmonic correction to inclination (rad)
    """
    C_is::Union{Nothing,Float64} = nothing
    """
    Amplitude of cosine harmonic correction to inclination (rad)
    """
    C_ic::Union{Nothing,Float64} = nothing
    """
    Amplitude of sine harmonic correction to orbit radius (m)
    """
    C_rs::Union{Nothing,Float64} = nothing
    """
    Amplitude of cosine harmonic correction to orbit radius (m)
    """
    C_rc::Union{Nothing,Float64} = nothing
    """
    Amplitude of sine harmonic correction to argument of latitude (rad)
    """
    C_us::Union{Nothing,Float64} = nothing
    """
    Amplitude of cosine harmonic correction to argument of latitude (rad)
    """
    C_uc::Union{Nothing,Float64} = nothing
    """
    SOW of the frame that delivered ephemeris II (pairing bookkeeping, see type docstring)
    """
    SOW_mt11::Union{Nothing,Int64} = nothing

    # ---- Clock correction (message types 30-34, §7.5) ----
    """
    Issue of data, clock (10 bits; a matched pair has IODE == IODC & 0xFF, §7.4.3)
    """
    IODC::Union{Nothing,Int64} = nothing
    """
    Clock correction reference time (s, ×300; the ICD writes `toc`)
    """
    t_0c::Union{Nothing,Int64} = nothing
    """
    SV clock bias (s; the ICD names it `a_0`)
    """
    a_f0::Union{Nothing,Float64} = nothing
    """
    SV clock drift (s/s; the ICD names it `a_1`)
    """
    a_f1::Union{Nothing,Float64} = nothing
    """
    SV clock drift rate (s/s²; the ICD names it `a_2`)
    """
    a_f2::Union{Nothing,Float64} = nothing

    # ---- Group delay (message type 30, §7.6) ----
    """
    Group delay differential of the B2a pilot component vs the B3I-referenced clock (s)
    """
    T_GD_B2ap::Union{Nothing,Float64} = nothing
    """
    Group delay differential between the B2a data and pilot components (s)
    """
    ISC_B2ad::Union{Nothing,Float64} = nothing
    """
    Group delay differential of the B1C pilot component (s)
    """
    T_GD_B1Cp::Union{Nothing,Float64} = nothing

    # ---- BDGIM ionosphere (message type 30, §7.8) ----
    """
    BDGIM parameter α₁ (TECu). The ICD names the nine `α_1` … `α_9`; the `bdgim`
    qualifier is this package's, keeping them distinct from the Klobuchar
    `α_0` … `α_3` of the other containers (see `b1c.jl`).
    """
    α_bdgim_1::Union{Nothing,Float64} = nothing
    """
    BDGIM parameter α₂ (TECu)
    """
    α_bdgim_2::Union{Nothing,Float64} = nothing
    """
    BDGIM parameter α₃ (TECu)
    """
    α_bdgim_3::Union{Nothing,Float64} = nothing
    """
    BDGIM parameter α₄ (TECu)
    """
    α_bdgim_4::Union{Nothing,Float64} = nothing
    """
    BDGIM parameter α₅ (TECu; broadcast with scale −2⁻³, Table 7-10)
    """
    α_bdgim_5::Union{Nothing,Float64} = nothing
    """
    BDGIM parameter α₆ (TECu)
    """
    α_bdgim_6::Union{Nothing,Float64} = nothing
    """
    BDGIM parameter α₇ (TECu)
    """
    α_bdgim_7::Union{Nothing,Float64} = nothing
    """
    BDGIM parameter α₈ (TECu)
    """
    α_bdgim_8::Union{Nothing,Float64} = nothing
    """
    BDGIM parameter α₉ (TECu)
    """
    α_bdgim_9::Union{Nothing,Float64} = nothing

    # ---- Earth orientation (message type 32, §7.11) ----
    """
    EOP data reference time (s, ×2⁴)
    """
    t_EOP::Union{Nothing,Int64} = nothing
    """
    X-axis polar motion at reference time (arc-seconds)
    """
    PM_X::Union{Nothing,Float64} = nothing
    """
    X-axis polar motion drift (arc-seconds/day)
    """
    PM_X_dot::Union{Nothing,Float64} = nothing
    """
    Y-axis polar motion at reference time (arc-seconds)
    """
    PM_Y::Union{Nothing,Float64} = nothing
    """
    Y-axis polar motion drift (arc-seconds/day)
    """
    PM_Y_dot::Union{Nothing,Float64} = nothing
    """
    UT1−UTC difference at reference time (s)
    """
    ΔUT1::Union{Nothing,Float64} = nothing
    """
    Rate of UT1−UTC difference (s/day)
    """
    ΔUT1_dot::Union{Nothing,Float64} = nothing

    # ---- BDT-UTC offset (message type 34, §7.12) ----
    """
    Bias coefficient of BDT relative to UTC (s)
    """
    A_0UTC::Union{Nothing,Float64} = nothing
    """
    Drift coefficient of BDT relative to UTC (s/s)
    """
    A_1UTC::Union{Nothing,Float64} = nothing
    """
    Drift rate coefficient of BDT relative to UTC (s/s²)
    """
    A_2UTC::Union{Nothing,Float64} = nothing
    """
    Current or past leap second count (s)
    """
    Δt_LS::Union{Nothing,Int64} = nothing
    """
    Reference time of week for the UTC parameters (s, ×2⁴; the ICD writes `t_ot`)
    """
    t_0t::Union{Nothing,Int64} = nothing
    """
    Reference week number for the UTC parameters (the ICD writes `WN_ot`)
    """
    WN_0t::Union{Nothing,Int64} = nothing
    """
    Leap second reference week number
    """
    WN_LSF::Union{Nothing,Int64} = nothing
    """
    Leap second reference day number (0-6)
    """
    DN::Union{Nothing,Int64} = nothing
    """
    Current or future leap second count (s)
    """
    Δt_LSF::Union{Nothing,Int64} = nothing

    # ---- BDT-GNSS offset (message type 33, §7.13) ----
    """
    GNSS type the BGTO parameters refer to (0 = unavailable, 1 = GPS, 2 = Galileo,
    3 = GLONASS, 4-7 reserved; §7.13.1)
    """
    GNSS_ID::Union{Nothing,Int64} = nothing
    """
    BGTO reference week number
    """
    WN_0BGTO::Union{Nothing,Int64} = nothing
    """
    BGTO reference time of week (s, ×2⁴)
    """
    t_0BGTO::Union{Nothing,Int64} = nothing
    """
    Bias coefficient of BDT relative to the identified GNSS time (s)
    """
    A_0BGTO::Union{Nothing,Float64} = nothing
    """
    Drift coefficient of BDT relative to the identified GNSS time (s/s)
    """
    A_1BGTO::Union{Nothing,Float64} = nothing
    """
    Drift rate coefficient of BDT relative to the identified GNSS time (s/s²)
    """
    A_2BGTO::Union{Nothing,Float64} = nothing

    # ---- Signal in space accuracy (message types 34/40, §7.16) ----
    # ICD v1.0 defers the SISAI definitions to a future update, so the raw
    # broadcast integers are stored unscaled.
    """
    Time of week for data prediction (s; 11-bit count, LSB 300 s — see
    `beidou_sisai_oc_block`)
    """
    t_op::Union{Nothing,Int64} = nothing
    """
    Satellite orbit radius & fixed clock bias accuracy index (raw, §7.16)
    """
    SISAI_ocb::Union{Nothing,Int64} = nothing
    """
    Satellite clock bias accuracy index (raw, §7.16)
    """
    SISAI_oc1::Union{Nothing,Int64} = nothing
    """
    Satellite clock drift accuracy index (raw, §7.16)
    """
    SISAI_oc2::Union{Nothing,Int64} = nothing
    """
    Satellite orbit along-track/cross-track accuracy index (raw, MT40, §7.16)
    """
    SISAI_oe::Union{Nothing,Int64} = nothing

    # ---- Almanacs (message types 31/33/40, §7.9 / §7.10) ----
    """
    Reduced almanacs keyed by PRN 1-63 (MT31: three per message; MT33: one per
    message); preallocated, overwritten in place by [`decode!`](@ref)
    """
    reduced_almanacs::Union{Nothing,SlotDictionary{BeiDouReducedAlmanac,64}} = nothing
    """
    Midi almanacs keyed by PRN 1-63 (MT40: one per message); preallocated,
    overwritten in place by [`decode!`](@ref)
    """
    midi_almanacs::Union{Nothing,SlotDictionary{BeiDouMidiAlmanac,64}} = nothing
end

@inline function BeiDouCNAV2Data(
    data::BeiDouCNAV2Data;
    last_message_type = data.last_message_type,
    SOW = data.SOW,
    WN = data.WN,
    HS = data.HS,
    DIF_B2a = data.DIF_B2a,
    SIF_B2a = data.SIF_B2a,
    AIF_B2a = data.AIF_B2a,
    SISMAI = data.SISMAI,
    DIF_B1C = data.DIF_B1C,
    SIF_B1C = data.SIF_B1C,
    AIF_B1C = data.AIF_B1C,
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
    SOW_mt10 = data.SOW_mt10,
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
    SOW_mt11 = data.SOW_mt11,
    IODC = data.IODC,
    t_0c = data.t_0c,
    a_f0 = data.a_f0,
    a_f1 = data.a_f1,
    a_f2 = data.a_f2,
    T_GD_B2ap = data.T_GD_B2ap,
    ISC_B2ad = data.ISC_B2ad,
    T_GD_B1Cp = data.T_GD_B1Cp,
    α_bdgim_1 = data.α_bdgim_1,
    α_bdgim_2 = data.α_bdgim_2,
    α_bdgim_3 = data.α_bdgim_3,
    α_bdgim_4 = data.α_bdgim_4,
    α_bdgim_5 = data.α_bdgim_5,
    α_bdgim_6 = data.α_bdgim_6,
    α_bdgim_7 = data.α_bdgim_7,
    α_bdgim_8 = data.α_bdgim_8,
    α_bdgim_9 = data.α_bdgim_9,
    t_EOP = data.t_EOP,
    PM_X = data.PM_X,
    PM_X_dot = data.PM_X_dot,
    PM_Y = data.PM_Y,
    PM_Y_dot = data.PM_Y_dot,
    ΔUT1 = data.ΔUT1,
    ΔUT1_dot = data.ΔUT1_dot,
    A_0UTC = data.A_0UTC,
    A_1UTC = data.A_1UTC,
    A_2UTC = data.A_2UTC,
    Δt_LS = data.Δt_LS,
    t_0t = data.t_0t,
    WN_0t = data.WN_0t,
    WN_LSF = data.WN_LSF,
    DN = data.DN,
    Δt_LSF = data.Δt_LSF,
    GNSS_ID = data.GNSS_ID,
    WN_0BGTO = data.WN_0BGTO,
    t_0BGTO = data.t_0BGTO,
    A_0BGTO = data.A_0BGTO,
    A_1BGTO = data.A_1BGTO,
    A_2BGTO = data.A_2BGTO,
    t_op = data.t_op,
    SISAI_ocb = data.SISAI_ocb,
    SISAI_oc1 = data.SISAI_oc1,
    SISAI_oc2 = data.SISAI_oc2,
    SISAI_oe = data.SISAI_oe,
    reduced_almanacs = data.reduced_almanacs,
    midi_almanacs = data.midi_almanacs,
)
    BeiDouCNAV2Data(
        last_message_type,
        SOW,
        WN,
        HS,
        DIF_B2a,
        SIF_B2a,
        AIF_B2a,
        SISMAI,
        DIF_B1C,
        SIF_B1C,
        AIF_B1C,
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
        SOW_mt10,
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
        SOW_mt11,
        IODC,
        t_0c,
        a_f0,
        a_f1,
        a_f2,
        T_GD_B2ap,
        ISC_B2ad,
        T_GD_B1Cp,
        α_bdgim_1,
        α_bdgim_2,
        α_bdgim_3,
        α_bdgim_4,
        α_bdgim_5,
        α_bdgim_6,
        α_bdgim_7,
        α_bdgim_8,
        α_bdgim_9,
        t_EOP,
        PM_X,
        PM_X_dot,
        PM_Y,
        PM_Y_dot,
        ΔUT1,
        ΔUT1_dot,
        A_0UTC,
        A_1UTC,
        A_2UTC,
        Δt_LS,
        t_0t,
        WN_0t,
        WN_LSF,
        DN,
        Δt_LSF,
        GNSS_ID,
        WN_0BGTO,
        t_0BGTO,
        A_0BGTO,
        A_1BGTO,
        A_2BGTO,
        t_op,
        SISAI_ocb,
        SISAI_oc1,
        SISAI_oc2,
        SISAI_oe,
        reduced_almanacs,
        midi_almanacs,
    )
end

# Field-by-field equality: the almanac `SlotDictionary` fields otherwise compare
# by identity through the default struct `==`.
Base.:(==)(a::BeiDouCNAV2Data, b::BeiDouCNAV2Data) = fields_equal(a, b)

# Every container field at its ICD size: one almanac slot per PRN (1-63).
preallocated_data(::Type{BeiDouCNAV2Data}) = BeiDouCNAV2Data(;
    reduced_almanacs = SlotDictionary{BeiDouReducedAlmanac,64}(),
    midi_almanacs = SlotDictionary{BeiDouMidiAlmanac,64}(),
)

# ---- Cache --------------------------------------------------------------------

"""
$(TYPEDEF)

Per-decoder cache for BeiDou B2a: the soft-symbol deque (624 = 600 frame symbols + 24
next-frame preamble symbols) and the Aff3ct LDPC BP decoder handle for the B-CNAV2 binary
image (`data/bcnv2.alist`), reused across frames. The keyword constructor rebuilds the
cache with a new `frame_decoded`, sharing every buffer.

# Fields

$(TYPEDFIELDS)
"""
struct BeiDouB2aCache <: AbstractGNSSCache
    """
    Soft-symbol buffer (624 = 600 frame + 24 next-frame preamble)
    """
    soft_buffer::CircularDeque{Float32}
    """
    LDPC decoder and its scratch buffers for the B-CNAV2 binary image (K=288, N=576)
    """
    ldpc::LDPCScratch
    """
    576-entry LLR scratch copied out of `soft_buffer` per frame
    """
    llr_scratch::Vector{Float32}
    """
    Whether the frame just synced cleared LDPC, CRC and the PRN gate
    (see `is_b2a_frame_decoded`)
    """
    frame_decoded::Bool
    """
    Preallocated almanac stores `raw_data` and `data` are decoded into,
    overwritten in place by [`decode!`](@ref)
    """
    storage::DataStorage{BeiDouCNAV2Data}
end

function BeiDouB2aCache()
    BeiDouB2aCache(
        CircularDeque{Float32}(B2A_WINDOW_SYMBOLS),
        committed_ldpc_scratch("bcnv2.alist"),
        Vector{Float32}(undef, B2A_ENCODED_SYMBOLS),
        false,
        DataStorage{BeiDouCNAV2Data}(),
    )
end

function BeiDouB2aCache(cache::BeiDouB2aCache; frame_decoded = cache.frame_decoded)
    BeiDouB2aCache(
        cache.soft_buffer,
        cache.ldpc,
        cache.llr_scratch,
        frame_decoded,
        cache.storage,
    )
end

"""
True when the frame that just synced was decoded into `raw_data`.

`decode_sync_unit` clears the flag on entry and sets it only once the frame has cleared
LDPC, CRC-24Q and the PRN gate. Otherwise `raw_data.SOW` still names the *previous* frame,
and promoting it would re-anchor the symbol counter to the frame that just synced,
freezing the time of transmission for one frame.
"""
is_b2a_frame_decoded(state::GNSSDecoderState{<:BeiDouCNAV2Data}) = state.cache.frame_decoded

# Two B2a caches are equal when their soft buffers match: the LDPC handle is stateless
# w.r.t. equality, and `frame_decoded` only carries one sync's outcome from
# `decode_sync_unit` to `validate_data`.
# `decode_sync_unit` to `validate_data`.
function Base.:(==)(a::BeiDouB2aCache, b::BeiDouB2aCache)
    deques_equal(a.soft_buffer, b.soft_buffer)
end

# No `packed_buffer_type` method: B2a reads the two preamble windows straight off the
# soft buffer, so the packed 624-bit window the default sync path would build per symbol
# is never needed (cf. Galileo I/NAV, E5a, E6-B and BeiDou B2b).

"""
    try_sync(state::GNSSDecoderState{<:BeiDouCNAV2Data}) -> Union{Nothing,Bool}

B-CNAV2 frame sync: the preamble `0xE24DE8` must appear at both ends of the 624-symbol
window (this frame's start and the next's), both upright or both inverted. Returns the
resolved polarity or `nothing`. The LDPC, CRC-24Q and own-PRN gates run in
`decode_sync_unit`.
"""
try_sync(state::GNSSDecoderState{<:BeiDouCNAV2Data}) = find_preamble_in_deque(
    soft_buffer(state),
    state.constants.preamble,
    B2A_PREAMBLE_SYMBOLS,
    B2A_FRAME_SYMBOLS,
)

# Record the polarity `try_sync` resolved; `decode_sync_unit` reads it back
# off the state when it copies the frame out of the deque (mirrors Galileo I/NAV).
function complement_buffer_if_necessary(
    state::GNSSDecoderState{<:BeiDouCNAV2Data},
    polarity_flipped::Bool,
)
    GNSSDecoderState(state; is_shifted_by_180_degrees = polarity_flipped), polarity_flipped
end

# ---- Constructors ---------------------------------------------------------------

"""
$(TYPEDSIGNATURES)

Create a decoder state for BeiDou B2a (B-CNAV2) navigation messages.

Initializes a [`GNSSDecoderState`](@ref) for the B-CNAV2 message from the 200 sps soft
symbols of the B2a data component (BDS-SIS-ICD-B2a-1.0). Each sync attempt matches the
24-symbol preamble `0xE24DE8` at both ends of the 600-symbol frame window (either
polarity); a matched frame is LDPC-decoded (64-ary LDPC(96,48) binary image), gated on
CRC-24Q and the broadcast PRN, and dispatched to per-message-type parsers (types 10, 11,
30-34, 40). Decoded fields land in a [`BeiDouCNAV2Data`](@ref).

# Arguments

  - `prn::Int`: Pseudo-Random Noise code identifier (1-63 for BeiDou satellites)

# Returns

  - `GNSSDecoderState{BeiDouCNAV2Data}`: Initialized decoder state for BeiDou B2a

# Example

```julia
state = BeiDouB2aDecoderState(19)         # PRN 19
state = decode!(state, soft_bits, num_bits)
if is_sat_healthy(state)
    # Use state.data for positioning
end
```

# See Also

  - [`GNSSDecoderState`](@ref): The underlying state structure
  - [`decode!`](@ref): Decode soft symbols using this state
  - [`reset_decoder_state!`](@ref): Reset after signal loss
  - [`is_sat_healthy`](@ref): Check satellite health status
"""
function BeiDouB2aDecoderState(prn)
    1 <= prn <= 63 || throw(ArgumentError("BeiDou PRN must be in 1..63"))
    GNSSDecoderState(
        prn,
        BeiDouCNAV2Data(),
        BeiDouCNAV2Data(),
        BeiDouB2aConstants(),
        BeiDouB2aCache(),
        nothing,
        false,
    )
end

# Dispatch from the GNSSSignals type: B-CNAV2 rides on the in-phase data component
# `BeiDouB2aI`; `BeiDouB2aQ` is the dataless pilot (BDS-SIS-ICD-B2a-1.0 §4.2).
function GNSSDecoderState(system::BeiDouB2aI, prn)
    BeiDouB2aDecoderState(prn)
end

# The signal this decoder demodulates (the B2a data component); signal metadata is
# forwarded through it (see `src/gnss.jl`).
get_signal_type(::BeiDouB2aConstants) = BeiDouB2aI

# ---- Completeness / health -----------------------------------------------------

function is_ephemeris_decoded(data::BeiDouCNAV2Data)
    # Message type 10 (ephemeris I + IODE) ...
    !isnothing(data.IODE) &&
        !isnothing(data.t_0e) &&
        is_known_sat_type(data.sat_type) &&
        !isnothing(data.ΔA) &&
        !isnothing(data.A_dot) &&
        !isnothing(data.Δn_0) &&
        !isnothing(data.Δn_0_dot) &&
        !isnothing(data.M_0) &&
        !isnothing(data.e) &&
        !isnothing(data.ω) &&
        # ... and message type 11 (ephemeris II) ...
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
        # ... delivered by adjacent frames. MT11 carries no IOD; the ICD pairs it with MT10
        # by requiring both to be broadcast continuously together (§6.2.3), i.e. in
        # consecutive 3 s frames (either order). Adjacency is evaluated modulo the week so
        # a pair straddling the week crossover (SOW 604797 ↔ 0) still pairs up.
        _b2a_sow_adjacent(data.SOW_mt10, data.SOW_mt11)
end

# `nothing` methods: Julia 1.10 does not narrow the field reads above through
# `isnothing`, so without them the call would dispatch dynamically.
_b2a_sow_adjacent(::Nothing, ::Union{Nothing,Integer}) = false
_b2a_sow_adjacent(::Integer, ::Nothing) = false

_b2a_sow_adjacent(sow_a::Integer, sow_b::Integer) =
    abs(GNSSSignals.fold_week_crossover(sow_a - sow_b)) == 3

function is_clock_correction_decoded(data::BeiDouCNAV2Data)
    # The clock block + IODC are carried identically in message types 30-34.
    #
    # T_GD/ISC are deliberately NOT required, unlike on B1C and B2b where the single-band
    # group delay rides in the clock's own block. Here `T_GD_B2ap`, `ISC_B2ad` and
    # `T_GD_B1Cp` are in MT30 *only* while clock and IODC are in all of MT30-34, so gating
    # on them would make a fix wait for MT30, and BDS-SIS-ICD-B2a-1.0 §6.2 does not
    # schedule it: "The broadcast order of the B-CNAV2 message types may be dynamically
    # adjusted, however Message Types 10 and 11 shall be broadcast continuously
    # together." Only the MT10/11 pair has a guaranteed cadence. Same reasoning as the GPS
    # CNAV decoder's T_GD.
    !isnothing(data.IODC) &&
        !isnothing(data.t_0c) &&
        !isnothing(data.a_f0) &&
        !isnothing(data.a_f1) &&
        !isnothing(data.a_f2)
end

# B-CNAV2 broadcasts the seconds of week as an 18-bit count in 3-second units,
# scaled at parse (§7.3, Table 7-2).
get_broadcast_time_of_week(data::BeiDouCNAV2Data) = data.SOW

# B-CNAV2's week number is 13 bits from BDT's 2006 origin, absolute until 2163.
get_broadcast_week(data::BeiDouCNAV2Data; approximate_year::Integer) = data.WN

"""
$(TYPEDSIGNATURES)

The BDT-to-`target` offset from message type 33, or `nothing` when this satellite is not
currently broadcasting one for `target`. One set at a time, tagged by `GNSS_ID`; see
`beidou_bgto_offset`.
"""
time_offset_parameters(state::GNSSDecoderState{<:BeiDouCNAV2Data}, target::TimeSystem) =
    beidou_bgto_offset(state, target)

function has_positioning_messages(data::BeiDouCNAV2Data)
    !isnothing(data.SOW) &&
        is_ephemeris_decoded(data) &&
        is_clock_correction_decoded(data) &&
        # Matched ephemeris/clock pair: IODE must equal the 8 LSBs of IODC
        # (BDS-SIS-ICD-B2a-1.0 §7.4.3).
        beidou_iode_matches_iodc(data.IODE, data.IODC)
end

"""
$(TYPEDSIGNATURES)

Check if the BeiDou B2a satellite is healthy and usable for positioning.

Examines the 2-bit satellite health status HS decoded from the most recent
message of types 11 or 30-34 or 40 (BDS-SIS-ICD-B2a-1.0 §7.14, Table 7-22):
the satellite is healthy iff `HS == 0` ("the satellite is healthy / provides
services"). HS = 1 means unhealthy or in test; 2-3 are reserved and treated
as unhealthy.

!!! warning

    Requires a health-carrying message to have been decoded and the
    positioning set to have been validated; returns `false` until then.

# Arguments

  - `state::GNSSDecoderState{<:BeiDouCNAV2Data}`: BeiDou B2a decoder state.

# Returns

  - `Bool`: `true` iff the broadcast health status indicates a usable satellite.
"""
function is_sat_healthy(state::GNSSDecoderState{<:BeiDouCNAV2Data})
    !isnothing(state.data.HS) && state.data.HS == 0
end

"""
$(TYPEDSIGNATURES)

Reset a BeiDou B2a decoder state after a signal loss or reacquisition.

Clears the soft-symbol buffer and the seconds-of-week field while preserving
the remaining decoded ephemeris and clock data in `raw_data`, so a receiver
can re-use the satellite after reacquisition without re-decoding all message
types. Mirrors the semantics of the GPS CNAV implementation.

# Arguments

  - `state::GNSSDecoderState{<:BeiDouCNAV2Data}`: Current BeiDou B2a decoder state

# Returns

  - `GNSSDecoderState{<:BeiDouCNAV2Data}`: Reset decoder state with cleared buffers

# See Also

  - [`BeiDouB2aDecoderState`](@ref): Create a fresh decoder state
  - [`decode!`](@ref): Continue decoding after reset
"""
function reset_decoder_state!(state::GNSSDecoderState{<:BeiDouCNAV2Data})
    empty!(state.cache.soft_buffer)
    GNSSDecoderState(
        state;
        raw_data = BeiDouCNAV2Data(state.raw_data; SOW = nothing),
        data = BeiDouCNAV2Data(),
        cache = BeiDouB2aCache(state.cache; frame_decoded = false),
        num_bits_since_time_reference = nothing,
        is_shifted_by_180_degrees = false,
    )
end

# ---- Frame pipeline -------------------------------------------------------------

"""
    decode_sync_unit(state::GNSSDecoderState{<:BeiDouCNAV2Data}, buffer)

Process one preamble-synchronized B-CNAV2 frame: LDPC-decode the 576 encoded soft symbols
(64-ary LDPC(96,48) binary image), gate on CRC-24Q over the 288-bit message, require the
broadcast PRN (bits 1-6) to match this decoder's (a preamble match on a cross-correlating
signal must not corrupt this satellite's data), parse the header (message type at bits
7-12, SOW at bits 13-30, ICD §6.2.1), then dispatch to the per-message-type parser.
Unknown or reserved types keep the decoded header only.

The 18-bit SOW counts 3 s units (Table 7-2) and denotes the rising edge of the *current*
frame's first preamble chip, 624 symbols (one frame plus one preamble) behind the newest
buffered symbol. The symbol counter is re-armed in `validate_data`, on promotion.

Sets the cache's `frame_decoded` flag (see `is_b2a_frame_decoded`) only when the frame
clears every gate and is parsed into `raw_data`; a frame dropped by LDPC/CRC or the PRN
gate leaves it false and `raw_data` untouched.
"""
function decode_sync_unit(state::GNSSDecoderState{<:BeiDouCNAV2Data}, ::Bool)
    # Clear the flag first: every path out of this function leaves it false
    # unless the frame is parsed into `raw_data` below.
    state =
        GNSSDecoderState(state; cache = BeiDouB2aCache(state.cache; frame_decoded = false))
    # The 576 encoded symbols sit between the leading preamble and the next frame's
    # (deque indices preamble_length+1 .. sync_unit_length). Negate the LLRs when the sync
    # hook flagged the frame as inverted (180° ambiguity).
    llr = copy_soft_window!(
        state.cache.llr_scratch,
        soft_buffer(state),
        B2A_PREAMBLE_SYMBOLS,
        B2A_ENCODED_SYMBOLS,
        state.is_shifted_by_180_degrees,
    )

    word = ldpc_decode_word(state.cache.ldpc, llr, UInt320)
    isnothing(word) && return state  # silently drop on CRC failure

    word_length = B2A_MESSAGE_BITS
    # PRN gate (ICD §7.1): the message states the transmitting satellite's PRN.
    prn = Int(get_bits(word, word_length, 1, 6))
    prn == state.prn || return state

    message_id = Int(get_bits(word, word_length, 7, 6))
    # SOW: 18 bits in units of 3 s (Table 7-2), epoch at the current frame's
    # first preamble chip (§7.3).
    SOW = Int64(get_bits(word, word_length, 13, 18)) * 3
    raw = BeiDouCNAV2Data(state.raw_data; last_message_type = message_id, SOW)

    raw = if message_id == 10
        parse_b2a_mt10(raw, word)
    elseif message_id == 11
        parse_b2a_mt11(raw, word)
    elseif message_id == 30
        parse_b2a_mt30(raw, word)
    elseif message_id == 31
        parse_b2a_mt31!(raw, word, state.cache.storage.raw)
    elseif message_id == 32
        parse_b2a_mt32(raw, word)
    elseif message_id == 33
        parse_b2a_mt33!(raw, word, state.cache.storage.raw)
    elseif message_id == 34
        parse_b2a_mt34(raw, word)
    elseif message_id == 40
        parse_b2a_mt40!(raw, word, state.cache.storage.raw)
    else
        raw  # unknown/reserved message type: header only
    end

    # The symbol counter is deliberately NOT re-armed here. `get_broadcast_time_of_week`
    # reads the *validated* `data.SOW`, which only advances when `validate_data` promotes,
    # and the MT10/MT11 adjacency gate skips promotion for the frame carrying a new
    # ephemeris-I block until its partner arrives. Re-arming on every CRC-valid frame while
    # the published SOW stood still made `tow + num_bits_since_time_reference/rate` read one
    # 3 s frame low during every skip: a ~9 km along-track error (measured on TEX-CUP as a
    # 2.2 km position excursion once per message cycle). The counter instead keeps counting
    # from the last *promoted* frame, consistent through skipped promotions and CRC
    # failures; `validate_data` re-anchors it on promotion.
    GNSSDecoderState(
        state;
        raw_data = raw,
        cache = BeiDouB2aCache(state.cache; frame_decoded = true),
    )
end

"""
    validate_data(state::GNSSDecoderState{<:BeiDouCNAV2Data})

Promote `raw_data` to `data` once the minimum positioning set is decoded and consistent:
the MT10+MT11 ephemeris pair from adjacent frames, a clock set from any of MT30-34, and
IODE == the 8 LSBs of IODC (BDS-SIS-ICD-B2a-1.0 §7.4.3). Only a frame that was actually
decoded is promoted (see `is_b2a_frame_decoded`).

`publish_data` overwrites the preallocated validated almanac stores, so `data` never
shares a container with `raw_data`.
"""
function validate_data(state::GNSSDecoderState{<:BeiDouCNAV2Data})
    # `decode!` calls this after *every* sync, including one whose frame failed LDPC/CRC or
    # the PRN gate and left `raw_data` holding the previous frame. Promoting that would
    # publish the previous SOW against this frame's sync epoch, one frame (3 s) low (#97).
    is_b2a_frame_decoded(state) || return state
    if has_positioning_messages(state.raw_data)
        # Promotion publishes this frame's SOW, so the symbol counter re-anchors with it;
        # the two must move together for `tow + num_bits_since_time_reference/rate` (see
        # `get_broadcast_time_of_week`) to stay on true time through frames the gate skips.
        # `raw_data.SOW` is the just-decoded frame's (checked above).
        return GNSSDecoderState(
            state;
            data = publish_data(state.cache.storage, state.raw_data),
            num_bits_since_time_reference = state.constants.sync_unit_length +
                                            state.constants.preamble_length,
        )
    end
    return state
end

# ---- Per-message-type bit-field extraction (ICD §6.2.3, Figures 6-3..6-20) -----
#
# `word` is the 288-bit message packed MSB-first into a `UInt320` (bit 1 = first PRN bit),
# read by 1-based start bit and length via `get_bits` / `get_twos_complement_num` /
# `get_bit`. Bit positions are from ICD Figures 6-3..6-20; scale factors and signedness
# from the §7 tables (Tables 7-2..7-21).
"""
Message type 10 — WN, integrity flags, IODE, ephemeris I (ICD Fig 6-3 / 6-11, Table 7-8).
"""
function parse_b2a_mt10(raw::BeiDouCNAV2Data, word::UInt320)
    word_length = B2A_MESSAGE_BITS
    SOW = raw.SOW
    # split: a Union keyword value takes the allocating kw path on Julia 1.10
    @split_nothing SOW BeiDouCNAV2Data(
        raw;
        WN = Int64(get_bits(word, word_length, 31, 13)),
        DIF_B2a = get_bit(word, word_length, 44),
        SIF_B2a = get_bit(word, word_length, 45),
        AIF_B2a = get_bit(word, word_length, 46),
        SISMAI = Int64(get_bits(word, word_length, 47, 4)),
        DIF_B1C = get_bit(word, word_length, 51),
        SIF_B1C = get_bit(word, word_length, 52),
        AIF_B1C = get_bit(word, word_length, 53),
        IODE = Int64(get_bits(word, word_length, 54, 8)),
        # Ephemeris I data block, bits 62-264 (Figure 6-11).
        t_0e = Int64(get_bits(word, word_length, 62, 11)) * 300,
        sat_type = Int64(get_bits(word, word_length, 73, 2)),
        ΔA = get_twos_complement_num(word, word_length, 75, 26) * 2.0^-9,
        A_dot = get_twos_complement_num(word, word_length, 101, 25) * 2.0^-21,
        Δn_0 = get_twos_complement_num(word, word_length, 126, 17) * 2.0^-44 * GNSS_PI,
        Δn_0_dot = get_twos_complement_num(word, word_length, 143, 23) * 2.0^-57 * GNSS_PI,
        M_0 = get_twos_complement_num(word, word_length, 166, 33) * 2.0^-32 * GNSS_PI,
        e = Int64(get_bits(word, word_length, 199, 33)) * 2.0^-34,
        ω = get_twos_complement_num(word, word_length, 232, 33) * 2.0^-32 * GNSS_PI,
        SOW_mt10 = SOW,
    )
end

"""
Message type 11 — HS, integrity flags, ephemeris II (ICD Fig 6-4 / 6-12, Table 7-8).
"""
function parse_b2a_mt11(raw::BeiDouCNAV2Data, word::UInt320)
    word_length = B2A_MESSAGE_BITS
    SOW = raw.SOW
    # split: a Union keyword value takes the allocating kw path on Julia 1.10
    @split_nothing SOW BeiDouCNAV2Data(
        raw;
        HS = Int64(get_bits(word, word_length, 31, 2)),
        DIF_B2a = get_bit(word, word_length, 33),
        SIF_B2a = get_bit(word, word_length, 34),
        AIF_B2a = get_bit(word, word_length, 35),
        SISMAI = Int64(get_bits(word, word_length, 36, 4)),
        DIF_B1C = get_bit(word, word_length, 40),
        SIF_B1C = get_bit(word, word_length, 41),
        AIF_B1C = get_bit(word, word_length, 42),
        # Ephemeris II data block, bits 43-264 (Figure 6-12).
        Ω_0 = get_twos_complement_num(word, word_length, 43, 33) * 2.0^-32 * GNSS_PI,
        i_0 = get_twos_complement_num(word, word_length, 76, 33) * 2.0^-32 * GNSS_PI,
        Ω_dot = get_twos_complement_num(word, word_length, 109, 19) * 2.0^-44 * GNSS_PI,
        i_dot = get_twos_complement_num(word, word_length, 128, 15) * 2.0^-44 * GNSS_PI,
        C_is = get_twos_complement_num(word, word_length, 143, 16) * 2.0^-30,
        C_ic = get_twos_complement_num(word, word_length, 159, 16) * 2.0^-30,
        C_rs = get_twos_complement_num(word, word_length, 175, 24) * 2.0^-8,
        C_rc = get_twos_complement_num(word, word_length, 199, 24) * 2.0^-8,
        C_us = get_twos_complement_num(word, word_length, 223, 21) * 2.0^-30,
        C_uc = get_twos_complement_num(word, word_length, 244, 21) * 2.0^-30,
        SOW_mt11 = SOW,
    )
end

# Common header block of message types 30-34 and 40: HS at bit 31, the six
# integrity flags and SISMAI at bits 33-42 (Figures 6-5..6-10).
function _parse_b2a_flags_block(raw::BeiDouCNAV2Data, word::UInt320)
    word_length = B2A_MESSAGE_BITS
    BeiDouCNAV2Data(
        raw;
        HS = Int64(get_bits(word, word_length, 31, 2)),
        DIF_B2a = get_bit(word, word_length, 33),
        SIF_B2a = get_bit(word, word_length, 34),
        AIF_B2a = get_bit(word, word_length, 35),
        SISMAI = Int64(get_bits(word, word_length, 36, 4)),
        DIF_B1C = get_bit(word, word_length, 40),
        SIF_B1C = get_bit(word, word_length, 41),
        AIF_B1C = get_bit(word, word_length, 42),
    )
end

# Clock correction data block (Figure 6-13, Table 7-5) at 1-based bit `start`. The IODC is
# not read here since its position varies: bit 112 in MT30/31/32, 218 in MT33 (behind the
# BGTO block and a reduced almanac), 134 in MT34 (behind the SISAIoc block; Figures
# 6-5..6-9). Each parser reads it itself.
_parse_b2a_clock_block(raw::BeiDouCNAV2Data, word::UInt320, start::Int) =
    BeiDouCNAV2Data(raw; beidou_clock_block(word, B2A_MESSAGE_BITS, start)...)

# SISAIoc data block (Figure 6-14) at 1-based bit `start`: t_op(11),
# SISAI_ocb(5), SISAI_oc1(3), SISAI_oc2(3). Raw values — ICD v1.0 defers
# their definition (§7.16).
_parse_b2a_sisai_oc_block(raw::BeiDouCNAV2Data, word::UInt320, start::Int) =
    BeiDouCNAV2Data(raw; beidou_sisai_oc_block(word, B2A_MESSAGE_BITS, start)...)

"""
Message type 30 — clock, IODC, group delay, BDGIM ionosphere (ICD Fig 6-5 / 6-13 / 6-15).
"""
function parse_b2a_mt30(raw::BeiDouCNAV2Data, word::UInt320)
    word_length = B2A_MESSAGE_BITS
    raw = _parse_b2a_flags_block(raw, word)
    raw = _parse_b2a_clock_block(raw, word, 43)
    BeiDouCNAV2Data(
        raw;
        IODC = Int64(get_bits(word, word_length, 112, 10)),
        T_GD_B2ap = get_twos_complement_num(word, word_length, 122, 12) * 2.0^-34,
        ISC_B2ad = get_twos_complement_num(word, word_length, 134, 12) * 2.0^-34,
        # BDGIM block, bits 146-219 (Figure 6-15, Table 7-10):
        beidou_bdgim_block(word, word_length, 146)...,
        T_GD_B1Cp = get_twos_complement_num(word, word_length, 220, 12) * 2.0^-34,
    )
end

"""
Message type 31 — clock, IODC, three reduced almanacs (ICD Fig 6-6 / 6-17).

Overwrites each broadcast PRN's slot of `raw.reduced_almanacs` in place — or,
while that is still `nothing`, of the preallocated `spare.reduced_almanacs`.
"""
function parse_b2a_mt31!(raw::BeiDouCNAV2Data, word::UInt320, spare::BeiDouCNAV2Data)
    word_length = B2A_MESSAGE_BITS
    raw = _parse_b2a_flags_block(raw, word)
    raw = _parse_b2a_clock_block(raw, word, 43)
    WN_a = Int(get_bits(word, word_length, 122, 13))
    t_0a = Int(get_bits(word, word_length, 135, 8)) * 2^12
    almanacs = raw.reduced_almanacs
    for start in (143, 181, 219)
        packet = beidou_reduced_almanac(word, word_length, start, WN_a, t_0a)
        isnothing(packet) && continue
        almanacs = writable_container(almanacs, spare.reduced_almanacs)
        set!(almanacs, packet.PRN_a, packet)
    end
    IODC = Int64(get_bits(word, word_length, 112, 10))
    # Only pass a concrete store: a Union keyword value takes the allocating kw
    # path on Julia 1.10 (and `nothing` would leave the field as is).
    almanacs === nothing && return BeiDouCNAV2Data(raw; IODC)
    BeiDouCNAV2Data(raw; IODC, reduced_almanacs = almanacs)
end

"""
Message type 32 — clock, IODC, Earth orientation parameters (ICD Fig 6-7 / 6-18).
"""
function parse_b2a_mt32(raw::BeiDouCNAV2Data, word::UInt320)
    word_length = B2A_MESSAGE_BITS
    raw = _parse_b2a_flags_block(raw, word)
    raw = _parse_b2a_clock_block(raw, word, 43)
    BeiDouCNAV2Data(
        raw;
        IODC = Int64(get_bits(word, word_length, 112, 10)),
        # EOP block, bits 122-259 (Figure 6-18, Table 7-18).
        beidou_eop_block(word, word_length, 122)...,
    )
end

"""
Message type 33 — clock, BGTO, one reduced almanac, IODC (ICD Fig 6-8 / 6-17 / 6-19).

Overwrites the almanac PRN's slot of `raw.reduced_almanacs` in place — or, while
that is still `nothing`, of the preallocated `spare.reduced_almanacs`.
"""
function parse_b2a_mt33!(raw::BeiDouCNAV2Data, word::UInt320, spare::BeiDouCNAV2Data)
    word_length = B2A_MESSAGE_BITS
    raw = _parse_b2a_flags_block(raw, word)
    raw = _parse_b2a_clock_block(raw, word, 43)
    WN_a = Int(get_bits(word, word_length, 228, 13))
    t_0a = Int(get_bits(word, word_length, 241, 8)) * 2^12
    packet = beidou_reduced_almanac(word, word_length, 180, WN_a, t_0a)
    raw = BeiDouCNAV2Data(
        raw;
        # BGTO block, bits 112-179 (Figure 6-19, Table 7-21).
        beidou_bgto_block(word, word_length, 112)...,
        IODC = Int64(get_bits(word, word_length, 218, 10)),
    )
    # Merged in a second step: splatting the BGTO block alongside a
    # `Union{Nothing,SlotDictionary}` keyword builds an abstractly typed
    # `NamedTuple`, which `juliac --trim` cannot resolve.
    isnothing(packet) && return raw
    reduced_almanacs = writable_container(raw.reduced_almanacs, spare.reduced_almanacs)
    set!(reduced_almanacs, packet.PRN_a, packet)
    BeiDouCNAV2Data(raw; reduced_almanacs)
end

"""
Message type 34 — SISAIoc, clock, IODC, BDT-UTC offset (ICD Fig 6-9 / 6-13 / 6-14 / 6-16).
"""
function parse_b2a_mt34(raw::BeiDouCNAV2Data, word::UInt320)
    word_length = B2A_MESSAGE_BITS
    raw = _parse_b2a_flags_block(raw, word)
    raw = _parse_b2a_sisai_oc_block(raw, word, 43)
    raw = _parse_b2a_clock_block(raw, word, 65)
    BeiDouCNAV2Data(
        raw;
        IODC = Int64(get_bits(word, word_length, 134, 10)),
        # BDT-UTC block, bits 144-240 (Figure 6-16, Table 7-20).
        beidou_bdt_utc_block(word, word_length, 144)...,
    )
end

"""
Message type 40 — SISAIoe, SISAIoc, one midi almanac (ICD Fig 6-10 / 6-14 / 6-20).

Overwrites the almanac PRN's slot of `raw.midi_almanacs` in place — or, while
that is still `nothing`, of the preallocated `spare.midi_almanacs`.
"""
function parse_b2a_mt40!(raw::BeiDouCNAV2Data, word::UInt320, spare::BeiDouCNAV2Data)
    word_length = B2A_MESSAGE_BITS
    raw = _parse_b2a_flags_block(raw, word)
    raw = BeiDouCNAV2Data(raw; SISAI_oe = Int64(get_bits(word, word_length, 43, 5)))
    raw = _parse_b2a_sisai_oc_block(raw, word, 48)
    # Midi almanac block, bits 70-225 (Figure 6-20, Table 7-13).
    alm = beidou_midi_almanac(word, word_length, 70)
    isnothing(alm) && return raw
    midi_almanacs = writable_container(raw.midi_almanacs, spare.midi_almanacs)
    set!(midi_almanacs, alm.PRN_a, alm)
    BeiDouCNAV2Data(raw; midi_almanacs)
end
