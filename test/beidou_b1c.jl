using Test
using Random
import Aff3ct
using GNSSDecoder
using Dictionaries
using GNSSDecoder: crc24q
using GNSSSignals

# ---------------------------------------------------------------------------
# Synthetic B-CNAV1 frame generator (test-only, mirrors the transmit chain of
# BDS-SIS-ICD-B1C-1.0 §6.2).
#
# A frame is 1800 channel symbols: 72 BCH-encoded subframe-1 symbols (21 PRN
# + 51 SOH) followed by the 1728-symbol block-interleaved concatenation of
# the 1200-symbol LDPC-encoded subframe 2 and the 528-symbol LDPC-encoded
# subframe 3. Symbols are emitted as ±1 Float32 (bit 0 ⇒ +1, bit 1 ⇒ -1).
#
# The LDPC encoding runs through the *non-binary* GF(2⁶) reference encoder
# from scripts/generate_beidou_alist.jl (`beidou_ldpc_encode`), i.e. a path
# independent of the binary-image `.alist` matrices the decoder uses. The
# interleaver is an independent reimplementation of the ICD §6.2.2.4
# staggered write (36×48, columns read out), against which the decoder's
# row-order constants are cross-checked.
# ---------------------------------------------------------------------------

const _B1C_REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
# GF(2⁶) arithmetic + reference encoder (gf64_*, beidou_ldpc_encode,
# BEIDOU_LDPC_CODES). Guarded: test/beidou_ldpc.jl may have included it already.
isdefined(@__MODULE__, :beidou_ldpc_encode) ||
    isdefined(Main, :beidou_ldpc_encode) ||
    include(joinpath(_B1C_REPO_ROOT, "scripts", "generate_beidou_alist.jl"))

"""
Hard bits (first-transmitted first) of a packed codeword.
"""
_cw_bits(cw::UInt64, len::Int) = [Int((cw >> i) & 1) for i = 0:(len-1)]

"""
72 subframe-1 bits for (prn, soh): BCH(21,6) PRN codeword then BCH(51,8) SOH codeword.
"""
_b1c_sf1_bits(prn::Int, soh::Int) = vcat(
    _cw_bits(GNSSDecoder.b1c_prn_codeword(prn), 21),
    _cw_bits(GNSSDecoder.b1c_soh_codeword(soh), 51),
)

"""
Sequentially pack `(len, value)` fields MSB-first into a `total`-bit Bool
vector, then append a CRC-24Q over everything before the 24 CRC bits.
Asserts the fields fill the block exactly.
"""
function _pack_block_with_crc(total::Int, fields)
    bits = falses(total)
    cursor = 1
    for (len, val) in fields
        v = UInt64(unsigned(val) & ((UInt64(1) << len) - 1))
        for i = 0:(len-1)
            bits[cursor+i] = ((v >> (len - 1 - i)) & 1) == 1
        end
        cursor += len
    end
    cursor == total - 24 + 1 || error("fields fill $(cursor-1) bits, expected $(total-24)")
    crc = crc24q(collect(bits[1:(total-24)]))
    for i = 0:23
        bits[total-24+1+i] = ((crc >> (23 - i)) & UInt32(1)) == 1
    end
    return collect(bits)
end

"""
GF(2⁶)-encode an info-bit block through the non-binary reference encoder.
"""
function _b1c_ldpc_encode_bits(code, info_bits::AbstractVector{Bool})
    codeword = beidou_ldpc_encode(code, gf64_bits_to_symbols(info_bits))
    Int.(gf64_symbols_to_bits(codeword))
end

"""
Independent ICD §6.2.2.4 interleaver: staggered row write (rows 3k+1, 3k+2
subframe 2, row 3k+3 subframe 3 for k = 0..10, then rows 34-36 the remaining
subframe-2 symbols), column read-out (`vec` of the matrix, Julia being
column-major, is exactly the top-to-bottom column read).
"""
function _b1c_interleave(sf2::Vector{Int}, sf3::Vector{Int})
    @assert length(sf2) == 1200 && length(sf3) == 528
    A = Matrix{Int}(undef, 36, 48)
    chunk(v, i) = v[((i-1)*48+1):(i*48)]
    sf2_next = 1
    sf3_next = 1
    for k = 0:10
        A[3k+1, :] = chunk(sf2, sf2_next);
        sf2_next += 1
        A[3k+2, :] = chunk(sf2, sf2_next);
        sf2_next += 1
        A[3k+3, :] = chunk(sf3, sf3_next);
        sf3_next += 1
    end
    for r = 34:36
        A[r, :] = chunk(sf2, sf2_next);
        sf2_next += 1
    end
    @assert sf2_next == 26 && sf3_next == 12
    return vec(A)
end

"""
One 1800-symbol frame as ±1 Float32 soft symbols.
"""
function _b1c_frame_symbols(
    prn::Int,
    soh::Int,
    sf2_bits::AbstractVector{Bool},
    sf3_bits::AbstractVector{Bool},
)
    sf2_enc = _b1c_ldpc_encode_bits(BEIDOU_LDPC_CODES.bcnv1_sf2, sf2_bits)
    sf3_enc = _b1c_ldpc_encode_bits(BEIDOU_LDPC_CODES.bcnv1_sf3, sf3_bits)
    bits = vcat(_b1c_sf1_bits(prn, soh), _b1c_interleave(sf2_enc, sf3_enc))
    Float32[b == 0 ? 1.0f0 : -1.0f0 for b in bits]
end

# --- Golden subframe-2 field values (ICD Figure 6-6 / 6-12 .. 6-14) ---------
# Raw (unscaled) field values, packed in figure order; expected decoded
# values asserted with explicit ICD scale factors below.
const _B1C_G = (
    WN = 800,
    HOW = 101,
    IODC = 517,
    IODE = 5,
    t_oe_raw = 1200,
    sat_type = 3,
    ΔA_raw = -8000,
    A_dot_raw = 1024,
    Δn_0_raw = -5000,
    Δn_0_dot_raw = 300000,
    M_0_raw = -1234567890,
    e_raw = 12345678,
    ω_raw = 987654321,
    Ω_0_raw = -1122334455,
    i_0_raw = 666777888,
    Ω_dot_raw = -100000,
    i_dot_raw = -4000,
    C_is_raw = -2000,
    C_ic_raw = 1500,
    C_rs_raw = -300000,
    C_rc_raw = 250000,
    C_us_raw = -80000,
    C_uc_raw = 70000,
    t_oc_raw = 1200,
    a_0_raw = -1000000,
    a_1_raw = 200000,
    a_2_raw = -500,
    T_GD_B2ap_raw = -100,
    ISC_B1Cd_raw = 50,
    T_GD_B1Cp_raw = -75,
)

"""
600-bit subframe 2 for the golden values (fields in ICD Figure 6-6 order).
"""
function _golden_sf2_bits(g = _B1C_G)
    _pack_block_with_crc(
        600,
        [
            (13, g.WN),
            (8, g.HOW),
            (10, g.IODC),
            (8, g.IODE),
            # Ephemeris I (Figure 6-12)
            (11, g.t_oe_raw),
            (2, g.sat_type),
            (26, g.ΔA_raw),
            (25, g.A_dot_raw),
            (17, g.Δn_0_raw),
            (23, g.Δn_0_dot_raw),
            (33, g.M_0_raw),
            (33, g.e_raw),
            (33, g.ω_raw),
            # Ephemeris II (Figure 6-13)
            (33, g.Ω_0_raw),
            (33, g.i_0_raw),
            (19, g.Ω_dot_raw),
            (15, g.i_dot_raw),
            (16, g.C_is_raw),
            (16, g.C_ic_raw),
            (24, g.C_rs_raw),
            (24, g.C_rc_raw),
            (21, g.C_us_raw),
            (21, g.C_uc_raw),
            # Clock correction (Figure 6-14)
            (11, g.t_oc_raw),
            (25, g.a_0_raw),
            (22, g.a_1_raw),
            (11, g.a_2_raw),
            # Group delays + Rev
            (12, g.T_GD_B2ap_raw),
            (12, g.ISC_B1Cd_raw),
            (12, g.T_GD_B1Cp_raw),
            (7, 0),
        ],
    )
end

# --- Golden subframe-3 pages (ICD Figures 6-8 .. 6-11) ----------------------

# Common page header after the PageID: HS DIF SIF AIF SISMAI.
_sf3_header(page; hs = 0) = [(6, page), (2, hs), (1, 0), (1, 0), (1, 1), (4, 5)]

"""
Page type 1: SISAI_oe + SISAI_oc + BDGIM iono + BDT-UTC (Figure 6-8).
"""
function _golden_sf3_page1_bits(; hs = 0)
    _pack_block_with_crc(
        264,
        vcat(
            _sf3_header(1; hs),
            [
                (5, 12),                               # SISAI_oe
                (11, 900),
                (5, 7),
                (3, 2),
                (3, 3),     # SISAI_oc: t_op ocb oc1 oc2
                # BDGIM α₁..α₉ (Figure 6-16)
                (10, 800),
                (8, -40),
                (8, 16),
                (8, 24),
                (8, 32),
                (8, -8),
                (8, 12),
                (8, -16),
                (8, 20),
                # BDT-UTC (Figure 6-17)
                (16, -12000),
                (13, 3000),
                (7, -50),
                (8, 4),
                (16, 30000),
                (13, 800),
                (13, 801),
                (3, 3),
                (8, 5),
                (27, 0),                               # Rev
            ],
        ),
    )
end

"""
Page type 2: SISAI_oc + WN_a/t_oa + four reduced almanacs (Figure 6-9).
"""
function _golden_sf3_page2_bits()
    red(prn) = [(6, prn), (2, 3), (8, -30), (7, 20), (7, -25), (8, 48)]
    _pack_block_with_crc(
        264,
        vcat(
            _sf3_header(2),
            [(11, 900), (5, 7), (3, 2), (3, 3)],       # SISAI_oc
            [(13, 812), (8, 100)],                     # WN_a, t_oa
            red(7),
            red(8),
            [(6, 0), (2, 0), (8, 0), (7, 0), (7, 0), (8, 0)],
            red(9),
            [(30, 0)],                                 # Rev
        ),
    )
end

"""
Page type 3: SISAI_oe + EOP + BGTO (Figure 6-10).
"""
function _golden_sf3_page3_bits()
    _pack_block_with_crc(
        264,
        vcat(
            _sf3_header(3),
            [
                (5, 12),                               # SISAI_oe
                # EOP (Figure 6-19)
                (16, 20000),
                (21, -300000),
                (15, 5000),
                (21, 250000),
                (15, -6000),
                (31, -20000000),
                (19, 30000),
                # BGTO (Figure 6-20): GNSS ID 1 = GPS
                (3, 1),
                (13, 810),
                (16, 20000),
                (16, -5000),
                (13, 1200),
                (7, -30),
                (14, 0),                               # Rev
            ],
        ),
    )
end

"""
Page type 4: SISAI_oc + one midi almanac (Figure 6-11).
"""
function _golden_sf3_page4_bits()
    _pack_block_with_crc(
        264,
        vcat(
            _sf3_header(4),
            [
                (11, 900),
                (5, 7),
                (3, 2),
                (3, 3),     # SISAI_oc
                # Midi almanac (Figure 6-21)
                (6, 25),
                (2, 2),
                (13, 812),
                (8, 100),
                (11, 500),
                (11, -256),
                (17, 103894),
                (16, -20000),
                (11, -60),
                (16, 15000),
                (16, -15000),
                (11, -800),
                (10, 100),
                (8, 64),
                (47, 0),                               # Rev
            ],
        ),
    )
end

@testset "BeiDou B1C (B-CNAV1)" begin
    prn = 30

    @testset "Subframe-1 BCH codeword tables (ICD §6.2.2.1)" begin
        # Systematic prefix: the first k emitted symbols spell the info value
        # MSB-first (consequence of the bit-reversed seeding).
        for p = 1:63
            bits = _cw_bits(GNSSDecoder.b1c_prn_codeword(p), 21)
            @test foldl((a, b) -> 2a + b, bits[1:6]) == p
        end
        for soh = 0:199
            bits = _cw_bits(GNSSDecoder.b1c_soh_codeword(soh), 51)
            @test foldl((a, b) -> 2a + b, bits[1:8]) == soh
        end
        # LFSR recurrence with the ICD characteristic polynomials: for
        # g21,6 = x⁶+x⁴+x²+x+1, c[i+6] = c[i+4] ⊻ c[i+2] ⊻ c[i+1] ⊻ c[i];
        # for g51,8 = x⁸+x⁷+x⁴+x³+x²+x+1,
        # c[i+8] = c[i+7] ⊻ c[i+4] ⊻ c[i+3] ⊻ c[i+2] ⊻ c[i+1] ⊻ c[i].
        for p = 1:63
            c = _cw_bits(GNSSDecoder.b1c_prn_codeword(p), 21)
            @test all(c[i+7] == c[i+5] ⊻ c[i+3] ⊻ c[i+2] ⊻ c[i+1] for i = 0:14)
        end
        for soh = 0:199
            c = _cw_bits(GNSSDecoder.b1c_soh_codeword(soh), 51)
            @test all(
                c[i+9] == c[i+8] ⊻ c[i+5] ⊻ c[i+4] ⊻ c[i+3] ⊻ c[i+2] ⊻ c[i+1] for i = 0:42
            )
        end
        # All PRN codewords distinct; complements are never codewords, so the
        # polarity of a subframe-1 match is unambiguous.
        prn_cws = [GNSSDecoder.b1c_prn_codeword(p) for p = 1:63]
        @test length(unique(prn_cws)) == 63
        mask21 = (UInt64(1) << 21) - 1
        @test isempty(intersect(Set(prn_cws), Set(cw ⊻ mask21 for cw in prn_cws)))
        soh_cws = collect(GNSSDecoder.B1C_SOH_CODEWORDS)
        @test length(unique(soh_cws)) == 200
    end

    @testset "Interleaver row orders match the ICD staggered write" begin
        # Interleave marker payloads through the independent transmitter-side
        # implementation, then invert with the decoder's constants.
        sf2 = collect(1:1200)
        sf3 = collect(2001:2528)
        interleaved = _b1c_interleave(sf2, sf3)
        rows = GNSSDecoder.deinterleave(interleaved, 36, 48)
        sf2_back =
            vcat((rows[((r-1)*48+1):(r*48)] for r in GNSSDecoder.B1C_SF2_ROW_ORDER)...)
        sf3_back =
            vcat((rows[((r-1)*48+1):(r*48)] for r in GNSSDecoder.B1C_SF3_ROW_ORDER)...)
        @test sf2_back == sf2
        @test sf3_back == sf3
    end

    @testset "Full frame decode (golden subframe 2 + paged subframe 3)" begin
        sf2 = Bool.(_golden_sf2_bits())
        frames = [
            _b1c_frame_symbols(prn, 57, sf2, Bool.(_golden_sf3_page1_bits())),
            _b1c_frame_symbols(prn, 58, sf2, Bool.(_golden_sf3_page2_bits())),
            _b1c_frame_symbols(prn, 59, sf2, Bool.(_golden_sf3_page3_bits())),
            _b1c_frame_symbols(prn, 60, sf2, Bool.(_golden_sf3_page4_bits())),
        ]
        stream = vcat(
            frames...,
            _b1c_frame_symbols(prn, 61, sf2, Bool.(_golden_sf3_page1_bits()))[1:72],
        )

        state = BeiDouB1CDecoderState(prn)
        @test !is_decoding_completed_for_positioning(state)
        state = decode(state, stream, length(stream))

        d = state.data
        @test is_decoding_completed_for_positioning(state)
        # A reserved orbit type (0) leaves the semi-major axis unknowable: it
        # is what selects the `A_ref` the broadcast `ΔA` corrects, so the
        # satellite must not be used (ICD Table 7-6).
        @test !GNSSDecoder.is_subframe2_decoded(BeiDouB1CData(d; sat_type = 0))
        @test all(
            GNSSDecoder.is_subframe2_decoded(BeiDouB1CData(d; sat_type = t)) for t = 1:3
        )
        @test !state.is_shifted_by_180_degrees
        # The SOH stamps the leading edge of *this* subframe 1 (ICD §7.3), so
        # the armed counter spans the whole frame plus the next frame's
        # subframe 1: 1800 + 72 = 1872 symbols = 18.72 s at 100 sps. (GPS
        # L1C-D's TOI stamps the *next* frame and anchors to 52 alone; using
        # that convention here would report the time a full 18 s early.)
        @test state.num_bits_after_valid_syncro_sequence == 1872
        # The last decoded frame carried SOH = 60.
        @test d.soh == 60
        @test d.WN == 800
        @test d.HOW == 101
        @test d.IODC == 517
        @test d.IODE == 5
        # Ephemeris (ICD Table 7-8 scale factors).
        @test d.t_oe == 1200 * 300
        @test d.sat_type == 3
        @test d.ΔA ≈ -8000 * 2.0^-9
        @test d.A_dot ≈ 1024 * 2.0^-21
        @test d.Δn_0 ≈ -5000 * 2.0^-44 * 3.1415926535898
        @test d.Δn_0_dot ≈ 300000 * 2.0^-57 * 3.1415926535898
        @test d.M_0 ≈ -1234567890 * 2.0^-32 * 3.1415926535898
        @test d.e ≈ 12345678 * 2.0^-34
        @test d.ω ≈ 987654321 * 2.0^-32 * 3.1415926535898
        @test d.Ω_0 ≈ -1122334455 * 2.0^-32 * 3.1415926535898
        @test d.i_0 ≈ 666777888 * 2.0^-32 * 3.1415926535898
        @test d.Ω_dot ≈ -100000 * 2.0^-44 * 3.1415926535898
        @test d.i_dot ≈ -4000 * 2.0^-44 * 3.1415926535898
        @test d.C_is ≈ -2000 * 2.0^-30
        @test d.C_ic ≈ 1500 * 2.0^-30
        @test d.C_rs ≈ -300000 * 2.0^-8
        @test d.C_rc ≈ 250000 * 2.0^-8
        @test d.C_us ≈ -80000 * 2.0^-30
        @test d.C_uc ≈ 70000 * 2.0^-30
        # Clock + group delay (Tables 7-5, 7-6).
        @test d.t_oc == 1200 * 300
        @test d.a_0 ≈ -1000000 * 2.0^-34
        @test d.a_1 ≈ 200000 * 2.0^-50
        @test d.a_2 ≈ -500 * 2.0^-66
        @test d.T_GD_B2ap ≈ -100 * 2.0^-34
        @test d.ISC_B1Cd ≈ 50 * 2.0^-34
        @test d.T_GD_B1Cp ≈ -75 * 2.0^-34
        # Page header fields (all pages carried HS=0, AIF=1, SISMAI=5).
        @test d.hs == 0
        @test d.dif === false
        @test d.sif === false
        @test d.aif === true
        @test d.sismai == 5
        @test d.sisai_oe == 12
        @test d.t_op == 900
        @test d.sisai_ocb == 7
        @test d.sisai_oc1 == 2
        @test d.sisai_oc2 == 3
        @test is_sat_healthy(state)
        # Page 1: BDGIM iono (Table 7-10; α₅ scale −2⁻³) + BDT-UTC (Table 7-20).
        @test d.α_1 ≈ 800 * 2.0^-3
        @test d.α_2 ≈ -40 * 2.0^-3
        @test d.α_3 ≈ 16 * 2.0^-3
        @test d.α_4 ≈ 24 * 2.0^-3
        @test d.α_5 ≈ -(32 * 2.0^-3)
        @test d.α_6 ≈ -8 * 2.0^-3
        @test d.α_7 ≈ 12 * 2.0^-3
        @test d.α_8 ≈ -16 * 2.0^-3
        @test d.α_9 ≈ 20 * 2.0^-3
        @test d.A_0UTC ≈ -12000 * 2.0^-35
        @test d.A_1UTC ≈ 3000 * 2.0^-51
        @test d.A_2UTC ≈ -50 * 2.0^-68
        @test d.Δt_LS == 4
        @test d.t_ot == 30000 * 2^4
        @test d.WN_ot == 800
        @test d.WN_LSF == 801
        @test d.DN == 3
        @test d.Δt_LSF == 5
        # Page 2: reduced almanacs for PRN 7, 8, 9 (the zero-PRN slot skipped).
        @test length(d.reduced_almanacs) == 3
        red = d.reduced_almanacs[7]
        @test red.PRN_a == 7
        @test red.sat_type == 3
        @test red.WN_a == 812
        @test red.t_oa == 100 * 2^12
        @test red.δA ≈ -30 * 2.0^9
        @test red.Ω_0 ≈ 20 * 2.0^-6 * 3.1415926535898
        @test red.Φ_0 ≈ -25 * 2.0^-6 * 3.1415926535898
        @test red.health == 48
        @test haskey(d.reduced_almanacs, 8) && haskey(d.reduced_almanacs, 9)
        # Page 3: EOP (Table 7-18) + BGTO for GPS (Table 7-21).
        @test d.t_EOP == 20000 * 2^4
        @test d.PM_X ≈ -300000 * 2.0^-20
        @test d.PM_X_dot ≈ 5000 * 2.0^-21
        @test d.PM_Y ≈ 250000 * 2.0^-20
        @test d.PM_Y_dot ≈ -6000 * 2.0^-21
        @test d.ΔUT1 ≈ -20000000 * 2.0^-24
        @test d.ΔUT1_dot ≈ 30000 * 2.0^-25
        @test length(d.bgtos) == 1
        bgto = d.bgtos[1]
        @test bgto.WN_0BGTO == 810
        @test bgto.t_0BGTO == 20000 * 2^4
        @test bgto.A_0BGTO ≈ -5000 * 2.0^-35
        @test bgto.A_1BGTO ≈ 1200 * 2.0^-51
        @test bgto.A_2BGTO ≈ -30 * 2.0^-68
        # Page 4: midi almanac for PRN 25 (Table 7-13).
        @test length(d.midi_almanacs) == 1
        alm = d.midi_almanacs[25]
        @test alm.PRN_a == 25
        @test alm.sat_type == 2
        @test alm.WN_a == 812
        @test alm.t_oa == 100 * 2^12
        @test alm.e ≈ 500 * 2.0^-16
        @test alm.δi ≈ -256 * 2.0^-14 * 3.1415926535898
        @test alm.sqrt_A ≈ 103894 * 2.0^-4
        @test alm.Ω_0 ≈ -20000 * 2.0^-15 * 3.1415926535898
        @test alm.Ω_dot ≈ -60 * 2.0^-33 * 3.1415926535898
        @test alm.ω ≈ 15000 * 2.0^-15 * 3.1415926535898
        @test alm.M_0 ≈ -15000 * 2.0^-15 * 3.1415926535898
        @test alm.a_f0 ≈ -800 * 2.0^-20
        @test alm.a_f1 ≈ 100 * 2.0^-37
        @test alm.health == 64
        @test d.num_sf3_pages_received == 4
        # Signal metadata forwarding.
        @test get_signal_type(state) == BeiDouB1C_D
        @test get_data_frequency(state) == get_data_frequency(BeiDouB1C_D)
        @test get_time_system_name(state) == "BeiDou Time"
        @test get_constellation_name(state) == "BeiDou"
    end

    @testset "180° polarity flip decodes identically" begin
        sf2 = Bool.(_golden_sf2_bits())
        stream = vcat(
            _b1c_frame_symbols(prn, 10, sf2, Bool.(_golden_sf3_page1_bits())),
            _b1c_frame_symbols(prn, 11, sf2, Bool.(_golden_sf3_page1_bits()))[1:72],
        )
        state = decode(BeiDouB1CDecoderState(prn), -stream, length(stream))
        @test state.is_shifted_by_180_degrees
        @test state.data.soh == 10
        @test state.data.WN == 800
        @test is_decoding_completed_for_positioning(state)
    end

    @testset "Unhealthy satellite (HS = 1) is reported" begin
        sf2 = Bool.(_golden_sf2_bits())
        stream = vcat(
            _b1c_frame_symbols(prn, 10, sf2, Bool.(_golden_sf3_page1_bits(hs = 1))),
            _b1c_frame_symbols(prn, 11, sf2, Bool.(_golden_sf3_page1_bits(hs = 1)))[1:72],
        )
        state = decode(BeiDouB1CDecoderState(prn), stream, length(stream))
        @test state.data.hs == 1
        @test is_decoding_completed_for_positioning(state)
        @test !is_sat_healthy(state)
    end

    @testset "Corrupted subframe-2 CRC is rejected, subframe 3 still decodes" begin
        # Encode a subframe 2 whose trailing CRC deliberately mismatches: the
        # LDPC decode converges cleanly to the (wrong-CRC) info block, and the
        # CRC gate must drop it.
        bad_sf2 = Bool.(_golden_sf2_bits())
        bad_sf2[600] = !bad_sf2[600]  # break the CRC itself
        stream = vcat(
            _b1c_frame_symbols(prn, 10, bad_sf2, Bool.(_golden_sf3_page1_bits())),
            _b1c_frame_symbols(prn, 11, bad_sf2, Bool.(_golden_sf3_page1_bits()))[1:72],
        )
        state = decode(BeiDouB1CDecoderState(prn), stream, length(stream))
        @test isnothing(state.raw_data.WN)          # SF2 dropped
        @test state.raw_data.soh == 10              # sync still locked
        @test state.raw_data.num_sf3_pages_received == 1   # SF3 decoded
        @test state.raw_data.hs == 0
        @test !is_decoding_completed_for_positioning(state)  # promotion gated
        @test !is_sat_healthy(state)
    end

    @testset "SOH discontinuity resets the decoder, then re-locks" begin
        sf2 = Bool.(_golden_sf2_bits())
        sf3 = Bool.(_golden_sf3_page1_bits())
        # Frames with SOH 10, 11, then a jump to 50, 51, 52 (+ next SF1).
        # Frame 10 decodes (10/11 pair); the buffer then never finds an
        # 11/12 pair, slides to the 50/51 pair, whose decode breaks the
        # monotonic-SOH check against the locked SOH=10 and resets the
        # decoder (emptying the buffer mid-frame-51); the stream finally
        # re-locks and decodes on the 52/53-tail pair.
        stream = vcat(
            _b1c_frame_symbols(prn, 10, sf2, sf3),
            _b1c_frame_symbols(prn, 11, sf2, sf3),
            _b1c_frame_symbols(prn, 50, sf2, sf3),
            _b1c_frame_symbols(prn, 51, sf2, sf3),
            _b1c_frame_symbols(prn, 52, sf2, sf3),
            _b1c_frame_symbols(prn, 53, sf2, sf3)[1:72],
        )
        state = decode(BeiDouB1CDecoderState(prn), stream, length(stream))
        @test state.raw_data.soh == 52
        @test is_decoding_completed_for_positioning(state)
    end

    @testset "Frames of another PRN never sync" begin
        sf2 = Bool.(_golden_sf2_bits())
        sf3 = Bool.(_golden_sf3_page1_bits())
        stream = vcat(
            _b1c_frame_symbols(prn, 10, sf2, sf3),
            _b1c_frame_symbols(prn, 11, sf2, sf3)[1:72],
        )
        state = decode(BeiDouB1CDecoderState(prn + 1), stream, length(stream))
        @test isnothing(state.raw_data.soh)
        @test state.raw_data.num_sf3_pages_received == 0
        @test !is_decoding_completed_for_positioning(state)
    end

    @testset "Noisy symbols still decode (LDPC margin)" begin
        rng = MersenneTwister(0xB1C)
        sf2 = Bool.(_golden_sf2_bits())
        stream = vcat(
            _b1c_frame_symbols(prn, 10, sf2, Bool.(_golden_sf3_page1_bits())),
            _b1c_frame_symbols(prn, 11, sf2, Bool.(_golden_sf3_page1_bits()))[1:72],
        )
        # Add noise only to the LDPC-protected payload symbols — the
        # subframe-1 sync match is a hard (error-free) codeword comparison.
        # The flooding sum-product LDPC decode is scale-sensitive (see the
        # soft-symbol convention note on `decode`), so the ±1 ± σ symbols are
        # fed on the matching LLR scale 2r/σ².
        σ = 0.4f0
        noisy = copy(stream)
        for i in eachindex(noisy)
            in_sf1 = mod1(i, 1800) <= 72 && i <= 1872
            in_sf1 || (noisy[i] += σ * randn(rng, Float32))
        end
        noisy .*= 2 / σ^2
        state = decode(BeiDouB1CDecoderState(prn), noisy, length(noisy))
        @test state.data.WN == 800
        @test is_decoding_completed_for_positioning(state)
    end

    @testset "Dispatch from the GNSSSignals type + reset" begin
        state = GNSSDecoderState(BeiDouB1C_D(), prn)
        @test state.prn == prn
        @test get_signal_type(state) == BeiDouB1C_D
        @test_throws ArgumentError BeiDouB1CDecoderState(0)
        @test_throws ArgumentError BeiDouB1CDecoderState(64)

        sf2 = Bool.(_golden_sf2_bits())
        stream = vcat(
            _b1c_frame_symbols(prn, 10, sf2, Bool.(_golden_sf3_page1_bits())),
            _b1c_frame_symbols(prn, 11, sf2, Bool.(_golden_sf3_page1_bits()))[1:72],
        )
        state = decode(BeiDouB1CDecoderState(prn), stream, length(stream))
        @test is_decoding_completed_for_positioning(state)
        state = reset_decoder_state(state)
        @test isnothing(state.raw_data.soh)
        @test state.raw_data.WN == 800   # ephemeris survives for warm re-lock
        @test isnothing(state.data.WN)   # validated data cleared
        @test GNSSDecoder.num_bits_buffered(state) == 0
    end
end
