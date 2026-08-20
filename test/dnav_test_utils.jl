# Transmit-chain helpers for the BeiDou B1I/B3I legacy navigation message
# tests (D1 and D2 NAV): BCH(15,11,1) encoding, the two-block bit
# interleaver, and content-domain subframe synthesis. Consumed by
# `beidou_b1i.jl` and `beidou_b3i.jl`, mirroring how `cnav_test_utils.jl`
# serves the GPS L5I/L2CM tests.
#
# The synthesis path is written independently of the decoder (plain integer
# arithmetic on the ICD definitions, BDS-SIS-ICD-B1I-3.0 §5.1.3) so that a
# decoder bug cannot cancel out in the round trip.

# 320-bit packed-subframe type. Guarded: under `runtests.jl` an identical
# `Main.UInt320` already exists (gpsl1.jl defines it), and re-running the
# macro would warn.
isdefined(Main, :UInt320) || BitIntegers.@define_integers 320

"""
    dnav_test_bch_encode(info11) -> UInt16

Systematically encode 11 information bits (MSB first) into a 15-bit
BCH(15,11,1) codeword with generator g(X) = 1 + X + X⁴
(BDS-SIS-ICD-B1I-3.0 §5.1.3 Figure 5-2): the 4 parity bits are the remainder
of `info · X⁴` modulo g(X).
"""
function dnav_test_bch_encode(info11::Integer)
    0 <= info11 < 1 << 11 || throw(ArgumentError("11-bit info field expected"))
    r = Int(info11) << 4
    for i = 14:-1:4
        if (r >> i) & 1 == 1
            r ⊻= Int(0b10011) << (i - 4)
        end
    end
    UInt16(Int(info11) << 4 | r & 0xF)
end

"""
    dnav_test_interleave(block1, block2) -> UInt32

Bit-interleave two 15-bit BCH codewords into one transmitted 30-bit word:
X1¹X2¹X1²X2²…X1¹¹X2¹¹P1¹P2¹…P1⁴P2⁴ (BDS-SIS-ICD-B1I-3.0 §5.1.3).
"""
function dnav_test_interleave(block1::Integer, block2::Integer)
    word = UInt32(0)
    for j = 1:11  # information bits, alternating block 1 / block 2
        word = word << 1 | UInt32((block1 >> (15 - j)) & 1)
        word = word << 1 | UInt32((block2 >> (15 - j)) & 1)
    end
    for m = 1:4   # parity bits, alternating likewise
        word = word << 1 | UInt32((block1 >> (4 - m)) & 1)
        word = word << 1 | UInt32((block2 >> (4 - m)) & 1)
    end
    word
end

"""
    dnav_test_encode_subframe(content) -> UInt320

Encode a 224-bit content word (26 information bits of word 1 + 9 × 22, MSB
first, as produced by [`dnav_test_content`](@ref)) into the transmitted
300-bit subframe: word 1 sends its first 15 bits raw plus one BCH block;
words 2-10 interleave two BCH blocks each.
"""
function dnav_test_encode_subframe(content::UInt320)
    cbits(start, len) = Int((content >> (224 - start - len + 1)) & (UInt320(1) << len - 1))
    subframe = UInt320(0)
    # Word 1: 15 raw bits + BCH(15,11) over content bits 16-26.
    word1 = UInt32(cbits(1, 15)) << 15 | UInt32(dnav_test_bch_encode(cbits(16, 11)))
    subframe = subframe << 30 | word1
    for w = 2:10
        start = 26 + (w - 2) * 22
        block1 = dnav_test_bch_encode(cbits(start + 1, 11))
        block2 = dnav_test_bch_encode(cbits(start + 12, 11))
        subframe = subframe << 30 | dnav_test_interleave(block1, block2)
    end
    subframe
end

"""
    dnav_test_content(fields...) -> UInt320

Assemble a 224-bit content word from `(value, num_bits)` pairs given in
transmission order. The pairs must sum to 224 bits; negative values are
encoded in two's complement of their field width.
"""
function dnav_test_content(fields::Vararg{Tuple{Integer,Integer}})
    content = UInt320(0)
    total = 0
    for (value, num_bits) in fields
        num_bits <= 64 ||
            value >= 0 ||
            throw(ArgumentError("negative values only supported up to 64-bit fields"))
        mask = UInt320(1) << num_bits - 1
        # Two's complement of the field width: reinterpret the 64-bit value
        # and mask (wide fields are only ever non-negative fillers).
        raw = UInt320(reinterpret(UInt64, Int64(value))) & mask
        content = content << num_bits | raw
        total += num_bits
    end
    total == 224 || throw(ArgumentError("content fields sum to $total bits, expected 224"))
    content
end

# The 11-bit D1/D2 preamble (BDS-SIS-ICD-B1I-3.0 §5.2.4.1).
const DNAV_TEST_PREAMBLE = 0b11100010010

"""
    dnav_test_soft_symbols(subframes...; polarity = 1.0f0) -> Vector{Float32}

Serialize 300-bit subframes (UInt320, MSB first) into ±1 soft symbols
(positive ⇒ bit 0), appending the 11-bit preamble of a following subframe so
the decoder's both-ends preamble check can fire on the last one.
"""
function dnav_test_soft_symbols(subframes::UInt320...; polarity::Float32 = 1.0f0)
    symbols = Float32[]
    for subframe in subframes
        for i = 1:300
            bit = (subframe >> (300 - i)) & 1
            push!(symbols, bit == 1 ? -polarity : polarity)
        end
    end
    for i = 1:11
        bit = (DNAV_TEST_PREAMBLE >> (11 - i)) & 1
        push!(symbols, bit == 1 ? -polarity : polarity)
    end
    symbols
end

# ---- Content builders for the D1 subframes used in the tests -------------------
#
# Field order and widths follow the bit-allocation figures (BDS-SIS-ICD-B1I-3.0
# Figures 5-8, 5-9, 5-10, 5-11-x), with parity bits removed (content domain).

function dnav_test_d1_subframe1_content(;
    SOW,
    SatH1,
    AODC,
    URAI,
    WN,
    t_0c_raw,
    T_GD1_raw,
    T_GD2_raw,
    α_raw,     # NTuple{4}
    β_raw,     # NTuple{4}
    a_2_raw,
    a_0_raw,
    a_1_raw,
    AODE,
)
    dnav_test_content(
        (DNAV_TEST_PREAMBLE, 11),
        (0, 4),          # Rev
        (1, 3),          # FraID = 1
        (SOW, 20),
        (SatH1, 1),
        (AODC, 5),
        (URAI, 4),
        (WN, 13),
        (t_0c_raw, 17),
        (T_GD1_raw, 10),
        (T_GD2_raw, 10),
        (α_raw[1], 8),
        (α_raw[2], 8),
        (α_raw[3], 8),
        (α_raw[4], 8),
        (β_raw[1], 8),
        (β_raw[2], 8),
        (β_raw[3], 8),
        (β_raw[4], 8),
        (a_2_raw, 11),
        (a_0_raw, 24),
        (a_1_raw, 22),
        (AODE, 5),
    )
end

function dnav_test_d1_subframe2_content(;
    SOW,
    Δn_raw,
    C_uc_raw,
    M_0_raw,
    e_raw,
    C_us_raw,
    C_rc_raw,
    C_rs_raw,
    sqrt_A_raw,
    t_0e_msb2,
)
    dnav_test_content(
        (DNAV_TEST_PREAMBLE, 11),
        (0, 4),
        (2, 3),          # FraID = 2
        (SOW, 20),
        (Δn_raw, 16),
        (C_uc_raw, 18),
        (M_0_raw, 32),
        (e_raw, 32),
        (C_us_raw, 18),
        (C_rc_raw, 18),
        (C_rs_raw, 18),
        (sqrt_A_raw, 32),
        (t_0e_msb2, 2),
    )
end

function dnav_test_d1_subframe3_content(;
    SOW,
    t_0e_lsb15,
    i_0_raw,
    C_ic_raw,
    Ω_dot_raw,
    C_is_raw,
    i_dot_raw,
    Ω_0_raw,
    ω_raw,
)
    dnav_test_content(
        (DNAV_TEST_PREAMBLE, 11),
        (0, 4),
        (3, 3),          # FraID = 3
        (SOW, 20),
        (t_0e_lsb15, 15),
        (i_0_raw, 32),
        (C_ic_raw, 18),
        (Ω_dot_raw, 24),
        (C_is_raw, 18),
        (i_dot_raw, 14),
        (Ω_0_raw, 32),
        (ω_raw, 32),
        (0, 1),          # Rev
    )
end

# Almanac page layout of subframe 4 pages 1-24 / subframe 5 pages 1-6 and — with
# the trailing 2-bit field read as AmID — subframe 5 pages 11-23 (Figure 5-11-1/-6).
function dnav_test_d1_almanac_page_content(;
    FraID,
    SOW,
    Pnum,
    sqrt_A_raw,
    a_1_raw,
    a_0_raw,
    Ω_0_raw,
    e_raw,
    δi_raw,
    t_oa_raw,
    Ω_dot_raw,
    ω_raw,
    M_0_raw,
    am_field,
)
    dnav_test_content(
        (DNAV_TEST_PREAMBLE, 11),
        (0, 4),
        (FraID, 3),
        (SOW, 20),
        (0, 1),          # Rev
        (Pnum, 7),
        (sqrt_A_raw, 24),
        (a_1_raw, 11),
        (a_0_raw, 11),
        (Ω_0_raw, 24),
        (e_raw, 17),
        (δi_raw, 16),
        (t_oa_raw, 8),
        (Ω_dot_raw, 17),
        (ω_raw, 24),
        (M_0_raw, 24),
        (am_field, 2),
    )
end

# Subframe 5 page 9 (Figure 5-11-4): BDT-GNSS time offsets.
function dnav_test_d1_subframe5_page9_content(;
    SOW,
    A_0GPS_raw,
    A_1GPS_raw,
    A_0Gal_raw,
    A_1Gal_raw,
    A_0GLO_raw,
    A_1GLO_raw,
)
    dnav_test_content(
        (DNAV_TEST_PREAMBLE, 11),
        (0, 4),
        (5, 3),          # FraID = 5
        (SOW, 20),
        (0, 1),          # Rev
        (9, 7),          # Pnum = 9
        (0, 30),         # Rev (2 + 22 + 6)
        (A_0GPS_raw, 14),
        (A_1GPS_raw, 16),
        (A_0Gal_raw, 14),
        (A_1Gal_raw, 16),
        (A_0GLO_raw, 14),
        (A_1GLO_raw, 16),
        (0, 58),         # Rev
    )
end

# Subframe 5 page 10 (Figure 5-11-5): BDT-UTC parameters.
function dnav_test_d1_subframe5_page10_content(;
    SOW,
    Δt_LS_raw,
    Δt_LSF_raw,
    WN_LSF,
    A_0UTC_raw,
    A_1UTC_raw,
    DN,
)
    dnav_test_content(
        (DNAV_TEST_PREAMBLE, 11),
        (0, 4),
        (5, 3),          # FraID = 5
        (SOW, 20),
        (0, 1),          # Rev
        (10, 7),         # Pnum = 10
        (Δt_LS_raw, 8),
        (Δt_LSF_raw, 8),
        (WN_LSF, 8),
        (A_0UTC_raw, 32),
        (A_1UTC_raw, 24),
        (DN, 8),
        (0, 90),         # Rev
    )
end

# Subframe 5 page 8 (Figure 5-11-3): health 20-30 + almanac reference time.
function dnav_test_d1_subframe5_page8_content(; SOW, health_codes, WN_a, t_oa_raw)
    length(health_codes) == 11 || throw(ArgumentError("11 health codes expected"))
    dnav_test_content(
        (DNAV_TEST_PREAMBLE, 11),
        (0, 4),
        (5, 3),
        (SOW, 20),
        (0, 1),
        (8, 7),          # Pnum = 8
        ntuple(i -> (Int(health_codes[i]), 9), 11)...,
        (WN_a, 8),
        (t_oa_raw, 8),
        (0, 63),         # Rev
    )
end

# ---- Content builders for the D2 subframe-1 pages -------------------------------
#
# Layouts per BDS-SIS-ICD-B1I-3.0 Figures 5-14-1..-10 (150-bit basic-nav halves;
# the trailing five words are integrity/differential data outside this decoder's
# scope and are filled with zeros here — words 6-10 of the content domain).

function dnav_test_d2_page_content(pnum1::Int, tail::Vector{Tuple{Int,Int}}; SOW)
    # Common header: Pre, Rev, FraID = 1, SOW, Pnum1; then the page-specific
    # fields of words 2-5 (68 content bits), then zeros for words 6-10.
    fields = Tuple{Integer,Integer}[
        (DNAV_TEST_PREAMBLE, 11),
        (0, 4),
        (1, 3),
        (SOW, 20),
        (pnum1, 4),
    ]
    append!(fields, tail)
    used = sum(f[2] for f in fields)
    push!(fields, (0, 224 - used))
    dnav_test_content(fields...)
end
