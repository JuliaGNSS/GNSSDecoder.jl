# Shared BeiDou legacy navigation message core (D1 NAV and D2 NAV), broadcast
# identically on B1I (BDS-SIS-ICD-B1I-3.0 §5) and B3I (BDS-SIS-ICD-B3I-1.0 §5).
# The signal layers in `b1i.jl` / `b3i.jl` are thin wrappers over this file,
# mirroring how `gps/cnav.jl` is shared by GPS L5I and L2CM.
#
# Message formats (BDS-SIS-ICD-B1I-3.0 §5.1.1):
#   - MEO/IGSO satellites (PRN 6-58) broadcast **D1** at 50 bps: five 300-bit
#     subframes per 30 s frame; subframes 1-3 carry the fundamental navigation
#     data of the broadcasting satellite, subframes 4-5 carry 24 pages each of
#     almanac, health, and time-offset data.
#   - GEO satellites (PRN 1-5, 59-63) broadcast **D2** at 500 bps: the same
#     300-bit subframe/word coding, but the fundamental navigation data is
#     spread over pages 1-10 of subframe 1 (one page per 3 s frame,
#     §5.3.2 Figures 5-14-1..-10). D2 subframes 2-4 (BDS integrity and
#     wide-area differential corrections) and subframe 5 (D2 almanac pages)
#     carry the legacy regional augmentation service and are intentionally
#     not decoded here — they are not needed for positioning readiness.
#
# Word coding (§5.1.3): every 30-bit word is BCH(15,11,1)-protected with
# g(X) = 1 + X + X⁴. Word 1 of a subframe sends its first 15 bits raw
# (preamble + FraID region) followed by one 15-bit BCH block; words 2-10
# interleave two BCH blocks bit-by-bit (X1¹X2¹X1²X2²…P1¹P2¹…P1⁴P2⁴), giving
# the pre-interleave layout [22 information bits][8 parity bits] that the
# ICD's bit-allocation figures index. The decoder de-interleaves, corrects
# single bit errors per block, and re-packs only the information bits into a
# 224-bit "content" word (26 bits of word 1 + 9 × 22 bits), in which every
# multi-word ICD field becomes contiguous — the parity-split field positions
# of the figures map 1:1 onto content positions (see `dnav_content_position`).

# ---- Frame geometry (BDS-SIS-ICD-B1I-3.0 §5.2.1 / §5.3.1) ------------------
#
# A subframe is 300 bits = 10 words of 30 bits, and the sync window is one
# subframe plus the next subframe's 11-bit preamble.

const DNAV_SUBFRAME_BITS = 300
const DNAV_WORD_BITS = 30
const DNAV_PREAMBLE_BITS = 11
const DNAV_WINDOW_BITS = DNAV_SUBFRAME_BITS + DNAV_PREAMBLE_BITS  # 311

"""
    BeiDouDNAVConstants{S}

BDCS constants and D1/D2 NAV message structure parameters for the BeiDou
B1I/B3I legacy signals, parameterized on the signal tag `S`
(`:BeiDouB1I` or `:BeiDouB3I`) like `GPSCNAVConstants`. The distinct tag
selects the signal identity ([`get_signal_type`](@ref)) — the message
structure and all constants are identical on both signals
(BDS-SIS-ICD-B3I-1.0 §5 mirrors BDS-SIS-ICD-B1I-3.0 §5).

# Fields

  - `syncro_sequence_length::Int`: Length of one subframe in bits (300)
  - `preamble::UInt16`: 11-bit preamble `11100010010` (modified Barker code, §5.2.4.1)
  - `preamble_length::Int`: Length of preamble in bits (11)
  - `word_length::Int`: Length of each word in bits (30)
  - `PI::Float64`: π = 3.1415926535898 (BDS-SIS-ICD-B1I-3.0 Table 5-11)
  - `Ω_dot_e::Float64`: BDCS Earth rotation rate = 7.2921150×10⁻⁵ rad/s (differs from WGS-84!)
  - `c::Float64`: Speed of light = 2.99792458×10⁸ m/s
  - `μ::Float64`: BDCS geocentric gravitational constant = 3.986004418×10¹⁴ m³/s²
  - `F::Float64`: Relativistic correction constant −2√μ/c² = -4.442807309×10⁻¹⁰ s/√m (§5.2.4.9)

# Reference

BDS-SIS-ICD-B1I-3.0, §5.1-§5.3 and Tables 5-7, 5-10, 5-11
"""
Base.@kwdef struct BeiDouDNAVConstants{S} <: AbstractGNSSConstants
    syncro_sequence_length::Int = DNAV_SUBFRAME_BITS
    preamble::UInt16 = 0b11100010010
    preamble_length::Int = DNAV_PREAMBLE_BITS
    word_length::Int = DNAV_WORD_BITS
    PI::Float64 = GNSS_PI
    Ω_dot_e::Float64 = BEIDOU_EARTH_ROTATION_RATE
    c::Float64 = SPEED_OF_LIGHT
    μ::Float64 = BEIDOU_μ
    F::Float64 = -4.442807309e-10
end

# ---- BCH(15,11,1) (BDS-SIS-ICD-B1I-3.0 §5.1.3) -------------------------------

# Generator polynomial g(X) = 1 + X + X⁴, MSB-first: 0b10011. Stored as an
# `Int` — a bare `0b10011` literal is a `UInt8`, whose left-shifts in the
# syndrome long division would overflow 8 bits.
const DNAV_BCH_GENERATOR = Int(0b10011)

# Codewords are 15 bits, MSB (= earliest transmitted bit, coefficient of X¹⁴)
# first. The syndrome is the remainder of the codeword polynomial modulo g(X).
function dnav_bch_syndrome(codeword::Integer)
    r = Int(codeword)
    for i = 14:-1:4
        if (r >> i) & 1 == 1
            r ⊻= DNAV_BCH_GENERATOR << (i - 4)
        end
    end
    r & 0xF
end

# Syndrome → 1-based error bit position (1 = MSB of the 15-bit codeword),
# reproducing the ICD's Figure 5-4 ROM list (e.g. syndrome 0001 flips bit 15,
# 1001 flips bit 1). Index [syndrome]; entry 0 is unused (zero syndrome =
# error-free).
const DNAV_BCH_ERROR_POSITION = let
    table = zeros(Int, 15)
    for pos = 1:15
        syndrome = dnav_bch_syndrome(1 << (15 - pos))
        table[syndrome] = pos
    end
    table
end

"""
    dnav_bch_decode(codeword) -> corrected codeword

Decode one BCH(15,11,1) block (15 bits, MSB first): correct the single bit
error indicated by a non-zero syndrome (BDS-SIS-ICD-B1I-3.0 §5.1.3, Figure
5-4). Two or more bit errors exceed the code's capability and are silently
mis-corrected — the cross-subframe data voting in `validate_data` is the
guard against such words reaching the validated data.
"""
function dnav_bch_decode(codeword::Integer)
    syndrome = dnav_bch_syndrome(codeword)
    syndrome == 0 && return Int(codeword)
    Int(codeword) ⊻ (1 << (15 - DNAV_BCH_ERROR_POSITION[syndrome]))
end

# De-interleave one transmitted 30-bit word (words 2-10) into its two 15-bit
# BCH blocks. Transmitted order is X1¹X2¹X1²X2²…X1¹¹X2¹¹P1¹P2¹P1²P2²P1³P2³P1⁴P2⁴
# (§5.1.3): the 22 information bits alternate block-1/block-2, then the 8
# parity bits alternate likewise.
function dnav_deinterleave_word(word::Integer)
    block1 = 0
    block2 = 0
    for i = 1:11  # information bits: transmitted bits 2i-1 (block 1), 2i (block 2)
        block1 = block1 << 1 | (Int(word) >> (30 - (2i - 1))) & 1
        block2 = block2 << 1 | (Int(word) >> (30 - 2i)) & 1
    end
    for m = 1:4   # parity bits: transmitted bits 21+2m (block 1), 22+2m (block 2)
        block1 = block1 << 1 | (Int(word) >> (30 - (21 + 2m))) & 1
        block2 = block2 << 1 | (Int(word) >> (30 - (22 + 2m))) & 1
    end
    block1, block2
end

# Extract transmitted word `word_number` (1-10, 30 bits) from the packed sync
# buffer (300 subframe bits + 11 trailing preamble bits, oldest at MSB).
function dnav_get_word(buffer, constants::BeiDouDNAVConstants, word_number::Int)
    num_words = constants.syncro_sequence_length ÷ constants.word_length
    shifted =
        buffer >>
        (constants.word_length * (num_words - word_number) + constants.preamble_length)
    UInt32(shifted & (UInt320(1) << constants.word_length - 1))
end

"""
    decode_dnav_subframe_content(buffer, constants) -> UInt320

BCH-decode all ten words of the 300-bit subframe held in the packed sync
buffer and concatenate the corrected *information* bits into the 224-bit
content word (word 1 bits 1-26, then bits 1-22 of each of words 2-10),
MSB-first and right-aligned. In the content domain every parity-split field
of the ICD bit-allocation figures is contiguous; positions are mapped from
the figures' 300-bit numbering by [`dnav_content_position`](@ref).
"""
function decode_dnav_subframe_content(buffer, constants::BeiDouDNAVConstants)
    content = UInt320(0)
    word1 = dnav_get_word(buffer, constants, 1)
    # Word 1: bits 1-15 raw (preamble, Rev, FraID, SOW MSB region starts at 19),
    # bits 16-30 one non-interleaved BCH(15,11,1) block (§5.2.2).
    corrected1 = dnav_bch_decode(word1 & 0x7FFF)
    # Content bits 1-26: 15 raw bits + the block's 11 information bits.
    content = content << 26 | UInt320(word1 >> 15) << 11 | UInt320(corrected1 >> 4)
    for w = 2:10
        word = dnav_get_word(buffer, constants, w)
        block1, block2 = dnav_deinterleave_word(word)
        info1 = dnav_bch_decode(block1) >> 4
        info2 = dnav_bch_decode(block2) >> 4
        content = content << 22 | UInt320(info1) << 11 | UInt320(info2)
    end
    content
end

# Number of information ("content") bits per subframe: 26 + 9 × 22.
const DNAV_CONTENT_BITS = 224

"""
    dnav_content_position(icd_position) -> Int

Map a bit position from the ICD's 300-bit subframe numbering (which counts
parity bits, BDS-SIS-ICD-B1I-3.0 Figures 5-8..5-14) to the corresponding
position in the 224-bit content word. Word 1 contributes bits 1-26 verbatim;
word w ≥ 2 contributes its 22 information bits. Positions inside a parity
field have no content counterpart and are rejected.
"""
function dnav_content_position(icd_position::Int)
    word = (icd_position - 1) ÷ 30 + 1
    offset = icd_position - (word - 1) * 30
    max_info = word == 1 ? 26 : 22
    1 <= offset <= max_info ||
        throw(ArgumentError("ICD bit $icd_position lies in a parity field"))
    word == 1 ? offset : 26 + (word - 2) * 22 + offset
end

# Field accessors on the content word. `start` is the position in the ICD's
# 300-bit numbering of the field's first (most significant) bit; thanks to the
# content mapping the whole field is contiguous even when the figures split it
# across parity boundaries.
dnav_bits(content, start::Int, len::Int) =
    get_bits(content, DNAV_CONTENT_BITS, dnav_content_position(start), len)
dnav_signed(content, start::Int, len::Int) =
    get_twos_complement_num(content, DNAV_CONTENT_BITS, dnav_content_position(start), len)
dnav_bit(content, start::Int) =
    get_bit(content, DNAV_CONTENT_BITS, dnav_content_position(start))

# ---- Decoded data containers --------------------------------------------------

"""
    BeiDouDNAVAlmanac

Almanac data for one BeiDou satellite, decoded from a single D1 almanac page
(subframe 4 pages 1-24 for SV 1-24, subframe 5 pages 1-6 for SV 25-30, and —
when the expanded-almanac identification `AmEpID` is `11` — subframe 5 pages
11-23 for SV 31-63 by time sharing, BDS-SIS-ICD-B1I-3.0 §5.2.4.13-§5.2.4.15).

Angles are stored in radians (the broadcast semicircle values scaled by the
ICD's π). The reference inclination the broadcast `δi` corrects is
i₀ = 0.3 semicircles for MEO/IGSO satellites and i₀ = 0 for GEO
(Table 5-15 note).

# Fields

  - `sqrt_A::Float64`: Square root of semi-major axis (√m)
  - `a_0::Float64`: Satellite clock bias (s)
  - `a_1::Float64`: Satellite clock rate (s/s)
  - `Ω_0::Float64`: Longitude of ascending node at reference time (rad)
  - `e::Float64`: Eccentricity (dimensionless)
  - `δi::Float64`: Correction of orbit reference inclination at reference time (rad)
  - `t_oa::Int`: Almanac reference time (s, scale 2¹²), from this page
  - `Ω_dot::Float64`: Rate of right ascension (rad/s)
  - `ω::Float64`: Argument of perigee (rad)
  - `M_0::Float64`: Mean anomaly at reference time (rad)
  - `WN_a::Int64`: Almanac reference week in force when this page was decoded,
    or `nothing` if subframe 5 page 8 had not been seen yet

!!! note "Use this record's own reference epoch"

    `t_oa` and `WN_a` here belong to *this* entry, and `BeiDouDNAVData` also
    carries a `t_oa`/`WN_a` pair — the global one most recently broadcast in
    subframe 5 page 8. The almanac cycle is 24 pages spread over 12 minutes,
    so the two can disagree across an almanac changeover; pairing an entry
    with the global epoch would then propagate the wrong reference time.
    Prefer the entry's own fields, as [`BeiDouReducedAlmanac`](@ref) forces
    by construction.

# Reference

BDS-SIS-ICD-B1I-3.0, Tables 5-12 and 5-14
"""
Base.@kwdef struct BeiDouDNAVAlmanac
    sqrt_A::Union{Nothing,Float64} = nothing
    a_0::Union{Nothing,Float64} = nothing
    a_1::Union{Nothing,Float64} = nothing
    Ω_0::Union{Nothing,Float64} = nothing
    e::Union{Nothing,Float64} = nothing
    δi::Union{Nothing,Float64} = nothing
    t_oa::Union{Nothing,Int} = nothing
    Ω_dot::Union{Nothing,Float64} = nothing
    ω::Union{Nothing,Float64} = nothing
    M_0::Union{Nothing,Float64} = nothing
    WN_a::Union{Nothing,Int64} = nothing
end

"""
    BeiDouDNAVData

Decoded BeiDou D1/D2 legacy navigation message data, shared by the B1I and
B3I decoders (the broadcast message is structurally identical on both
signals; see `src/beidou/b1i.jl` for the group-delay semantics that differ).

All parameters conform to BDS-SIS-ICD-B1I-3.0 §5.2 (D1) / §5.3 (D2). Angles
broadcast in semicircles are stored in radians (scaled by the ICD π); times
are in seconds of BeiDou Time (BDT).

# Frame Fields

  - `last_subframe_id::Int`: FraID of the last decoded subframe (1-5)
  - `SOW::Int64`: Seconds of week at the leading edge of the current subframe's
    preamble (D1) or of subframe 1 of the current frame (D2), §5.2.4.3/§5.3.3.1
  - `num_bits_after_valid_syncro_sequence_after_last_SOW::Int`: Symbol-counter
    value when `SOW` was decoded (drives the SOW plausibility screen)

# Subframe 1 (D1) / Subframe 1 Pages 1-2 (D2) - Clock, Health, Iono

  - `sat_h1::Bool`: Autonomous satellite health flag (0 = good, §5.2.4.6)
  - `AODC::Int64`: Age of data, clock (§5.2.4.8)
  - `urai::Int64`: User range accuracy index (0-15, §5.2.4.5)
  - `ura::Float64`: User range accuracy (m); `nothing` while URAI = 15 (no prediction)
  - `WN::Int64`: BDT week number (0-8191, weeks since 2006-01-01, §5.2.4.4)
  - `t_0c::Int64`: Clock correction reference time (s, scale 2³)
  - `a_0::Float64`: Clock bias (s, scale 2⁻³³)
  - `a_1::Float64`: Clock rate (s/s, scale 2⁻⁵⁰)
  - `a_2::Float64`: Clock drift rate (s/s², scale 2⁻⁶⁶)
  - `T_GD1::Float64`: B1I equipment group delay differential (s, broadcast in 0.1 ns)
  - `T_GD2::Float64`: B2I equipment group delay differential (s, broadcast in 0.1 ns)
  - `α_0..α_3, β_0..β_3::Float64`: Klobuchar ionospheric model parameters (§5.2.4.7)
  - `AODE::Int64`: Age of data, ephemeris (§5.2.4.11)

# Subframes 2-3 (D1) / Subframe 1 Pages 3-10 (D2) - Ephemeris

  - `t_0e::Int64`: Ephemeris reference time (s, scale 2³; split across
    subframes 2 and 3 in D1 — assembled once both parts are present)
  - `sqrt_A::Float64`: Square root of semi-major axis (√m, scale 2⁻¹⁹)
  - `e::Float64`: Eccentricity (scale 2⁻³³)
  - `ω::Float64`: Argument of perigee (rad)
  - `Δn::Float64`: Mean motion difference (rad/s)
  - `M_0::Float64`: Mean anomaly at reference time (rad)
  - `Ω_0::Float64`: Longitude of ascending node at weekly epoch (rad)
  - `Ω_dot::Float64`: Rate of right ascension (rad/s)
  - `i_0::Float64`: Inclination at reference time (rad)
  - `i_dot::Float64`: Rate of inclination (IDOT, rad/s)
  - `C_uc, C_us::Float64`: Harmonic corrections to argument of latitude (rad, scale 2⁻³¹)
  - `C_rc, C_rs::Float64`: Harmonic corrections to orbit radius (m, scale 2⁻⁶)
  - `C_ic, C_is::Float64`: Harmonic corrections to inclination (rad, scale 2⁻³¹)

# D1 Subframes 4-5 - Almanac, Health, Time Offsets

  - `almanac::Dictionary{Int,BeiDouDNAVAlmanac}`: Per-SVID almanac (SV 1-30,
    plus SV 31-63 when the expanded almanac is broadcast)
  - `health::Dictionary{Int,UInt16}`: Per-SVID 9-bit satellite health
    information words (Table 5-16; 0 = fully healthy)
  - `AmEpID::Int64`: Identification of expanded almanacs (§5.2.4.14)
  - `WN_a::Int64`: Almanac week number (modulo 256, §5.2.4.16)
  - `t_oa::Int64`: Almanac reference time from subframe 5 page 8 (s, scale 2¹²)
  - `A_0GPS, A_1GPS::Float64`: BDT-GPS time offset (s, s/s; §5.2.4.19, not broadcast temporarily)
  - `A_0Gal, A_1Gal::Float64`: BDT-Galileo time offset (s, s/s; §5.2.4.20)
  - `A_0GLO, A_1GLO::Float64`: BDT-GLONASS time offset (s, s/s; §5.2.4.21)
  - `A_0UTC, A_1UTC::Float64`: BDT-UTC offset polynomial (s, s/s; §5.2.4.18)
  - `Δt_LS, Δt_LSF::Int64`: Leap seconds before/after the new leap second (s)
  - `WN_LSF, DN::Int64`: Week number and day number of the new leap second

# Reference

BDS-SIS-ICD-B1I-3.0 §5.2.4, §5.3.3 (and identically BDS-SIS-ICD-B3I-1.0 §5)
"""
Base.@kwdef struct BeiDouDNAVData <: AbstractBeiDouData
    last_subframe_id::Int = 0
    SOW::Union{Nothing,Int64} = nothing
    num_bits_after_valid_syncro_sequence_after_last_SOW::Union{Nothing,Int} = nothing

    # Fundamental navigation information: health, clock, group delay, iono
    sat_h1::Union{Nothing,Bool} = nothing
    AODC::Union{Nothing,Int64} = nothing
    urai::Union{Nothing,Int64} = nothing
    ura::Union{Nothing,Float64} = nothing
    WN::Union{Nothing,Int64} = nothing
    t_0c::Union{Nothing,Int64} = nothing
    a_0::Union{Nothing,Float64} = nothing
    a_1::Union{Nothing,Float64} = nothing
    a_2::Union{Nothing,Float64} = nothing
    T_GD1::Union{Nothing,Float64} = nothing
    T_GD2::Union{Nothing,Float64} = nothing
    α_0::Union{Nothing,Float64} = nothing
    α_1::Union{Nothing,Float64} = nothing
    α_2::Union{Nothing,Float64} = nothing
    α_3::Union{Nothing,Float64} = nothing
    β_0::Union{Nothing,Float64} = nothing
    β_1::Union{Nothing,Float64} = nothing
    β_2::Union{Nothing,Float64} = nothing
    β_3::Union{Nothing,Float64} = nothing
    AODE::Union{Nothing,Int64} = nothing

    # Ephemeris. In D1 t_0e is split: 2 MSBs in subframe 2, 15 LSBs in
    # subframe 3 (Figures 5-9/5-10); the raw staging halves live here so the
    # assembled `t_0e` only appears once both subframes contributed.
    t_0e_msb2::Union{Nothing,Int64} = nothing
    t_0e_lsb15::Union{Nothing,Int64} = nothing
    t_0e::Union{Nothing,Int64} = nothing
    sqrt_A::Union{Nothing,Float64} = nothing
    e::Union{Nothing,Float64} = nothing
    ω::Union{Nothing,Float64} = nothing
    Δn::Union{Nothing,Float64} = nothing
    M_0::Union{Nothing,Float64} = nothing
    Ω_0::Union{Nothing,Float64} = nothing
    Ω_dot::Union{Nothing,Float64} = nothing
    i_0::Union{Nothing,Float64} = nothing
    i_dot::Union{Nothing,Float64} = nothing
    C_uc::Union{Nothing,Float64} = nothing
    C_us::Union{Nothing,Float64} = nothing
    C_rc::Union{Nothing,Float64} = nothing
    C_rs::Union{Nothing,Float64} = nothing
    C_ic::Union{Nothing,Float64} = nothing
    C_is::Union{Nothing,Float64} = nothing

    # D1 subframes 4/5: almanac, health, time offsets
    almanac::Union{Nothing,Dictionary{Int,BeiDouDNAVAlmanac}} = nothing
    health::Union{Nothing,Dictionary{Int,UInt16}} = nothing
    AmEpID::Union{Nothing,Int64} = nothing
    WN_a::Union{Nothing,Int64} = nothing
    t_oa::Union{Nothing,Int64} = nothing
    A_0GPS::Union{Nothing,Float64} = nothing
    A_1GPS::Union{Nothing,Float64} = nothing
    A_0Gal::Union{Nothing,Float64} = nothing
    A_1Gal::Union{Nothing,Float64} = nothing
    A_0GLO::Union{Nothing,Float64} = nothing
    A_1GLO::Union{Nothing,Float64} = nothing
    A_0UTC::Union{Nothing,Float64} = nothing
    A_1UTC::Union{Nothing,Float64} = nothing
    Δt_LS::Union{Nothing,Int64} = nothing
    Δt_LSF::Union{Nothing,Int64} = nothing
    WN_LSF::Union{Nothing,Int64} = nothing
    DN::Union{Nothing,Int64} = nothing
end

function BeiDouDNAVData(
    data::BeiDouDNAVData;
    last_subframe_id = data.last_subframe_id,
    SOW = data.SOW,
    num_bits_after_valid_syncro_sequence_after_last_SOW = data.num_bits_after_valid_syncro_sequence_after_last_SOW,
    sat_h1 = data.sat_h1,
    AODC = data.AODC,
    urai = data.urai,
    ura = data.ura,
    WN = data.WN,
    t_0c = data.t_0c,
    a_0 = data.a_0,
    a_1 = data.a_1,
    a_2 = data.a_2,
    T_GD1 = data.T_GD1,
    T_GD2 = data.T_GD2,
    α_0 = data.α_0,
    α_1 = data.α_1,
    α_2 = data.α_2,
    α_3 = data.α_3,
    β_0 = data.β_0,
    β_1 = data.β_1,
    β_2 = data.β_2,
    β_3 = data.β_3,
    AODE = data.AODE,
    t_0e_msb2 = data.t_0e_msb2,
    t_0e_lsb15 = data.t_0e_lsb15,
    t_0e = data.t_0e,
    sqrt_A = data.sqrt_A,
    e = data.e,
    ω = data.ω,
    Δn = data.Δn,
    M_0 = data.M_0,
    Ω_0 = data.Ω_0,
    Ω_dot = data.Ω_dot,
    i_0 = data.i_0,
    i_dot = data.i_dot,
    C_uc = data.C_uc,
    C_us = data.C_us,
    C_rc = data.C_rc,
    C_rs = data.C_rs,
    C_ic = data.C_ic,
    C_is = data.C_is,
    almanac = data.almanac,
    health = data.health,
    AmEpID = data.AmEpID,
    WN_a = data.WN_a,
    t_oa = data.t_oa,
    A_0GPS = data.A_0GPS,
    A_1GPS = data.A_1GPS,
    A_0Gal = data.A_0Gal,
    A_1Gal = data.A_1Gal,
    A_0GLO = data.A_0GLO,
    A_1GLO = data.A_1GLO,
    A_0UTC = data.A_0UTC,
    A_1UTC = data.A_1UTC,
    Δt_LS = data.Δt_LS,
    Δt_LSF = data.Δt_LSF,
    WN_LSF = data.WN_LSF,
    DN = data.DN,
)
    BeiDouDNAVData(
        last_subframe_id,
        SOW,
        num_bits_after_valid_syncro_sequence_after_last_SOW,
        sat_h1,
        AODC,
        urai,
        ura,
        WN,
        t_0c,
        a_0,
        a_1,
        a_2,
        T_GD1,
        T_GD2,
        α_0,
        α_1,
        α_2,
        α_3,
        β_0,
        β_1,
        β_2,
        β_3,
        AODE,
        t_0e_msb2,
        t_0e_lsb15,
        t_0e,
        sqrt_A,
        e,
        ω,
        Δn,
        M_0,
        Ω_0,
        Ω_dot,
        i_0,
        i_dot,
        C_uc,
        C_us,
        C_rc,
        C_rs,
        C_ic,
        C_is,
        almanac,
        health,
        AmEpID,
        WN_a,
        t_oa,
        A_0GPS,
        A_1GPS,
        A_0Gal,
        A_1Gal,
        A_0GLO,
        A_1GLO,
        A_0UTC,
        A_1UTC,
        Δt_LS,
        Δt_LSF,
        WN_LSF,
        DN,
    )
end

# The default struct `==` falls back to `===` (reference equality), which fails
# for the mutable `almanac` and `health` `Dictionary` fields even when their
# contents match. Compare field-by-field (mirrors `GPSL1CAData`, the LNAV
# decoder this one is modelled on).
#
# This is equality of the *published* data only. Subframe voting compares
# candidates with `dnav_compare_data` below, which weighs the ephemeris fields
# the vote is about rather than every field, and is deliberately not this.
Base.:(==)(a::BeiDouDNAVData, b::BeiDouDNAVData) = fields_equal(a, b)

struct VotedBeiDouDNAVData
    vote::Int
    data::BeiDouDNAVData
end

# Needed on top of the `BeiDouDNAVData` method above, not implied by it: Julia's
# default `==` for a struct is `===`, which compares fields with `===` rather
# than dispatching to their `==`. So without this the wrapper — and hence
# `BeiDouDNAVCache`'s `old_data` vector, and hence the whole decoder state —
# would still compare by reference.
Base.:(==)(a::VotedBeiDouDNAVData, b::VotedBeiDouDNAVData) = fields_equal(a, b)

# One collected D2 subframe-1 page: its SOW (frame epoch) and content word.
struct BeiDouD2Page
    sow::Int64
    content::UInt320
end

Base.:(==)(a::BeiDouD2Page, b::BeiDouD2Page) = a.sow == b.sow && a.content == b.content

"""
$(TYPEDEF)

Per-decoder cache for the BeiDou B1I/B3I legacy signals.

Holds the soft-symbol `CircularDeque{Float32}` (capacity = 300 + 11 = 311),
the data-voting tally used by `confirm_data` (mirroring GPS L1 C/A), and —
for GEO satellites — the D2 subframe-1 page collection: the fundamental
navigation data of D2 is spread over pages 1-10 (one page per 3 s frame,
BDS-SIS-ICD-B1I-3.0 §5.3.2), so decoded pages are staged here until all ten
are present with a consistent SOW chain. Following the framework convention,
only the soft buffer is mutated in place; tally and pages are rebuilt
immutably and threaded through a new cache.

# Fields

$(TYPEDFIELDS)
"""
struct BeiDouDNAVCache <: AbstractGNSSCache
    """
    Soft-symbol buffer (`DNAV_WINDOW_BITS` = 300 syncro + 11 preamble)
    """
    soft_buffer::CircularDeque{Float32}
    """
    Voting tally used by `confirm_data` for subframe-level data validation
    """
    old_data::Vector{VotedBeiDouDNAVData}
    """
    D2 only: subframe-1 pages collected so far, keyed by Pnum1 (1-10)
    """
    d2_pages::Dictionary{Int,BeiDouD2Page}
    """
    Whether the subframe just decoded delivered a SOW that passed the screen
    (see `is_dnav_SOW_from_this_subframe`)
    """
    sow_is_fresh::Bool
end

function BeiDouDNAVCache()
    BeiDouDNAVCache(
        CircularDeque{Float32}(DNAV_WINDOW_BITS),
        Vector{VotedBeiDouDNAVData}(),
        Dictionary{Int,BeiDouD2Page}(),
        false,
    )
end

function BeiDouDNAVCache(
    cache::BeiDouDNAVCache;
    soft_buffer = cache.soft_buffer,
    old_data = cache.old_data,
    d2_pages = cache.d2_pages,
    sow_is_fresh = cache.sow_is_fresh,
)
    BeiDouDNAVCache(soft_buffer, old_data, d2_pages, sow_is_fresh)
end

function Base.:(==)(a::BeiDouDNAVCache, b::BeiDouDNAVCache)
    deques_equal(a.soft_buffer, b.soft_buffer) &&
        a.old_data == b.old_data &&
        a.d2_pages == b.d2_pages
end

packed_buffer_type(::GNSSDecoderState{<:BeiDouDNAVData}) = UInt320

# ---- Completeness checks ------------------------------------------------------

# `T_GD1` is required and `T_GD2` is not, which is the package-wide group-delay
# rule (`is_decoding_completed_for_positioning` in src/gnss.jl) read onto this
# message. The broadcast clock polynomial is referenced to B3I on both signals
# (§5.2.4.10), so a B1I user must apply `T_GD1` to reach its own signal — a
# metre-scale bias, not a refinement — and it rides in the same subframe 1
# (D1) / page 1 (D2) as the clock this gate already waits for, so requiring it
# costs no time to first fix. `T_GD2` is the *B2I* differential: the field is
# decoded and published like any other, but never gated on — the same treatment
# B1C's `T_GD_B2ap` gets — because nothing in this stack can consume it. There is
# no B2I signal in GNSSSignals and none is planned
# (JuliaGNSS/GNSSSignals.jl#156, closed won't-implement): BDS-2 retired 13 of
# its last 15 satellites in April 2026, leaving two IGSO transmitters and no
# successor carrying the signal, so a decoded `T_GD2` has no range to correct.
#
# A B3I user needs neither (its clock is already B3I-referenced), but B1I and
# B3I share one `BeiDouDNAVData` and this check dispatches on the data, so one
# gate serves both. Requiring `T_GD1` on B3I too is free — same subframe — and
# is the conservative side of the trade: the alternative would let a B1I
# consumer see a "ready" decoder with no group delay at all.
function is_dnav_clock_and_health_decoded(data::BeiDouDNAVData)
    !isnothing(data.sat_h1) &&
        !isnothing(data.AODC) &&
        !isnothing(data.urai) &&
        !isnothing(data.WN) &&
        !isnothing(data.t_0c) &&
        !isnothing(data.a_0) &&
        !isnothing(data.a_1) &&
        !isnothing(data.a_2) &&
        !isnothing(data.T_GD1) &&
        !isnothing(data.AODE)
end

function is_dnav_ephemeris_decoded(data::BeiDouDNAVData)
    !isnothing(data.t_0e) &&
        !isnothing(data.sqrt_A) &&
        !isnothing(data.e) &&
        !isnothing(data.ω) &&
        !isnothing(data.Δn) &&
        !isnothing(data.M_0) &&
        !isnothing(data.Ω_0) &&
        !isnothing(data.Ω_dot) &&
        !isnothing(data.i_0) &&
        !isnothing(data.i_dot) &&
        !isnothing(data.C_uc) &&
        !isnothing(data.C_us) &&
        !isnothing(data.C_rc) &&
        !isnothing(data.C_rs) &&
        !isnothing(data.C_ic) &&
        !isnothing(data.C_is)
end

"""
$(TYPEDSIGNATURES)

Check if the BeiDou satellite is healthy and usable for positioning, from the
legacy D1/D2 navigation message.

Examines the autonomous satellite health flag (SatH1) broadcast in D1
subframe 1 / D2 subframe 1 page 1: `0` means the broadcasting satellite is
good, `1` means not (BDS-SIS-ICD-B1I-3.0 / -B3I-1.0 §5.2.4.6).

One method covers both B1I and B3I: the two signals carry the same SatH1 bit
of the same message, so unlike GPS L5I and L2CM — which select *different*
health bits out of a shared CNAV container — there is nothing per-signal to
dispatch on.

!!! warning

    Requires the fundamental navigation data to have been decoded and
    validated; returns `false` until then.

# Arguments

  - `state::GNSSDecoderState{<:BeiDouDNAVData}`: BeiDou B1I or B3I decoder state

# Returns

  - `Bool`: `true` iff SatH1 indicates a healthy satellite
"""
function is_sat_healthy(state::GNSSDecoderState{<:BeiDouDNAVData})
    state.data.sat_h1 === false
end

function is_decoding_completed_for_positioning(data::BeiDouDNAVData)
    # Klobuchar ionosphere is deliberately *not* gated on: `src/gnss.jl` states
    # that second-order corrections are excluded from readiness and should be
    # applied when present rather than waited for. It happens to ride in the
    # same subframe as the clock set here (D1 subframe 1 / D2 page 2), so this
    # is a contract fix rather than a behaviour change — but the contract is
    # what other consumers rely on.
    !isnothing(data.SOW) &&
        is_dnav_clock_and_health_decoded(data) &&
        is_dnav_ephemeris_decoded(data)
end

# ---- SOW plausibility screen ---------------------------------------------------
#
# Frame sync is only an 11-bit preamble match, so — exactly as for GPS L1 C/A
# (see `is_plausible_TOW`) — a false lock can reach the SOW field with an
# arbitrary count behind valid BCH parity. The same two screens apply: once
# the symbol counter runs across two SOW decodes, the elapsed symbols predict
# the next SOW exactly (D1: 50 sps and a 300-symbol subframe grid; D2: 500 sps
# and SOW readings from subframe 1 only, 1500 symbols apart); before that, a
# bounded forward step.

const DNAV_MAX_SOW_GAP = 600

function is_plausible_dnav_SOW(
    SOW_count,
    prev_SOW,
    prev_SOW_anchor,
    num_bits_after_valid_syncro_sequence,
    symbols_per_second,
)
    SOW_count < SECONDS_PER_WEEK || return false
    isnothing(prev_SOW) && return true
    if !isnothing(prev_SOW_anchor) && !isnothing(num_bits_after_valid_syncro_sequence)
        elapsed = num_bits_after_valid_syncro_sequence - prev_SOW_anchor
        (elapsed > 0 && elapsed % 300 == 0 && elapsed % symbols_per_second == 0) ||
            return false
        return Int64(SOW_count) ==
               mod(prev_SOW + elapsed ÷ symbols_per_second, SECONDS_PER_WEEK)
    end
    ΔSOW = mod(Int64(SOW_count) - prev_SOW, SECONDS_PER_WEEK)
    return 0 < ΔSOW <= DNAV_MAX_SOW_GAP
end

# D1 runs at 50 symbols per second, D2 (GEO satellites) at 500 (§5.1.1). The
# format is a property of the satellite, not of the signal, so it is selected
# per PRN.
dnav_symbols_per_second(prn::Integer) = is_beidou_geo(prn) ? 500 : 50

"""
    GNSSSignals.get_data_frequency(state::GNSSDecoderState{<:BeiDouDNAVData})

Navigation-message symbol rate of a B1I/B3I decoder state: 50 Hz for the
MEO/IGSO satellites (D1 NAV) and 500 Hz for the GEO satellites (D2 NAV,
BDS-SIS-ICD-B1I-3.0 §5.1.1).

This specializes the generic forward through `get_signal_type` (see
`src/gnss.jl`): the type-level `GNSSSignals.get_data_frequency(BeiDouB1I)`
answers 50 Hz because a *signal type* cannot know the satellite, but the D1/D2
rate is a property of the PRN — which the decoder state does know — so the
state-level accessor reports the rate this state actually demodulates.
"""
GNSSSignals.get_data_frequency(state::GNSSDecoderState{<:BeiDouDNAVData}) =
    dnav_symbols_per_second(state.prn) * Hz

# Screen and store a freshly decoded SOW (with its symbol-counter anchor) in
# `raw_data`.
#
# A SOW that fails the screen leaves the previously accepted (SOW, anchor)
# pair in place rather than clearing it. Clearing would defeat the screen
# itself: `is_plausible_dnav_SOW` waves through any SOW once `prev_SOW` is
# `nothing`, so a rejection would leave the very next subframe unscreened.
# Keeping the pair costs nothing — the elapsed-symbol prediction spans however
# many subframes were skipped, since `elapsed % 300 == 0` holds for any whole
# number of them.
#
# Freshness is therefore not carried by `SOW` being non-`nothing` but by the
# `sow_is_fresh` cache flag this raises on acceptance — see
# `is_dnav_SOW_from_this_subframe`.
function store_dnav_SOW(state::GNSSDecoderState{<:BeiDouDNAVData}, content)
    SOW_count = Int64(dnav_bits(content, 19, 20))
    prev_SOW = state.raw_data.SOW
    prev_anchor = state.raw_data.num_bits_after_valid_syncro_sequence_after_last_SOW
    num_bits = state.num_bits_after_valid_syncro_sequence
    is_plausible = is_plausible_dnav_SOW(
        SOW_count,
        prev_SOW,
        prev_anchor,
        num_bits,
        dnav_symbols_per_second(state.prn),
    )
    is_plausible || return state
    GNSSDecoderState(
        state;
        raw_data = BeiDouDNAVData(
            state.raw_data;
            SOW = SOW_count,
            num_bits_after_valid_syncro_sequence_after_last_SOW = num_bits,
        ),
        cache = BeiDouDNAVCache(state.cache; sow_is_fresh = true),
    )
end

"""
True when `raw_data.SOW` was decoded from the subframe that just synced.

`decode_syncro_sequence` clears the flag as it starts and `store_dnav_SOW`
sets it only on acceptance, so it distinguishes a SOW read *now* from one
carried over. The symbol counter cannot serve here: it stays `nothing` until
the first promotion, so an anchor-versus-counter comparison would read
`nothing == nothing` and call every subframe fresh.

Two paths depend on this. Promotion re-anchors the counter to the SOW epoch,
which is only meaningful for a SOW read from this very subframe: D2 reads a
SOW in subframe 1 only (§5.3.3.1(2)), so subframes 2-5 must not re-anchor a
SOW up to four subframes old. And because BCH(15,11,1) is a perfect code —
`dnav_bch_decode` always returns *some* codeword and can never reject — the
SOW screen is the only integrity gate this message has, so a subframe whose
SOW failed it must not write its payload fields either.
"""
is_dnav_SOW_from_this_subframe(state::GNSSDecoderState{<:BeiDouDNAVData}) =
    state.cache.sow_is_fresh

# ---- D1 subframe parsers -------------------------------------------------------
#
# Field positions cite the first (most significant) bit in the ICD's 300-bit
# numbering (BDS-SIS-ICD-B1I-3.0 Figures 5-8, 5-9, 5-10, 5-11-x); scale
# factors cite Tables 5-4..5-18. Parity-split fields are contiguous in the
# content domain, so a single (start, length) pair reads the whole field.

# URAI → URA in meters (Table 5-4): N < 6 ⇒ 2^(N/2+1) with the ICD's rounding
# for odd N; 6 ≤ N < 15 ⇒ 2^(N-2); N = 15 ⇒ no accuracy prediction (`nothing`).
function dnav_ura_from_urai(urai::Integer)
    urai == 15 && return nothing
    urai == 1 && return 2.8
    urai == 3 && return 5.7
    urai == 5 && return 11.3
    urai < 6 ? 2.0^(urai / 2 + 1) : 2.0^(urai - 2)
end

function decode_d1_subframe1(state::GNSSDecoderState{<:BeiDouDNAVData}, content)
    urai = Int64(dnav_bits(content, 49, 4))
    data = BeiDouDNAVData(
        state.raw_data;
        sat_h1 = dnav_bit(content, 43),
        AODC = Int64(dnav_bits(content, 44, 5)),
        urai,
        ura = dnav_ura_from_urai(urai),
        WN = Int64(dnav_bits(content, 61, 13)),
        t_0c = Int64(dnav_bits(content, 74, 17)) << 3,
        T_GD1 = dnav_signed(content, 99, 10) * 1.0e-10,   # 0.1 ns → s
        T_GD2 = dnav_signed(content, 109, 10) * 1.0e-10,  # 0.1 ns → s
        α_0 = dnav_signed(content, 127, 8) / (1 << 30),
        α_1 = dnav_signed(content, 135, 8) / (1 << 27),
        α_2 = dnav_signed(content, 151, 8) / (1 << 24),
        α_3 = dnav_signed(content, 159, 8) / (1 << 24),
        β_0 = dnav_signed(content, 167, 8) * Float64(1 << 11),
        β_1 = dnav_signed(content, 183, 8) * Float64(1 << 14),
        β_2 = dnav_signed(content, 191, 8) * Float64(1 << 16),
        β_3 = dnav_signed(content, 199, 8) * Float64(1 << 16),
        a_2 = dnav_signed(content, 215, 11) / 2.0^66,
        a_0 = dnav_signed(content, 226, 24) / 2.0^33,
        a_1 = dnav_signed(content, 258, 22) / 2.0^50,
        AODE = Int64(dnav_bits(content, 288, 5)),
    )
    GNSSDecoderState(state; raw_data = data)
end

# Assemble the 17-bit t_0e (× 2³ s) from its subframe-2 MSBs and subframe-3
# LSBs, but only when the two halves provably belong to one ephemeris.
#
# D1 splits t_0e across two subframes and, having no issue-of-data stamp to
# pair them with, offers only their adjacency: subframe 3 follows subframe 2
# 6 s later within the same frame, and the ICD changes the ephemeris at frame
# boundaries. So the halves match iff subframe 2 was the immediately preceding
# decoded subframe. If subframe 2 was missed or failed its SOW screen, the
# stored MSBs belong to an earlier frame and pairing them with fresh LSBs
# would fabricate a t_0e that is in neither ephemeris — plausible enough to
# key the voting dataset and pass every other check.
function with_assembled_t_0e(data::BeiDouDNAVData, prev_subframe_id::Int)
    prev_subframe_id == 2 || return data
    (isnothing(data.t_0e_msb2) || isnothing(data.t_0e_lsb15)) && return data
    BeiDouDNAVData(data; t_0e = (data.t_0e_msb2 << 15 | data.t_0e_lsb15) << 3)
end

function decode_d1_subframe2(state::GNSSDecoderState{<:BeiDouDNAVData}, content)
    PI = state.constants.PI
    data = BeiDouDNAVData(
        state.raw_data;
        Δn = dnav_signed(content, 43, 16) * PI / 2.0^43,
        C_uc = dnav_signed(content, 67, 18) / 2.0^31,
        M_0 = dnav_signed(content, 93, 32) * PI / 2.0^31,
        e = Int64(dnav_bits(content, 133, 32)) / 2.0^33,
        C_us = dnav_signed(content, 181, 18) / 2.0^31,
        C_rc = dnav_signed(content, 199, 18) / (1 << 6),
        C_rs = dnav_signed(content, 225, 18) / (1 << 6),
        sqrt_A = Int64(dnav_bits(content, 251, 32)) / 2.0^19,
        t_0e_msb2 = Int64(dnav_bits(content, 291, 2)),
        # These MSBs open a new pairing: the LSBs still held here are the
        # previous frame's, so drop them (and the t_0e they assembled) rather
        # than let subframe 3 pair across the frame boundary.
        t_0e_lsb15 = nothing,
        t_0e = nothing,
    )
    GNSSDecoderState(state; raw_data = data)
end

function decode_d1_subframe3(state::GNSSDecoderState{<:BeiDouDNAVData}, content)
    PI = state.constants.PI
    data = BeiDouDNAVData(
        state.raw_data;
        # Bits 3-17 of t_0e: the 10 middle bits (subframe 3 word 2) and the 5
        # LSBs (word 3) are contiguous in the content domain (Figure 5-10).
        t_0e_lsb15 = Int64(dnav_bits(content, 43, 15)),
        i_0 = dnav_signed(content, 66, 32) * PI / 2.0^31,
        C_ic = dnav_signed(content, 106, 18) / 2.0^31,
        Ω_dot = dnav_signed(content, 132, 24) * PI / 2.0^43,
        C_is = dnav_signed(content, 164, 18) / 2.0^31,
        i_dot = dnav_signed(content, 190, 14) * PI / 2.0^43,
        Ω_0 = dnav_signed(content, 212, 32) * PI / 2.0^31,
        ω = dnav_signed(content, 252, 32) * PI / 2.0^31,
    )
    GNSSDecoderState(
        state;
        raw_data = with_assembled_t_0e(data, state.raw_data.last_subframe_id),
    )
end

# ---- D1 subframe 4/5 page parsers ----------------------------------------------

# Almanac page layout shared by subframe 4 pages 1-24, subframe 5 pages 1-6,
# and subframe 5 pages 11-23 (Figures 5-11-1 / 5-11-6); trailing 2-bit field
# is AmEpID on the former two, AmID on the expanded pages.
function decode_d1_almanac_page(state::GNSSDecoderState{<:BeiDouDNAVData}, content, sv_id)
    PI = state.constants.PI
    sqrt_A_raw = Int64(dnav_bits(content, 51, 24))
    # An all-zero √A marks an empty/dummy almanac slot (no satellite has a
    # zero semi-major axis); skip it rather than store zeros.
    sqrt_A_raw == 0 && return state
    entry = BeiDouDNAVAlmanac(;
        sqrt_A = sqrt_A_raw / (1 << 11),
        a_1 = dnav_signed(content, 91, 11) / 2.0^38,
        a_0 = dnav_signed(content, 102, 11) / (1 << 20),
        Ω_0 = dnav_signed(content, 121, 24) * PI / (1 << 23),
        e = Int64(dnav_bits(content, 153, 17)) / (1 << 21),
        δi = dnav_signed(content, 170, 16) * PI / (1 << 19),
        t_oa = Int64(dnav_bits(content, 194, 8)) << 12,
        Ω_dot = dnav_signed(content, 202, 17) * PI / 2.0^38,
        ω = dnav_signed(content, 227, 24) * PI / (1 << 23),
        M_0 = dnav_signed(content, 259, 24) * PI / (1 << 23),
        # Snapshot the reference week in force now: the page carries its own
        # t_oa but no week, and the global `WN_a` moves on at the next
        # almanac changeover.
        WN_a = state.raw_data.WN_a,
    )
    almanac = _merge_keyed(state.raw_data.almanac, Int(sv_id), entry)
    GNSSDecoderState(state; raw_data = BeiDouDNAVData(state.raw_data; almanac))
end

# Store `count` consecutive 9-bit health words (Table 5-16) starting at the
# content position of ICD bit 51, for SV IDs `first_sv_id .. first_sv_id+count-1`.
function decode_d1_health_page(
    state::GNSSDecoderState{<:BeiDouDNAVData},
    content,
    first_sv_id::Int,
    count::Int,
)
    health =
        isnothing(state.raw_data.health) ? Dictionary{Int,UInt16}() :
        copy(state.raw_data.health)
    start = dnav_content_position(51)
    for k = 0:(count-1)
        code = UInt16(get_bits(content, DNAV_CONTENT_BITS, start + 9k, 9))
        set!(health, Int(first_sv_id + k), code)
    end
    GNSSDecoderState(state; raw_data = BeiDouDNAVData(state.raw_data; health))
end

# Record the AmEpID broadcast at the tail of every subframe 4 page and
# subframe 5 pages 1-6 (Figure 5-11-1); it gates the expanded-almanac pages.
function with_am_ep_id(state::GNSSDecoderState{<:BeiDouDNAVData}, content)
    GNSSDecoderState(
        state;
        raw_data = BeiDouDNAVData(
            state.raw_data;
            AmEpID = Int64(dnav_bits(content, 291, 2)),
        ),
    )
end

# Map an expanded-almanac page (subframe 5 pages 11-23, valid only while
# AmEpID = 11) to the SV ID it carries, per the AmID time-sharing scheme of
# Table 5-13; `nothing` for reserved combinations.
function expanded_almanac_sv_id(am_id, pnum)
    am_id == 0b01 && return 31 + (pnum - 11)
    am_id == 0b10 && return 44 + (pnum - 11)
    am_id == 0b11 && return pnum <= 17 ? 57 + (pnum - 11) : nothing
    nothing
end

function decode_d1_subframe5_page7(state, content)
    decode_d1_health_page(state, content, 1, 19)
end

function decode_d1_subframe5_page8(state, content)
    state = decode_d1_health_page(state, content, 20, 11)
    data = BeiDouDNAVData(
        state.raw_data;
        WN_a = Int64(dnav_bits(content, 190, 8)),
        # toa here is split 5 MSBs / 3 LSBs across words 7/8 (Figure 5-11-3),
        # contiguous in the content domain.
        t_oa = Int64(dnav_bits(content, 198, 8)) << 12,
    )
    GNSSDecoderState(state; raw_data = data)
end

function decode_d1_subframe5_page9(state, content)
    data = BeiDouDNAVData(
        state.raw_data;
        A_0GPS = dnav_signed(content, 97, 14) * 1.0e-10,   # 0.1 ns → s
        A_1GPS = dnav_signed(content, 111, 16) * 1.0e-10,  # 0.1 ns/s → s/s
        A_0Gal = dnav_signed(content, 135, 14) * 1.0e-10,
        A_1Gal = dnav_signed(content, 157, 16) * 1.0e-10,
        A_0GLO = dnav_signed(content, 181, 14) * 1.0e-10,
        A_1GLO = dnav_signed(content, 195, 16) * 1.0e-10,
    )
    GNSSDecoderState(state; raw_data = data)
end

function decode_d1_subframe5_page10(state, content)
    data = BeiDouDNAVData(
        state.raw_data;
        Δt_LS = Int64(dnav_signed(content, 51, 8)),
        Δt_LSF = Int64(dnav_signed(content, 67, 8)),
        WN_LSF = Int64(dnav_bits(content, 75, 8)),
        A_0UTC = dnav_signed(content, 91, 32) / 2.0^30,
        A_1UTC = dnav_signed(content, 131, 24) / 2.0^50,
        DN = Int64(dnav_bits(content, 163, 8)),
    )
    GNSSDecoderState(state; raw_data = data)
end

function decode_d1_subframe4(state::GNSSDecoderState{<:BeiDouDNAVData}, content)
    pnum = Int(dnav_bits(content, 44, 7))
    pnum in 1:24 || return state
    # Subframe 4 pages 1-24 carry the almanac of SV 1-24 (§5.2.4.13) and the
    # AmEpID flag (Figure 5-11-1).
    state = with_am_ep_id(state, content)
    decode_d1_almanac_page(state, content, pnum)
end

function decode_d1_subframe5(state::GNSSDecoderState{<:BeiDouDNAVData}, content)
    pnum = Int(dnav_bits(content, 44, 7))
    if pnum in 1:6
        # Almanac of SV 25-30, plus AmEpID (§5.2.4.13, Figure 5-11-1).
        state = with_am_ep_id(state, content)
        state = decode_d1_almanac_page(state, content, 24 + pnum)
    elseif pnum == 7
        state = decode_d1_subframe5_page7(state, content)
    elseif pnum == 8
        state = decode_d1_subframe5_page8(state, content)
    elseif pnum == 9
        state = decode_d1_subframe5_page9(state, content)
    elseif pnum == 10
        state = decode_d1_subframe5_page10(state, content)
    elseif pnum in 11:23
        # Expanded almanac for SV 31-63, valid only while AmEpID = 11
        # (§5.2.4.14); AmID (bits 291-292 of these pages) selects the SV group.
        if state.raw_data.AmEpID == 0b11
            am_id = Int(dnav_bits(content, 291, 2))
            sv_id = expanded_almanac_sv_id(am_id, pnum)
            if !isnothing(sv_id)
                state = decode_d1_almanac_page(state, content, sv_id)
            end
        end
    elseif pnum == 24
        # Expanded health for SV 31-63, same AmEpID/AmID gating (Table 5-17).
        if state.raw_data.AmEpID == 0b11
            am_id = Int(dnav_bits(content, 216, 2))
            if am_id == 0b01
                state = decode_d1_health_page(state, content, 31, 13)
            elseif am_id == 0b10
                state = decode_d1_health_page(state, content, 44, 13)
            elseif am_id == 0b11
                state = decode_d1_health_page(state, content, 57, 7)
            end
        end
    end
    state
end

# ---- D2 subframe-1 page collection and parsing ----------------------------------
#
# D2 spreads the fundamental navigation data over pages 1-10 of subframe 1,
# one page per 3-second frame (Figures 5-14-1..-10); most ephemeris fields
# are split across pages. Decoded pages are staged in the cache keyed by
# Pnum1 together with their frame SOW; once all ten are present with the
# expected 3-seconds-per-page SOW chain (i.e. they belong to one 30-second
# broadcast cycle), the complete field set is assembled into `raw_data`.

# Content positions (ICD 300-bit numbering) of the per-page fields. All D2
# scale factors equal their D1 counterparts (§5.3.3.1).
function parse_d2_pages(state::GNSSDecoderState{<:BeiDouDNAVData}, pages)
    PI = state.constants.PI
    c(p) = pages[p].content
    urai = Int64(dnav_bits(c(1), 61, 4))
    # Split fields, assembled MSBs-first from their per-page pieces before
    # sign extension (positions per Figures 5-14-3..-10).
    a_1_raw =
        dnav_bits(c(3), 133, 4) << 18 | dnav_bits(c(4), 47, 6) << 12 |
        dnav_bits(c(4), 61, 12)
    a_2_raw = dnav_bits(c(4), 73, 10) << 1 | dnav_bits(c(4), 91, 1)
    C_uc_raw = dnav_bits(c(4), 121, 14) << 4 | dnav_bits(c(5), 47, 4)
    C_us_raw = dnav_bits(c(5), 99, 14) << 4 | dnav_bits(c(5), 121, 4)
    e_raw =
        dnav_bits(c(5), 125, 10) << 22 | dnav_bits(c(6), 47, 6) << 16 |
        dnav_bits(c(6), 61, 16)
    sqrt_A_raw = dnav_bits(c(6), 77, 6) << 26 | dnav_bits(c(6), 91, 26)
    C_ic_raw =
        dnav_bits(c(6), 125, 10) << 8 | dnav_bits(c(7), 47, 6) << 2 | dnav_bits(c(7), 61, 2)
    t_0e_raw = dnav_bits(c(7), 81, 2) << 15 | dnav_bits(c(7), 91, 15)
    i_0_raw =
        dnav_bits(c(7), 106, 7) << 25 | dnav_bits(c(7), 121, 14) << 11 |
        dnav_bits(c(8), 47, 6) << 5 | dnav_bits(c(8), 61, 5)
    C_rc_raw = dnav_bits(c(8), 66, 17) << 1 | dnav_bits(c(8), 91, 1)
    Ω_dot_raw =
        dnav_bits(c(8), 110, 3) << 21 | dnav_bits(c(8), 121, 16) << 5 |
        dnav_bits(c(9), 47, 5)
    Ω_0_raw =
        dnav_bits(c(9), 52, 1) << 31 | dnav_bits(c(9), 61, 22) << 9 | dnav_bits(c(9), 91, 9)
    ω_raw =
        dnav_bits(c(9), 100, 13) << 19 | dnav_bits(c(9), 121, 14) << 5 |
        dnav_bits(c(10), 47, 5)
    i_dot_raw = dnav_bits(c(10), 52, 1) << 13 | dnav_bits(c(10), 61, 13)
    sext(raw, width) = get_twos_complement_num(UInt64(raw), width, 1, width)
    data = BeiDouDNAVData(
        state.raw_data;
        sat_h1 = dnav_bit(c(1), 47),
        AODC = Int64(dnav_bits(c(1), 48, 5)),
        urai,
        ura = dnav_ura_from_urai(urai),
        WN = Int64(dnav_bits(c(1), 65, 13)),
        t_0c = Int64(dnav_bits(c(1), 78, 17)) << 3,
        T_GD1 = dnav_signed(c(1), 103, 10) * 1.0e-10,
        T_GD2 = dnav_signed(c(1), 121, 10) * 1.0e-10,
        α_0 = sext(dnav_bits(c(2), 47, 6) << 2 | dnav_bits(c(2), 61, 2), 8) / (1 << 30),
        α_1 = dnav_signed(c(2), 63, 8) / (1 << 27),
        α_2 = dnav_signed(c(2), 71, 8) / (1 << 24),
        α_3 = sext(dnav_bits(c(2), 79, 4) << 4 | dnav_bits(c(2), 91, 4), 8) / (1 << 24),
        β_0 = dnav_signed(c(2), 95, 8) * Float64(1 << 11),
        β_1 = dnav_signed(c(2), 103, 8) * Float64(1 << 14),
        β_2 = sext(dnav_bits(c(2), 111, 2) << 6 | dnav_bits(c(2), 121, 6), 8) *
              Float64(1 << 16),
        β_3 = dnav_signed(c(2), 127, 8) * Float64(1 << 16),
        a_0 = sext(dnav_bits(c(3), 101, 12) << 12 | dnav_bits(c(3), 121, 12), 24) / 2.0^33,
        a_1 = sext(a_1_raw, 22) / 2.0^50,
        a_2 = sext(a_2_raw, 11) / 2.0^66,
        AODE = Int64(dnav_bits(c(4), 92, 5)),
        Δn = dnav_signed(c(4), 97, 16) * PI / 2.0^43,
        C_uc = sext(C_uc_raw, 18) / 2.0^31,
        M_0 = sext(
            dnav_bits(c(5), 51, 2) << 30 | dnav_bits(c(5), 61, 22) << 8 |
            dnav_bits(c(5), 91, 8),
            32,
        ) * PI / 2.0^31,
        C_us = sext(C_us_raw, 18) / 2.0^31,
        e = Float64(e_raw) / 2.0^33,
        sqrt_A = Float64(sqrt_A_raw) / 2.0^19,
        C_ic = sext(C_ic_raw, 18) / 2.0^31,
        C_is = dnav_signed(c(7), 63, 18) / 2.0^31,
        t_0e = Int64(t_0e_raw) << 3,
        t_0e_msb2 = Int64(t_0e_raw >> 15),
        t_0e_lsb15 = Int64(t_0e_raw & 0x7FFF),
        i_0 = sext(i_0_raw, 32) * PI / 2.0^31,
        C_rc = sext(C_rc_raw, 18) / (1 << 6),
        C_rs = dnav_signed(c(8), 92, 18) / (1 << 6),
        Ω_dot = sext(Ω_dot_raw, 24) * PI / 2.0^43,
        Ω_0 = sext(Ω_0_raw, 32) * PI / 2.0^31,
        ω = sext(ω_raw, 32) * PI / 2.0^31,
        i_dot = sext(i_dot_raw, 14) * PI / 2.0^43,
    )
    GNSSDecoderState(state; raw_data = data)
end

function decode_d2_subframe1(state::GNSSDecoderState{<:BeiDouDNAVData}, content)
    # The staged page's SOW must be this frame's SOW (already screened); an
    # implausible SOW also invalidates the page for the completeness chain.
    sow = state.raw_data.SOW
    isnothing(sow) && return state
    pnum1 = Int(dnav_bits(content, 43, 4))
    pnum1 in 1:10 || return state
    pages = copy(state.cache.d2_pages)
    set!(pages, pnum1, BeiDouD2Page(sow, content))
    state = GNSSDecoderState(state; cache = BeiDouDNAVCache(state.cache; d2_pages = pages))
    # Parse only when all ten pages of one 30-second broadcast cycle are
    # present: page p must carry SOW(page 1) + 3(p-1) (one page per 3 s frame).
    all(p -> haskey(pages, p), 1:10) || return state
    sow1 = pages[1].sow
    all(p -> pages[p].sow == mod(sow1 + 3 * (p - 1), SECONDS_PER_WEEK), 1:10) ||
        return state
    parse_d2_pages(state, pages)
end

# ---- Frame dispatch --------------------------------------------------------------

function decode_syncro_sequence(state::GNSSDecoderState{<:BeiDouDNAVData}, buffer)
    # Clear the freshness flag first: every path out of this function leaves
    # it false unless `store_dnav_SOW` accepts a SOW from this very subframe.
    state =
        GNSSDecoderState(state; cache = BeiDouDNAVCache(state.cache; sow_is_fresh = false))
    content = decode_dnav_subframe_content(buffer, state.constants)
    fra_id = Int(dnav_bits(content, 16, 3))
    fra_id in 1:5 || return state  # FraID 110/111 are reserved (Table 5-3)
    if is_beidou_geo(state.prn)
        # D2: the fundamental navigation data lives in subframe 1 pages 1-10;
        # subframes 2-4 (BDS integrity / wide-area differential corrections)
        # and subframe 5 (D2 almanac & time offsets, 120 pages) belong to the
        # legacy regional augmentation service and are intentionally not
        # decoded (§5.3.3.1). The D2 SOW refers to subframe 1 of the frame
        # (§5.3.3.1 (2)), so it is read only there.
        fra_id == 1 || return state
        state = store_dnav_SOW(state, content)
        # The SOW screen is this message's only integrity gate (see
        # `is_dnav_SOW_from_this_subframe`), so a subframe that failed it
        # contributes nothing — not its payload and not its FraID.
        is_dnav_SOW_from_this_subframe(state) || return state
        state = decode_d2_subframe1(state, content)
        state = GNSSDecoderState(
            state;
            raw_data = BeiDouDNAVData(state.raw_data; last_subframe_id = fra_id),
        )
    else
        # D1: every subframe carries its own SOW (§5.2.4.3).
        state = store_dnav_SOW(state, content)
        is_dnav_SOW_from_this_subframe(state) || return state
        # `last_subframe_id` is updated *after* the parsers run, so that they
        # observe the FraID of the preceding subframe. `decode_d1_subframe3`
        # needs it to prove that subframe 2 came immediately before, which is
        # what makes the two halves of `t_0e` one ephemeris.
        if fra_id == 1
            state = decode_d1_subframe1(state, content)
        elseif fra_id == 2
            state = decode_d1_subframe2(state, content)
        elseif fra_id == 3
            state = decode_d1_subframe3(state, content)
        elseif fra_id == 4
            state = decode_d1_subframe4(state, content)
        elseif fra_id == 5
            state = decode_d1_subframe5(state, content)
        end
        state = GNSSDecoderState(
            state;
            raw_data = BeiDouDNAVData(state.raw_data; last_subframe_id = fra_id),
        )
    end
    state
end

# ---- Data validation (voting) -----------------------------------------------------
#
# The legacy message has no issue-of-data stamps: AODC/AODE are broadcast
# *ages*, not version numbers, so cross-subframe consistency cannot be gated
# on an IODC/IODE match the way GPS L1 C/A does. The dataset key is instead
# the (t_0c, t_0e) reference-time pair — the ICD guarantees both change
# whenever their parameter sets change (§5.2.4.9/§5.2.4.12: "the value of
# toc/toe shall monotonically increase over the week and shall change if any
# of the parameters change") — and the same broadcast-repetition voting as
# GPS L1 C/A guards against BCH(15,11,1) mis-corrections (the code corrects
# one error per 15-bit block; two errors are silently mis-corrected).

function dnav_compare_data(data::BeiDouDNAVData, new_data::BeiDouDNAVData)
    data.WN == new_data.WN &&
        data.t_0c == new_data.t_0c &&
        data.a_0 == new_data.a_0 &&
        data.a_1 == new_data.a_1 &&
        data.a_2 == new_data.a_2 &&
        data.T_GD1 == new_data.T_GD1 &&
        data.T_GD2 == new_data.T_GD2 &&
        data.t_0e == new_data.t_0e &&
        data.sqrt_A == new_data.sqrt_A &&
        data.e == new_data.e &&
        data.ω == new_data.ω &&
        data.Δn == new_data.Δn &&
        data.M_0 == new_data.M_0 &&
        data.Ω_0 == new_data.Ω_0 &&
        data.Ω_dot == new_data.Ω_dot &&
        data.i_0 == new_data.i_0 &&
        data.i_dot == new_data.i_dot &&
        data.C_uc == new_data.C_uc &&
        data.C_us == new_data.C_us &&
        data.C_rc == new_data.C_rc &&
        data.C_rs == new_data.C_rs &&
        data.C_ic == new_data.C_ic &&
        data.C_is == new_data.C_is
end

dnav_dataset_key(data::BeiDouDNAVData) = (data.t_0c, data.t_0e)

# Thread an updated voting tally through a new cache (see GPS L1 C/A's
# `with_old_data`).
dnav_with_old_data(state, new_old_data; kwargs...) = GNSSDecoderState(
    state;
    cache = BeiDouDNAVCache(state.cache; old_data = new_old_data),
    kwargs...,
)

# Promote `raw_data` to validated `data`, re-anchoring the symbol counter to
# the epoch the SOW stamps, so that a consumer reads the current time as
# `SOW + num_bits_after_valid_syncro_sequence / get_data_frequency(state)`.
#
# That epoch is the *leading edge of this subframe's preamble*
# (BDS-SIS-ICD-B1I-3.0 §5.2.4.3 for D1, §5.3.3.1(2) for D2), and at promotion
# time the whole subframe plus the next preamble sit in the window — so the
# anchor is `syncro_sequence_length + preamble_length` (311 symbols; 6.22 s at
# D1's 50 sps, 0.622 s at D2's 500 sps), matching `b2a.jl` and `b2b.jl`.
#
# This is where BeiDou parts company with GPS: LNAV's HOW TOW stamps the
# *next* subframe, so `gps/l1ca.jl` anchors to `preamble_length` alone. Using
# the GPS anchor here would report the time one whole subframe early.
#
# The SOW anchor is rebased into the same counting frame: promotion only
# happens with a SOW decoded in this very subframe (`validate_data` enforces
# that via `is_dnav_SOW_from_this_subframe`), so the anchor being rebased
# equals the counter being replaced.
function dnav_promote_data(state, new_old_data)
    anchor = state.constants.syncro_sequence_length + state.constants.preamble_length
    promoted = BeiDouDNAVData(
        state.raw_data;
        num_bits_after_valid_syncro_sequence_after_last_SOW = anchor,
    )
    dnav_with_old_data(
        state,
        new_old_data;
        raw_data = promoted,
        data = promoted,
        num_bits_after_valid_syncro_sequence = anchor,
    )
end

function dnav_confirm_data(state, max_vote = 20)
    old_data = state.cache.old_data
    key = dnav_dataset_key(state.raw_data)

    matching_idx = findfirst(old_data) do entry
        dnav_dataset_key(entry.data) == key && dnav_compare_data(entry.data, state.raw_data)
    end

    if isnothing(matching_idx)
        # First sighting of this dataset — stage it and wait for the broadcast
        # to repeat it. GPS L1 C/A promotes here on first sight, which is safe
        # only because its word parity detects errors; BCH(15,11,1) is a
        # perfect code, so `dnav_bch_decode` always returns *some* codeword and
        # a two-error word is silently mis-corrected. Repetition voting is the
        # only thing standing between a mis-correction and `state.data`, so
        # first-fix promotion would walk straight past it. The cost is one
        # extra broadcast cycle before the first fix.
        new_old_data = push!(copy(old_data), VotedBeiDouDNAVData(0, state.raw_data))
        return dnav_with_old_data(state, new_old_data; raw_data = BeiDouDNAVData())
    end

    curr_score = old_data[matching_idx].vote
    new_vote = increment_voting(curr_score, max_vote)
    best_score = maximum(e.vote for e in old_data if dnav_dataset_key(e.data) == key)

    if best_score > curr_score
        new_old_data = copy(old_data)
        new_old_data[matching_idx] = VotedBeiDouDNAVData(new_vote, state.raw_data)
        return dnav_with_old_data(state, new_old_data; raw_data = BeiDouDNAVData())
    end

    new_old_data = if new_vote == max_vote && length(old_data) > 1
        [VotedBeiDouDNAVData(new_vote, state.raw_data)]
    else
        updated = copy(old_data)
        updated[matching_idx] = VotedBeiDouDNAVData(new_vote, state.raw_data)
        updated
    end

    dnav_promote_data(state, new_old_data)
end

function validate_data(state::GNSSDecoderState{<:BeiDouDNAVData})
    # Promotion re-anchors the symbol counter to this subframe's SOW epoch, so
    # it may only run for a SOW decoded from the subframe that just synced.
    # `decode` calls this hook after *every* sync, including the D2 subframes
    # 2-5 that carry no SOW at all.
    is_dnav_SOW_from_this_subframe(state) || return state
    if is_decoding_completed_for_positioning(state.raw_data)
        state = dnav_confirm_data(state)
    end
    return state
end

"""
$(TYPEDSIGNATURES)

Reset a BeiDou B1I/B3I decoder state after a signal loss or reacquisition.

Clears the soft-symbol buffer, the seconds-of-week field (and its
symbol-counter anchor), and the staged D2 pages, while preserving the other
decoded parameters in `raw_data`. This allows faster recovery after brief
signal outages without a full re-decode of all subframes.

!!! note

    Like the GPS L1 C/A reset, `WN` is intentionally not cleared as it is
    broadcast only in subframe 1 (D1) / page 1 (D2); a week rollover during
    the outage yields a briefly stale week number.
"""
function reset_decoder_state(state::GNSSDecoderState{<:BeiDouDNAVData})
    empty!(state.cache.soft_buffer)
    GNSSDecoderState(
        state;
        raw_data = BeiDouDNAVData(
            state.raw_data;
            SOW = nothing,
            num_bits_after_valid_syncro_sequence_after_last_SOW = nothing,
        ),
        data = BeiDouDNAVData(),
        cache = BeiDouDNAVCache(
            state.cache;
            d2_pages = Dictionary{Int,BeiDouD2Page}(),
            sow_is_fresh = false,
        ),
        num_bits_after_valid_syncro_sequence = nothing,
    )
end
