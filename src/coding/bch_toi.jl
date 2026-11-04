# BCH codeword generation from a Fibonacci LFSR, and the GPS L1C-D
# Time-Of-Interval (TOI) table built with it.
#
# Two signals synchronise on a BCH codeword rather than a preamble, generated the same
# way: an LFSR of `width` stages seeded with the bit-reversed information value and
# clocked `num_symbols` times, emitting the register LSB. GPS L1C-D's TOI is BCH(51,8);
# BeiDou B1C's subframe 1 is BCH(21,6) over the PRN plus BCH(51,8) over the SOH
# (`beidou/b1c.jl`). `bch_lfsr_codeword` is that construction; only the L1C-D table is
# precomputed, because its 400 entries are searched exhaustively on every symbol.
#
# Per IS-GPS-800G §3.2.3.2 the 52-symbol subframe 1 is built from a 9-bit TOI `t ∈ 0..399`:
#
#   - bit 0 (transmitted first) is the MSB of TOI, `(t >> 8) & 1`.
#   - bits 1..51 are the 51 output symbols of an 8-stage LFSR initialised with the
#     bit-reversed low 8 bits of `t`, clocked 51 times. Feedback polynomial
#     `1 + x^3 + x^5 + x^6 + x^7 + x^8` (taps 0b10011111 on the binary state register).
#     The symbols are XORed with the MSB, so a polarity flip of the whole segment inverts
#     both the MSB and the BCH output bits.
#
# Codewords are packed into `UInt64` (52 bits) for fast equality comparison against
# hard-sliced soft symbols. The table is `const`, so sync search never allocates.
#
# The polarity flip of a Costas-ambiguous lock inverts every bit: the bitwise NOT of the
# codeword masked to 52 bits.

const TOI_BCH_REGISTER_WIDTH = 8
const TOI_BCH_TAP = 0b10011111  # IS-GPS-800G §3.2.3.2, octal 0o237
const TOI_BCH_CODEWORD_LEN = 52
const TOI_BCH_DATA_LEN = 51
const TOI_BCH_MASK52 = (UInt64(1) << TOI_BCH_CODEWORD_LEN) - UInt64(1)
const TOI_RANGE = 400  # 9-bit modulo-400 counter
# The TOI's MSB (weight 256) is XORed over the whole codeword, so
# `codeword(t + 256) == ~codeword(t)`: a hard-decision match cannot tell (TOI, polarity)
# from (TOI + 256, flipped polarity) when both are in range. The consumer resolves the
# pair; see `decode_sync_unit` in `src/gps/l1c_d.jl`.
const TOI_COMPLEMENT_OFFSET = 256

# Reverse the low `n` bits of `x` (matches PocketSDR's `sdr_code.rev_reg`).
function reverse_low_bits(x::Integer, n::Int)
    r = UInt(0)
    @inbounds for i = 0:(n-1)
        r = (r << 1) | ((UInt(x) >> i) & UInt(1))
    end
    return r
end

"""
    bch_lfsr_codeword(value, num_symbols, width, tap) -> UInt64

Emit the `num_symbols` output symbols of a `width`-stage Fibonacci LFSR seeded with the
bit-reversed `value`, packed into a `UInt64` with **bit 0 = the first emitted symbol**.
`tap` is the feedback mask: the non-leading generator-polynomial coefficients on the
binary state register.

Seeding with the reversed value and emitting the register LSB makes the code systematic
(the first `width` symbols spell `value` MSB-first), as both ICDs assume.
PocketSDR implements it as `LFSR(n, rev_reg(v, k), taps, k)`; its `(code + 1) / 2` is the
register LSB.

Callers: the GPS L1C-D TOI table (`num_symbols = 51`, `width = 8`) and BeiDou B1C's
subframe-1 PRN and SOH codewords (21/6 and 51/8, `beidou/b1c.jl`).
"""
function bch_lfsr_codeword(value::Integer, num_symbols::Int, width::Int, tap::Integer)
    R = reverse_low_bits(UInt(value) & ((UInt(1) << width) - UInt(1)), width) % UInt32
    out = UInt64(0)
    @inbounds for i = 0:(num_symbols-1)
        out |= UInt64(R & UInt32(0x1)) << i
        # Galois-style update: new MSB = parity of R AND tap; shift right.
        feedback = UInt32(count_ones(R & UInt32(tap)) & 1)
        R = (feedback << (width - 1)) | (R >> 1)
    end
    return out
end

# Full 52-bit codeword for TOI `t ∈ 0..399`; bit 0 is the first transmitted symbol.
# Verified against a Spirent GSS7000 L1C recording: the first 8 emitted bits spell the low
# 8 TOI bits MSB-first (systematic, as assumed).
function _toi_codeword(t::Int)
    0 <= t < TOI_RANGE || error("TOI out of range: $t")
    msb = UInt64((t >> 8) & 1)
    bch = bch_lfsr_codeword(t & 0xff, TOI_BCH_DATA_LEN, TOI_BCH_REGISTER_WIDTH, TOI_BCH_TAP)
    # Each BCH output bit is XORed with the MSB (PocketSDR's `((code+1)//2) ^ bit9`,
    # IS-GPS-800G §3.2.3.2).
    mask51 = (UInt64(1) << TOI_BCH_DATA_LEN) - UInt64(1)
    bch_xored = bch ⊻ (msb == 0 ? UInt64(0) : mask51)
    return (bch_xored << 1) | msb
end

const BCH_TOI_CODEWORDS = ntuple(i -> _toi_codeword(i - 1), TOI_RANGE)

# ---- Sync detector ----------------------------------------------------------

"""
    BCHToiSync(TOI::Int, polarity_flipped::Bool)

Result of a successful two-subframe BCH(51,8) sync. `TOI` is that of the *first* matched
subframe. `polarity_flipped == true` means a Costas lock 180° off: invert every bit
downstream.
"""
struct BCHToiSync
    TOI::Int
    polarity_flipped::Bool
end

"""
    pack_hard_codeword(bits) -> UInt64

Pack 52 hard-decision symbols (any iterable of `Bool`-castable values) into a `UInt64`
codeword with the first symbol at bit 0. Errors if `length(bits) != 52`.

Encoder- and test-side only: `try_sync` hard-slices straight off the soft-symbol deque
with `pack_bits_lsb_first` (`src/gnss.jl`), same bit order.
"""
function pack_hard_codeword(bits)
    length(bits) == TOI_BCH_CODEWORD_LEN ||
        error("expected $(TOI_BCH_CODEWORD_LEN) symbols, got $(length(bits))")
    word = UInt64(0)
    @inbounds for (i, b) in enumerate(bits)
        if Bool(b)
            word |= UInt64(1) << (i - 1)
        end
    end
    return word
end

"""
    sync_bch_toi(first52, next52) -> Union{BCHToiSync, Nothing}

Two-subframe BCH(51,8) match used by GPS L1C-D frame sync. `first52` and `next52` are
52-symbol windows holding the BCH-encoded TOI of two *consecutive* subframes. Each is a
packed `UInt64` hard codeword (what the receiver passes, via `pack_bits_lsb_first` in
`src/gnss.jl`) or an iterable of 52 `Bool`-castable hard decisions, packed by
[`pack_hard_codeword`](@ref). Soft symbols are not accepted.

Returns a [`BCHToiSync`](@ref) for the *lowest* `TOI ∈ 0..399` that makes either:

  - `first52 == BCH_TOI_CODEWORDS[TOI]` and `next52 == BCH_TOI_CODEWORDS[(TOI+1) mod 400]`
    — reported as `polarity_flipped == false`, or
  - the bitwise complements of `first52` and `next52` matching the same pair — reported
    as `polarity_flipped == true` (Costas lock 180° off).

Otherwise returns `nothing`.

Inherent ambiguity: since the MSB of the TOI is XORed over the 51 LFSR bits, the codeword
for `t + 256` is the complement of that for `t` (when both are in range). "TOI=`t`, no flip"
and "TOI=`t + 256`, flipped" are indistinguishable for `t ∈ 0..143`; this function always
reports the *low* branch (PocketSDR's `sync_CNV2_frame` order). The caller must resolve the
pair: `decode_sync_unit` in `src/gps/l1c_d.jl` uses TOI continuity with the previous frame,
else tries subframe 2's LDPC+CRC under both branches (the CRC is the only in-band oracle).
Taking the reported branch at face value is wrong for every frame with true TOI ≥ 256, i.e.
43 minutes of every 2-hour interval.
"""
function sync_bch_toi(first52, next52)
    first_word = _as_hard_codeword(first52)
    next_word = _as_hard_codeword(next52)
    first_word_inverted = first_word ⊻ TOI_BCH_MASK52
    next_word_inverted = next_word ⊻ TOI_BCH_MASK52
    @inbounds for TOI = 0:(TOI_RANGE-1)
        codeword_toi = BCH_TOI_CODEWORDS[TOI+1]
        codeword_next_toi = BCH_TOI_CODEWORDS[((TOI+1)%TOI_RANGE)+1]
        if first_word == codeword_toi && next_word == codeword_next_toi
            return BCHToiSync(TOI, false)
        elseif first_word_inverted == codeword_toi &&
               next_word_inverted == codeword_next_toi
            return BCHToiSync(TOI, true)
        end
    end
    return nothing
end

_as_hard_codeword(x::UInt64) = x & TOI_BCH_MASK52
_as_hard_codeword(x) = pack_hard_codeword(x)
