# BCH codeword generation from a Fibonacci LFSR, and the GPS L1C-D
# Time-Of-Interval (TOI) table built with it.
#
# Two signals synchronise on a BCH codeword rather than a preamble, and both
# generate it the same way — an LFSR of `width` stages seeded with the
# bit-reversed information value and clocked `num_symbols` times, emitting the
# register LSB. GPS L1C-D's TOI is BCH(51,8) (this file, below); BeiDou B1C's
# subframe 1 is BCH(21,6) over the PRN concatenated with BCH(51,8) over the SOH
# (`beidou/b1c.jl`). `bch_lfsr_codeword` is that construction, taking the three
# numbers that differ; only the L1C-D table is precomputed here, because only
# its 400 entries are searched exhaustively on every symbol.
#
# Per IS-GPS-800G §3.2.3.2 the 52-symbol subframe 1 is built from a 9-bit
# TOI value `t ∈ 0..399` as follows:
#
#   - bit 0 (transmitted first) is the MSB of TOI, i.e. `(t >> 8) & 1`.
#   - bits 1..51 are the 51 output symbols of an 8-stage LFSR initialised
#     with the bit-reversed low 8 bits of `t`, clocked 51 times. The LFSR
#     feedback polynomial is `1 + x^3 + x^5 + x^6 + x^7 + x^8` (taps
#     0b10011111 on the binary state register). The resulting `{0,1}`
#     symbols are XORed with the MSB so that a polarity flip of the entire
#     subframe-1 segment can be expressed by inverting *both* the MSB and
#     the BCH output bits simultaneously.
#
# Each codeword fits in 52 bits — we pack them MSB-first into `UInt64` for
# fast equality + Hamming-distance comparison against hard-sliced soft
# symbols. The 400-entry table is `const` so it is materialised once at
# precompile time and never allocates during sync search.
#
# The full polarity flip seen at the receiver (Costas-ambiguous lock)
# corresponds to inverting every transmitted bit, which is the bitwise NOT
# of the codeword masked to 52 bits.

const TOI_BCH_REGISTER_WIDTH = 8
const TOI_BCH_TAP = 0b10011111  # IS-GPS-800G §3.2.3.2, octal 0o237
const TOI_BCH_CODEWORD_LEN = 52
const TOI_BCH_DATA_LEN = 51
const TOI_BCH_MASK52 = (UInt64(1) << TOI_BCH_CODEWORD_LEN) - UInt64(1)
const TOI_RANGE = 400  # 9-bit modulo-400 counter

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

Emit the `num_symbols` output symbols of a `width`-stage Fibonacci LFSR seeded
with the bit-reversed `value`, packed into a `UInt64` with **bit 0 = the first
emitted symbol**. `tap` is the feedback mask: the non-leading coefficients of the
generator polynomial, on the binary state register.

Seeding with the *reversed* value and emitting the register LSB is what makes the
code systematic — the first `width` emitted symbols spell `value` MSB-first —
which is the convention both ICDs that use this construction rely on, and the one
PocketSDR implements as `LFSR(n, rev_reg(v, k), taps, k)`. (Its `LFSR` emits
`CHIP[R & 1]` with `CHIP = (-1, 1)`, so its `(code + 1) / 2` recovers exactly the
register LSB this returns.)

Callers: the GPS L1C-D TOI table below (`num_symbols = 51`, `width = 8`) and
BeiDou B1C's subframe-1 PRN and SOH codewords (21/6 and 51/8, `beidou/b1c.jl`).
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

# Compute the full 52-bit codeword for a given TOI `t ∈ 0..399`.
# Bit 0 (LSB of the returned `UInt64`) is the first transmitted symbol.
#
# The 51 BCH symbols are the shared construction above with L1C-D's three
# numbers. Verified against a Spirent GSS7000 L1C recording: the first 8 emitted
# bits for TOI `t` spell the low 8 TOI bits MSB-first, i.e. the code really is
# systematic in the direction assumed here.
function _toi_codeword(t::Int)
    0 <= t < TOI_RANGE || error("TOI out of range: $t")
    msb = UInt64((t >> 8) & 1)
    bch = bch_lfsr_codeword(t & 0xff, TOI_BCH_DATA_LEN, TOI_BCH_REGISTER_WIDTH, TOI_BCH_TAP)
    # Each BCH output bit is XORed with the MSB before transmission (matches
    # PocketSDR's `((code+1)//2) ^ bit9` and IS-GPS-800G §3.2.3.2 wiring,
    # in which the MSB rides on top of the LFSR output).
    mask51 = (UInt64(1) << TOI_BCH_DATA_LEN) - UInt64(1)
    bch_xored = bch ⊻ (msb == 0 ? UInt64(0) : mask51)
    return (bch_xored << 1) | msb
end

const BCH_TOI_CODEWORDS = ntuple(i -> _toi_codeword(i - 1), TOI_RANGE)

# ---- Sync detector ----------------------------------------------------------

"""
    BCHToiSync(toi::Int, polarity_flipped::Bool)

Result of a successful multi-subframe BCH(51,8) sync. `toi` is the TOI value
of the *first* of the two subframes that matched. `polarity_flipped == true`
means the receiver is Costas-locked 180° off and every bit must be
inverted before downstream processing.
"""
struct BCHToiSync
    toi::Int
    polarity_flipped::Bool
end

"""
    pack_hard_codeword(bits) -> UInt64

Pack 52 hard-decision symbols (any iterable of `Bool`-castable values, e.g.
`Vector{Bool}`, `Vector{UInt8}`, `BitVector`) into a `UInt64` codeword with
the first symbol at bit 0. Errors if `length(bits) != 52`.

This is the encoder-side and test-side packer. The receiver never reaches it:
`try_sync` hard-slices its two windows straight off the soft-symbol deque with
`pack_bits_lsb_first` (`src/gnss.jl`), which produces the same bit order without
materialising the 52 symbols first.
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

Run the multi-subframe BCH(51,8) match used by GPS L1C-D frame sync.
`first52` and `next52` are 52-symbol windows that, if the receiver is
synchronised, hold the BCH-encoded TOI of two *consecutive* subframes. Both
inputs are accepted as either packed `UInt64` hard codewords — which is what
the receiver passes, via `pack_bits_lsb_first` (`src/gnss.jl`) — or an iterable of 52
`Bool`-castable hard decisions, which [`pack_hard_codeword`](@ref) packs.
*Soft* symbols are not accepted: hard-slicing them is the deque reader's job,
one step earlier.

Returns a [`BCHToiSync`](@ref) for the *lowest* `toi ∈ 0..399` that makes
either:

  - `first52 == BCH_TOI_CODEWORDS[toi]` and `next52 == BCH_TOI_CODEWORDS[(toi+1) mod 400]`
    — reported as `polarity_flipped == false`, or
  - the bitwise complement of `first52` and `next52` matching the same
    pair — reported as `polarity_flipped == true` (Costas-lock 180° off).

If neither holds for any TOI, returns `nothing`.

Note on inherent ambiguity: because the BCH(51,8) construction XORs the
51 LFSR bits with the MSB of the 9-bit TOI, the codeword for `t + 256`
is the bitwise complement of the codeword for `t` (whenever both are in
range). The receiver therefore cannot tell apart "TOI=`t`, no flip" from
"TOI=`t + 256`, with flip" for `t ∈ 0..143`. This function follows
PocketSDR's `sync_CNV2_frame` policy and returns the lowest-TOI match
first; downstream code uses an out-of-band check (e.g. the SF2 WN MSB)
to break the tie.

This mirrors PocketSDR's `sync_CNV2_frame` algorithm — see
`/home/schoenbrod/Code/PocketSDR/python/sdr_nav.py`.
"""
function sync_bch_toi(first52, next52)
    first_word = _as_hard_codeword(first52)
    next_word = _as_hard_codeword(next52)
    first_word_inverted = first_word ⊻ TOI_BCH_MASK52
    next_word_inverted = next_word ⊻ TOI_BCH_MASK52
    @inbounds for toi = 0:(TOI_RANGE-1)
        codeword_toi = BCH_TOI_CODEWORDS[toi+1]
        codeword_next_toi = BCH_TOI_CODEWORDS[((toi+1)%TOI_RANGE)+1]
        if first_word == codeword_toi && next_word == codeword_next_toi
            return BCHToiSync(toi, false)
        elseif first_word_inverted == codeword_toi &&
               next_word_inverted == codeword_next_toi
            return BCHToiSync(toi, true)
        end
    end
    return nothing
end

_as_hard_codeword(x::UInt64) = x & TOI_BCH_MASK52
_as_hard_codeword(x) = pack_hard_codeword(x)
