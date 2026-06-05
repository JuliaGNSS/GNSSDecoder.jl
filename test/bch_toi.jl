using Test
using Random
using GNSSDecoder
using GNSSDecoder: BCH_TOI_CODEWORDS, BCHToiSync, sync_bch_toi,
                    pack_hard_codeword, soft_to_hard_codeword,
                    TOI_RANGE, TOI_BCH_CODEWORD_LEN, TOI_BCH_MASK52

# Independently re-derive the 400-entry CNAV-2 subframe 1 codeword table,
# matching the algorithm in PocketSDR's `sync_CNV2_frame`
# (`python/sdr_nav.py`). This is the authoritative cross-check: if our
# table matches this, it matches PocketSDR's.
function pocketsdr_reference_codeword(t::Int)
    # Bit-reverse the low 8 bits of `t` into the LFSR register.
    R = 0
    for i in 0:7
        R = (R << 1) | ((t >> i) & 1)
    end
    R = UInt32(R)
    tap = UInt32(0b10011111)
    bit9 = UInt64((t >> 8) & 1)
    cw = bit9                  # symbol 0
    for i in 0:50
        lsb = R & UInt32(1)
        # PocketSDR maps LFSR LSB via CHIP=(1,-1) then (code+1)//2 -> complement.
        bit = UInt64(1) ⊻ UInt64(lsb)
        bit ⊻= bit9
        cw |= bit << (i + 1)
        feedback_bits = R & tap
        feedback = UInt32(count_ones(feedback_bits) & 1)
        R = (feedback << 7) | (R >> 1)
    end
    return cw
end

@testset "BCH(51,8) TOI codec" begin
    @testset "Codeword table" begin
        @test length(BCH_TOI_CODEWORDS) == TOI_RANGE
        for t in 0:(TOI_RANGE - 1)
            @test BCH_TOI_CODEWORDS[t + 1] == pocketsdr_reference_codeword(t)
        end
    end

    @testset "Frozen golden codewords" begin
        # The `pocketsdr_reference_codeword` cross-check above re-derives the
        # codewords with the *same* LFSR algorithm as the implementation, so
        # it only catches an inconsistency between the two — a refactor that
        # changed both in lock-step would slip through. These literals are an
        # external anchor: they were computed by an independent pure-Python
        # reimplementation of PocketSDR's `sdr_code.LFSR`/`rev_reg`
        # (`sync_CNV2_frame`, IS-GPS-800G §3.2.3.2) and frozen here. They are
        # packed first-symbol-at-bit-0, matching `BCH_TOI_CODEWORDS`.
        #
        # Coverage: both TOI-MSB polarities (t<256 ⇒ bit9=0, t≥256 ⇒ bit9=1),
        # the polarity boundary (255/256), the complement pair (0 and 256 are
        # bitwise complements within 52 bits — the documented inherent
        # ambiguity), and the modulo-400 wrap endpoint (399).
        golden = (
            0   => UInt64(0x000FFFFFFFFFFFFE),
            1   => UInt64(0x000382D25F6A30FE),
            42  => UInt64(0x000DCB2A94C5C756),
            255 => UInt64(0x00085636C0E68A00),
            256 => UInt64(0x0000000000000001),
            399 => UInt64(0x000975FEF53927E3),
        )
        for (t, cw) in golden
            @test BCH_TOI_CODEWORDS[t + 1] == cw
        end
        # The complement pair from the docstring's ambiguity note: codeword
        # for t=256 is the 52-bit bitwise complement of the codeword for t=0.
        @test BCH_TOI_CODEWORDS[0 + 1] ⊻ TOI_BCH_MASK52 == BCH_TOI_CODEWORDS[256 + 1]
    end

    @testset "Codewords are 52 bits" begin
        for cw in BCH_TOI_CODEWORDS
            @test (cw & ~TOI_BCH_MASK52) == 0
        end
    end

    @testset "Codewords are distinct" begin
        @test length(unique(BCH_TOI_CODEWORDS)) == TOI_RANGE
    end

    @testset "MSB symbol equals (TOI >> 8) & 1" begin
        # First transmitted symbol holds the MSB of the 9-bit TOI counter.
        for t in 0:(TOI_RANGE - 1)
            @test (BCH_TOI_CODEWORDS[t + 1] & 0x1) == ((t >> 8) & 0x1)
        end
    end

    @testset "pack_hard_codeword + soft_to_hard_codeword" begin
        t = 17
        cw = BCH_TOI_CODEWORDS[t + 1]
        bits = [Bool((cw >> i) & 1) for i in 0:51]
        @test pack_hard_codeword(bits) == cw
        soft = Float32[ b ? -1.0f0 : +1.0f0 for b in bits ]
        @test soft_to_hard_codeword(soft) == cw
    end

    @testset "sync detector finds matching TOI pair (unambiguous range)" begin
        # TOI values in 144..255 have no inverted twin in the codebook
        # (their would-be twin at toi+256 falls outside 0..399). They
        # yield a single unambiguous interpretation.
        for t in (144, 200, 250)
            a = BCH_TOI_CODEWORDS[t + 1]
            b = BCH_TOI_CODEWORDS[((t + 1) % TOI_RANGE) + 1]
            hit = sync_bch_toi(a, b)
            @test hit isa BCHToiSync
            @test hit.toi == t
            @test hit.polarity_flipped == false
        end
    end

    @testset "sync wraps modulo 400" begin
        # toi = 399 is in 144..255-equivalent unambiguous territory:
        # its twin would be at 655 (out of range). next TOI is 0.
        t = TOI_RANGE - 1
        a = BCH_TOI_CODEWORDS[t + 1]
        b = BCH_TOI_CODEWORDS[1]       # toi = 0
        hit = sync_bch_toi(a, b)
        @test hit isa BCHToiSync
        @test hit.toi == t
        @test hit.polarity_flipped == false
    end

    @testset "polarity-flip detection (unambiguous range)" begin
        # In the 144..255 unambiguous range, an inverted pair has no
        # normal-polarity interpretation — the detector must report a
        # flip.
        for t in (144, 200, 250)
            a = BCH_TOI_CODEWORDS[t + 1] ⊻ TOI_BCH_MASK52
            b = BCH_TOI_CODEWORDS[((t + 1) % TOI_RANGE) + 1] ⊻ TOI_BCH_MASK52
            hit = sync_bch_toi(a, b)
            @test hit isa BCHToiSync
            @test hit.toi == t
            @test hit.polarity_flipped == true
        end
    end

    @testset "Documented BCH polarity ambiguity in 0..143" begin
        # CW[t+256] == ~CW[t] for t ∈ 0..143 (the BCH construction XORs
        # the 51 LFSR bits with the MSB of the 9-bit TOI). The receiver
        # cannot tell apart "TOI=t, no flip" from "TOI=t+256, with flip"
        # without an out-of-band tie-breaker. The detector follows
        # PocketSDR's policy of returning the lowest-TOI match first.
        # Range t ∈ 0..142 keeps both `t + 256` and `t + 256 + 1` inside
        # the 0..399 codebook; at t = 143 only the t+256 codeword is in
        # range, but its inverted twin pair would need a t+257 codeword
        # in the table — t+257 = 400 is *not* in the table, so t = 143
        # is *unambiguous* even though CW[399] = ~CW[143].
        for t in (0, 17, 100, 142)
            cw_t  = BCH_TOI_CODEWORDS[t + 1]
            cw_t1 = BCH_TOI_CODEWORDS[t + 2]
            @test BCH_TOI_CODEWORDS[t + 256 + 1] == (cw_t  ⊻ TOI_BCH_MASK52)
            @test BCH_TOI_CODEWORDS[t + 256 + 2] == (cw_t1 ⊻ TOI_BCH_MASK52)
            # Detector picks the lower-TOI interpretation when both fit.
            hit = sync_bch_toi(cw_t, cw_t1)
            @test hit isa BCHToiSync
            @test hit.toi == t
            @test hit.polarity_flipped == false
        end
    end

    @testset "sync accepts iterable inputs" begin
        # t=200: unambiguous interpretation.
        t = 200
        a_bits = [Bool((BCH_TOI_CODEWORDS[t + 1] >> i) & 1) for i in 0:51]
        b_bits = [Bool((BCH_TOI_CODEWORDS[t + 2] >> i) & 1) for i in 0:51]
        hit = sync_bch_toi(a_bits, b_bits)
        @test hit isa BCHToiSync
        @test hit.toi == t
    end

    @testset "sync rejects mismatched second window" begin
        # A correct subframe N codeword followed by a NON-consecutive
        # subframe N+5 codeword must not sync.
        t = 200
        a = BCH_TOI_CODEWORDS[t + 1]
        b = BCH_TOI_CODEWORDS[t + 6]
        @test sync_bch_toi(a, b) === nothing
    end

    @testset "Monte-Carlo: no false positive on random 52-bit pairs" begin
        rng = Random.MersenneTwister(0xC0FFEE)
        trials = 1000
        false_positives = 0
        for _ in 1:trials
            a = rand(rng, UInt64) & TOI_BCH_MASK52
            b = rand(rng, UInt64) & TOI_BCH_MASK52
            if sync_bch_toi(a, b) !== nothing
                false_positives += 1
            end
        end
        @test false_positives == 0
    end
end
