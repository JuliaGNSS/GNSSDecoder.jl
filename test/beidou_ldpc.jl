using Test
using Random
import Aff3ct
using Aff3ct: LDPCMatrix, LDPCBPDecoder
using GNSSDecoder

# The BeiDou B-CNAV 64-ary LDPC pipeline: ICD-transcribed GF(2^6) matrices
# (`scripts/beidou_ldpc_coords.jl`) → binary-image `.alist` artefacts in
# `data/` → Aff3ct binary BP decoding through the shared
# `GNSSDecoder.load_ldpc_decoder` / `GNSSDecoder.ldpc_decode_word` helpers.
@testset "BeiDou LDPC alist files" begin
    repo_root = normpath(joinpath(@__DIR__, ".."))
    # Also brings the GF(2^6) reference arithmetic + systematic encoder into
    # scope (`gf64_mul`, `beidou_ldpc_encode`, `gf64_symbols_to_bits`, …); the
    # per-signal BeiDou tests reuse them to synthesize transmit chains.
    isdefined(Main, :beidou_ldpc_encode) ||
        include(joinpath(repo_root, "scripts", "generate_beidou_alist.jl"))

    committed(name) = joinpath(repo_root, "data", "$(name).alist")

    @testset "GF(2^6) arithmetic (p(x) = 1 + x + x^6)" begin
        # α = x = 2: α^6 = x + 1 = 3 under p(x) = x^6 + x + 1.
        α6 = foldl((a, _) -> gf64_mul(a, 2), 1:6; init = 1)
        @test α6 == 3
        # The multiplicative group has order 63: α^63 = 1.
        @test foldl((a, _) -> gf64_mul(a, 2), 1:63; init = 1) == 1
        # Every unit has an inverse consistent with multiplication.
        @test all(gf64_mul(v, gf64_inv(v)) == 1 for v = 1:63)
        # ICD bit mapping: symbol 1 ⇔ 000001, MSB first.
        @test gf64_symbols_to_bits([1]) == [false, false, false, false, false, true]
        @test gf64_bits_to_symbols(gf64_symbols_to_bits([35, 13])) == [35, 13]
    end

    @testset "Committed files exist and have expected dimensions" begin
        for (name, n, k) in (
            ("bcnv1_sf2", 200, 100),
            ("bcnv1_sf3", 88, 44),
            ("bcnv2", 96, 48),
            ("bcnv3", 162, 81),
        )
            path = committed(name)
            @test isfile(path)
            first_line = open(io -> readline(io), path)
            N, M = parse.(Int, split(first_line))
            @test N == 6n
            @test M == 6 * (n - k)
        end
    end

    # The generator is fully self-contained — the ICD matrix coordinates live
    # in `scripts/beidou_ldpc_coords.jl`, so regeneration runs everywhere (CI,
    # fresh clones) with no external dependency. This pins the committed
    # `.alist` artefacts to the spec coordinates byte-for-byte. It also
    # re-runs the internal cross-check that the binary image annihilates
    # codewords of the independent GF(2^6) systematic encoder.
    @testset "Generator is reproducible (byte-compare)" begin
        mktempdir() do tmp
            generate_beidou_alist(tmp)
            for name in ("bcnv1_sf2", "bcnv1_sf3", "bcnv2", "bcnv3")
                regen = joinpath(tmp, "$(name).alist")
                @test isfile(regen)
                @test read(regen) == read(committed(name))
            end
        end
    end

    @testset "GF(2^6) encode -> binary-image BP decode round-trip" begin
        rng = Random.MersenneTwister(0xBD5)
        for (name, code) in pairs(BEIDOU_LDPC_CODES)
            decoder = GNSSDecoder.load_ldpc_decoder(committed(String(name)))
            for _ = 1:5
                message = rand(rng, 0:63, code.k)
                bits = gf64_symbols_to_bits(beidou_ldpc_encode(code, message))
                # Noisy LLRs (positive ⇒ bit 0), well inside BP's margin.
                llr = Float32[(b ? -2.0f0 : 2.0f0) for b in bits]
                llr .+= 0.8f0 .* randn(rng, Float32, length(llr))
                info = Aff3ct.decode(decoder, llr)
                @test Vector{Bool}(info .== 1) == bits[1:(6*code.k)]
            end
        end
    end
end
