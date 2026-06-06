using Test
using Random
# Import the module for the qualified `Aff3ct.decode` call, and bring only the
# LDPC names into scope. `using Aff3ct` would also export `decode` and
# `ViterbiDecoder`, clashing with `GNSSDecoder.decode` and the `ViterbiDecoder`
# package already in scope from runtests.jl.
import Aff3ct
using Aff3ct: LDPCMatrix, LDPCEncoder, LDPCBPDecoder, encode
using GNSSDecoder

# Re-run the generator helper to materialise the SF2 / SF3 `.alist` files
# in a temporary directory, then byte-compare against the committed files.
# This guarantees the committed artefacts are reproducible from the
# PocketSDR cross-reference and the generator helper that ships with the
# repo.
@testset "LDPC alist files" begin
    repo_root = normpath(joinpath(@__DIR__, ".."))
    committed_sf2 = joinpath(repo_root, "data", "cnv2_sf2.alist")
    committed_sf3 = joinpath(repo_root, "data", "cnv2_sf3.alist")

    @testset "Committed files exist and have expected dimensions" begin
        @test isfile(committed_sf2)
        @test isfile(committed_sf3)
        # First line of an alist is "N M".
        for (path, n_expected, m_expected) in (
            (committed_sf2, 1200, 600),
            (committed_sf3, 548,  274),
        )
            first_line = open(io -> readline(io), path)
            n, m = parse.(Int, split(first_line))
            @test n == n_expected
            @test m == m_expected
        end
    end

    # The generator is fully self-contained — the IS-GPS-800G coordinates live
    # in `scripts/cnv2_ldpc_coords.jl`, so regeneration runs everywhere (CI,
    # fresh clones) with no external dependency. This pins the committed
    # `.alist` artefacts to the spec coordinates byte-for-byte.
    @testset "Generator is reproducible (byte-compare)" begin
        include(joinpath(repo_root, "scripts", "generate_alist.jl"))
        mktempdir() do tmp
            generate_alist(tmp)
            for (fname, committed) in (
                ("cnv2_sf2.alist", committed_sf2),
                ("cnv2_sf3.alist", committed_sf3),
            )
                regen = joinpath(tmp, fname)
                @test isfile(regen)
                @test read(regen) == read(committed)
            end
        end
    end

    @testset "AFF3CT round-trip: SF2 encode -> BP decode" begin
        H = LDPCMatrix(committed_sf2)
        @test H.M == 600
        @test H.N == 1200
        @test H.K == 600
        encoder = LDPCEncoder(H)
        decoder = LDPCBPDecoder(H; num_iterations = 50)
        rng = Random.MersenneTwister(0x5F2)
        U_K = Int32.(rand(rng, Bool, H.K))
        X_N = encode(encoder, U_K)
        # BPSK soft mapping matches AFF3CT.jl's own test_ldpc.jl.
        Y_N = Float32[100.0f0 * (1.0f0 - 2.0f0 * x) for x in X_N]
        V_K = Aff3ct.decode(decoder, Y_N)
        @test V_K == U_K
    end

    @testset "AFF3CT round-trip: SF3 encode -> BP decode" begin
        H = LDPCMatrix(committed_sf3)
        @test H.M == 274
        @test H.N == 548
        @test H.K == 274
        encoder = LDPCEncoder(H)
        decoder = LDPCBPDecoder(H; num_iterations = 50)
        rng = Random.MersenneTwister(0x5F3)
        U_K = Int32.(rand(rng, Bool, H.K))
        X_N = encode(encoder, U_K)
        Y_N = Float32[100.0f0 * (1.0f0 - 2.0f0 * x) for x in X_N]
        V_K = Aff3ct.decode(decoder, Y_N)
        @test V_K == U_K
    end
end
