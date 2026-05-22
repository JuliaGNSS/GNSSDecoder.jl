using Test
using Random
using GNSSDecoder
using GNSSDecoder: deinterleave!, interleave!

@testset "Block deinterleaver" begin
    @testset "Round-trip identity (38×46, Float32)" begin
        rng = Random.MersenneTwister(0xDEC0DE)
        src = randn(rng, Float32, 38 * 46)
        tx  = similar(src)
        rx  = similar(src)
        interleave!(tx, src, 38, 46)
        deinterleave!(rx, tx, 38, 46)
        @test rx == src
    end

    @testset "Round-trip identity (other sizes, Int8)" begin
        rng = Random.MersenneTwister(0xCAFE)
        for (rows, cols) in ((4, 5), (16, 32), (38, 46), (1, 100), (100, 1))
            src = rand(rng, Int8, rows * cols)
            tx  = similar(src)
            rx  = similar(src)
            interleave!(tx, src, rows, cols)
            deinterleave!(rx, tx, rows, cols)
            @test rx == src
        end
    end

    @testset "Round-trip identity (Bool)" begin
        rng = Random.MersenneTwister(7)
        src = rand(rng, Bool, 38 * 46)
        tx  = similar(src)
        rx  = similar(src)
        interleave!(tx, src, 38, 46)
        deinterleave!(rx, tx, 38, 46)
        @test rx == src
    end

    @testset "Known fixture: 4×3 example (12 symbols)" begin
        # Hand-traced: pre-interleaver (row-major) [1..12] becomes
        # column-major-read [1,4,7,10, 2,5,8,11, 3,6,9,12] on the wire.
        src = collect(1:12)
        wire = similar(src)
        interleave!(wire, src, 4, 3)
        @test wire == [1, 4, 7, 10, 2, 5, 8, 11, 3, 6, 9, 12]
        # Deinterleaver recovers the original.
        rx = similar(src)
        deinterleave!(rx, wire, 4, 3)
        @test rx == src
    end

    @testset "Dimension mismatch errors" begin
        a = zeros(Float32, 10)
        b = zeros(Float32, 12)
        @test_throws DimensionMismatch deinterleave!(a, b, 3, 4)
        @test_throws DimensionMismatch interleave!(a, b, 3, 4)
    end
end
