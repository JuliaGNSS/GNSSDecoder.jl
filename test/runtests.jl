using Test, GNSSDecoder, BitIntegers, ViterbiDecoder, GNSSSignals
using Aqua
using Dictionaries

BitIntegers.@define_integers 4000

@testset "GNSSDecoder.jl" begin
    @testset "Aqua" begin
        Aqua.test_all(GNSSDecoder)
    end

    include("bit_fiddling.jl")
    include("gpsl1.jl")
    include("galileo_e1b.jl")
end