using Test, GNSSDecoder, BitIntegers, ViterbiDecoder, GNSSSignals
using Aqua
using Dictionaries

@testset "GNSSDecoder.jl" begin
    @testset "Aqua" begin
        Aqua.test_all(GNSSDecoder)
    end

    include("bit_fiddling.jl")
    include("gpsl1.jl")
    include("gpsl5.jl")
    include("galileo_e1b.jl")
end