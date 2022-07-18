module GNSSDecoder
    
    using DocStringExtensions, GNSSSignals, Parameters, BitIntegers

    export decode,
        GPSL1DecoderState

    include("gnss.jl")
    include("bit_fiddling.jl")
    include("gpsl1.jl")
end# end of module