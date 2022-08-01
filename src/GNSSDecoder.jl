module GNSSDecoder
    
    using DocStringExtensions, GNSSSignals, BitIntegers, ViterbiDecoder

    export decode,
        GPSL1DecoderState,
        GalileoE1BDecoderState

    include("gnss.jl")
    include("bit_fiddling.jl")
    include("gpsl1.jl")
    include("galileo_e1b.jl")
end