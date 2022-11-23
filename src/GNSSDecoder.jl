module GNSSDecoder
    
    using DocStringExtensions, GNSSSignals, BitIntegers, ViterbiDecoder, CRC

    export decode,
        GPSL1DecoderState,
        GalileoE1BDecoderState,
        is_sat_healthy,
        GNSSDecoderState

    galCRC24 = crc(spec(24, 0x864cfb, 0x000000, false, false, 0x000000, 0xcde703))

    include("gnss.jl")
    include("bit_fiddling.jl")
    include("gpsl1.jl")
    include("galileo_e1b.jl")
end