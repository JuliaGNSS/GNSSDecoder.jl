module GNSSDecoder

using DocStringExtensions,
    GNSSSignals, BitIntegers, ViterbiDecoder, CRC, Dictionaries, DataStructures

export decode,
    GPSL1CADecoderState,
    GalileoE1BDecoderState,
    is_sat_healthy,
    GNSSDecoderState,
    reset_decoder_state

const galCRC24 = crc(spec(24, 0x864cfb, 0x000000, false, false, 0x000000, 0xcde703))

include("gnss.jl")
include("bit_fiddling.jl")
include("gps/l1ca.jl")
include("galileo/e1b.jl")
end
