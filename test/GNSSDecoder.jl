@testset "Initialization" begin
    decode = GNSSDecoder.init_decode()
    decode(0xf3d7b3701b55108d,64)
    decode(0xd3c1fd964174aeec,64)
    decode(0x2a5aff6f0939fc0b,64)
    decode(0x1bfe95ef202a27c9,64)
    decode(0xc1da1ca476336869,64)
    decode(0x806d20013add3c1f,64)

@test decode.parameters1.prev_30 == false

    decode(0xd96416dafc1fbbfd,64)
    decode(0xfaa8e8763b98a420,64)
    decode(0xd91ab37d9d2bcc58,64)
    decode(0x48a3977fc00962f5,64)
    decode(0x2afaadd3c1fd9641,64)
    decode(0x6575e13f0d31cd07,64)
    decode(0xd87c6e48fb5fc47f,64)
    decode(0xbe76803dfe0d1ea1,64)
    decode(0xd0353409e75c67c7,64)
    decode(0xdd3c1fd96415d3ac,64)
    decode(0x003ac969788a8b00,64)
    decode(0x10c9d69b5c3ec884,64)
    decode(0x407160b52c6bc016,64)
    decode(0xb0b47b06c78dd3c1,64)
    decode(0xfd964154ebe1ba3c,64)
    decode(0xc42689616af27476,64)
    decode(0xc1a020e395295dca,64)
    decode(0x9843056807b98a3e,64)
    decode(0xe11312dd3c1fd964,64)

@test decode.found_preambles.found_inverted_preamble == true
@test decode.found_preambles.preamble_pos == 62
#@show decode.data.t_oe
GNSSDecoder.calcSatPosition(478800,decode.data,decode.satPosition_param)
@show decode.satPosition_param
@show decode.found_preambles.preamble_pos
end
