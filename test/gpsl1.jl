BitIntegers.@define_integers 1536
BitIntegers.@define_integers 320
const GPSL1DATA = UInt1536(0x8b010c06ef056f410d000004def4351756ed43ed2357f4afe163920d2f0bff00295b9ebf9f88b010c06ef035644808d518abf98cb9d2094bfe1c3e92dbb769706d42853f19a83172d0e0748b010c06ef014ecffa302d1c9d38430077d85c9473c27eef4e6ac5a0325b0050c0cedc3c47c8b010c06eeff39075aaaac5555556aaaaaaa55555556aaaaaaa55555556aaaaaaa55555559c8b010c06eefd21840aaaab2aaaaabcaaaaaaf2aaaaabcaaaaaaf2aaaaabcaaaaaaf2aaaaabc8b)

@testset "GPS L1 decoding" begin
    decoder = GPSL1DecoderState(1)

    @test decoder.prn == 1
    @test decoder.num_bits_buffered == 0
    @test isnothing(decoder.num_bits_after_valid_subframe)
    @test decoder.raw_buffer == UInt(0)
    @test decoder.buffer == UInt(0)
    @test GNSSDecoder.calc_preamble_mask(decoder) == 0b11111111

    state = GNSSDecoder.push_bit(decoder, UInt(1))
    @test state.num_bits_buffered == 1
    @test state.raw_buffer == UInt(1)
    @test state.buffer == UInt(0)
    @test isnothing(state.num_bits_after_valid_subframe)
    @test GNSSDecoder.is_enough_buffered_bits_to_decode(state) == false
    @test GNSSDecoder.find_preamble(state) == false

    constants = GNSSDecoder.GPSL1Constants()
    @test constants.preamble == 0b10001011
    @test constants.preamble_length == 8
    @test constants.word_length == 30
    @test constants.subframe_length == 300

    raw_buffer = UInt320(constants.preamble) << UInt(300) + UInt320(constants.preamble)
    state = GNSSDecoder.GNSSDecoderState(1, raw_buffer, UInt320(0), GNSSDecoder.GPSL1Data(), GNSSDecoder.GPSL1Data(), GNSSDecoder.GPSL1Constants(), 308, nothing)
    @test GNSSDecoder.find_preamble(state) == true
    @test GNSSDecoder.complement_buffer_if_necessary(state) == GNSSDecoder.GNSSDecoderState(state, buffer = raw_buffer)
    @test GNSSDecoder.is_enough_buffered_bits_to_decode(state) == true

    raw_buffer = UInt320(~constants.preamble) << UInt(300) + UInt320(~constants.preamble)
    state = GNSSDecoder.GNSSDecoderState(1, raw_buffer, UInt320(0), GNSSDecoder.GPSL1Data(), GNSSDecoder.GPSL1Data(), GNSSDecoder.GPSL1Constants(), 308, nothing)
    @test GNSSDecoder.find_preamble(state) == true
    @test GNSSDecoder.complement_buffer_if_necessary(state) == GNSSDecoder.GNSSDecoderState(state, buffer = ~raw_buffer)

    buffer = UInt320(constants.preamble) << UInt(300) + UInt320(constants.preamble) + UInt320(1) << UInt(8)
    state = GNSSDecoder.GNSSDecoderState(1, buffer, buffer, GNSSDecoder.GPSL1Data(), GNSSDecoder.GPSL1Data(), GNSSDecoder.GPSL1Constants(), 308, nothing)
    @test GNSSDecoder.get_word(state, 10) == 1
end

@testset "GPS L1 test data decoding" begin
    decoder = GPSL1DecoderState(1)

    test_data = GNSSDecoder.GPSL1Data(
        integrity_status_flag = false,
        TOW = 34945,
        alert_flag = false,
        anti_spoof_flag = true,
        trans_week = 67,
        codeonl2 = 1,
        ura =  2.0,
        svhealth = "000000",
        IODC = "0001001000",
        l2pcode = false,
        T_GD = -1.0710209608078003e-8,
        t_oc = 216000,
        a_f2 = 0.0,
        a_f1 = -4.774847184307873e-12,
        a_f0 = -0.00018549291417002678,
        IODE_Sub_2 = "01001000",
        C_rs = 70.65625,
        Δn = 3.930878022562108e-9,
        M_0 = 2.4393048719362045,
        C_uc = 3.604218363761902e-6,
        e = 0.01144192845094949,
        C_us = 1.3023614883422852e-5,
        sqrt_A = 5153.7995529174805,
        t_oe = 216000,
        fit_interval = false,
        AODO = 31,
        C_ic = -1.73225998878479e-7,
        Ω_0 = 0.0600607702978756,
        C_is = -2.2351741790771484e-7,
        i_0 = 0.9781895349147778,
        C_rc = 136.34375,
        ω = 0.635978551768012,
        Ω_dot = -7.383521839035659e-9,
        IODE_Sub_3 = "01001000",
        IDOT = -3.4465721349922174e-10,
    )

    state = decode(decoder, GPSL1DATA, 1508)
    @test state.data == test_data

    state = decode(decoder, ~GPSL1DATA, 1508)
    @test state.data == test_data
end