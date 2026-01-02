BitIntegers.@define_integers 1536
BitIntegers.@define_integers 320
const GPSL1DATA =
    uint1536"0x8b010c06ef056f410d000004def4351756ed43ed2357f4afe163920d2f0bff00295b9ebf9f88b010c06ef035644808d518abf98cb9d2094bfe1c3e92dbb769706d42853f19a83172d0e0748b010c06ef014ecffa302d1c9d38430077d85c9473c27eef4e6ac5a0325b0050c0cedc3c47c8b010c06eeff39075aaaac5555556aaaaaaa55555556aaaaaaa55555556aaaaaaa55555559c8b010c06eefd21840aaaab2aaaaabcaaaaaaf2aaaaabcaaaaaaf2aaaaabcaaaaaaf2aaaaabc8b"

@testset "GPS L1 constructor" begin
    gpsl1 = GPSL1()

    state_1 = GPSL1DecoderState(21)
    state_2 = GNSSDecoderState(gpsl1, 21)
    @test state_1.prn == state_2.prn &&
        state_1.raw_buffer == state_2.raw_buffer &&
        state_1.buffer == state_2.buffer &&
        state_1.raw_data == state_2.raw_data &&
        state_1.data == state_2.data &&
        state_1.constants == state_2.constants &&
        state_1.num_bits_buffered == state_2.num_bits_buffered &&
        state_1.num_bits_after_valid_syncro_sequence == state_2.num_bits_after_valid_syncro_sequence &&
        state_1.is_shifted_by_180_degrees == state_2.is_shifted_by_180_degrees
end

@testset "GPS L1 decoding" begin
    decoder = GPSL1DecoderState(1)

    @test decoder.prn == 1
    @test decoder.num_bits_buffered == 0
    @test isnothing(decoder.num_bits_after_valid_syncro_sequence)
    @test decoder.raw_buffer == UInt(0)
    @test decoder.buffer == UInt(0)
    @test GNSSDecoder.calc_preamble_mask(decoder) == 0b11111111

    state = GNSSDecoder.push_bit(decoder, UInt(1))
    @test state.num_bits_buffered == 1
    @test state.raw_buffer == UInt(1)
    @test state.buffer == UInt(0)
    @test isnothing(state.num_bits_after_valid_syncro_sequence)
    @test GNSSDecoder.is_enough_buffered_bits_to_decode(state) == false
    @test GNSSDecoder.find_preamble(state) == false

    constants = GNSSDecoder.GPSL1Constants()
    @test constants.preamble == 0b10001011
    @test constants.preamble_length == 8
    @test constants.word_length == 30
    @test constants.syncro_sequence_length == 300

    raw_buffer = UInt320(constants.preamble) << UInt(300) + UInt320(constants.preamble)
    state = GNSSDecoder.GNSSDecoderState(
        1,
        raw_buffer,
        UInt320(0),
        GNSSDecoder.GPSL1Data(),
        GNSSDecoder.GPSL1Data(),
        GNSSDecoder.GPSL1Constants(),
        GNSSDecoder.GPSL1Cache(),
        308,
        nothing,
        false,
    )
    @test GNSSDecoder.find_preamble(state) == true
    @test GNSSDecoder.complement_buffer_if_necessary(state) == GNSSDecoder.GNSSDecoderState(
        state;
        buffer = raw_buffer,
        is_shifted_by_180_degrees = false,
    )
    @test GNSSDecoder.is_enough_buffered_bits_to_decode(state) == true

    raw_buffer = UInt320(~constants.preamble) << UInt(300) + UInt320(~constants.preamble)
    state = GNSSDecoder.GNSSDecoderState(
        1,
        raw_buffer,
        UInt320(0),
        GNSSDecoder.GPSL1Data(),
        GNSSDecoder.GPSL1Data(),
        GNSSDecoder.GPSL1Constants(),
        GNSSDecoder.GPSL1Cache(),
        308,
        nothing,
        false,
    )
    @test GNSSDecoder.find_preamble(state) == true
    @test GNSSDecoder.complement_buffer_if_necessary(state) == GNSSDecoder.GNSSDecoderState(
        state;
        buffer = ~raw_buffer,
        is_shifted_by_180_degrees = true,
    )

    buffer =
        UInt320(constants.preamble) << UInt(300) +
        UInt320(constants.preamble) +
        UInt320(1) << UInt(8)
    state = GNSSDecoder.GNSSDecoderState(
        1,
        buffer,
        buffer,
        GNSSDecoder.GPSL1Data(),
        GNSSDecoder.GPSL1Data(),
        GNSSDecoder.GPSL1Constants(),
        GNSSDecoder.GPSL1Cache(),
        308,
        nothing,
        false,
    )
    @test GNSSDecoder.get_word(state, 10) == 1
end

@testset "GPS L1 test data decoding" begin
    decoder = GPSL1DecoderState(1)

    test_data = GNSSDecoder.GPSL1Data(;
        last_subframe_id = 5,
        integrity_status_flag = false,
        TOW = 34945 * 6,
        alert_flag = false,
        anti_spoof_flag = true,
        trans_week = 67,
        codeonl2 = 1,
        ura = 2.0,
        svhealth = "000000",
        IODC = "0001001000",
        l2pcode = false,
        T_GD = -1.0710209608078003e-8,
        t_0c = 216000,
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
        t_0e = 216000,
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
        i_dot = -3.4465721349922174e-10,
    )

    state = decode(decoder, GPSL1DATA, 1508)
    @test state.data == test_data
    @test is_sat_healthy(state) == true

    state = decode(decoder, ~GPSL1DATA, 1508)
    @test state.data == test_data
    @test is_sat_healthy(state) == true

    # test confirm_data
    state = GNSSDecoder.GNSSDecoderState(
        state;
        raw_data = GNSSDecoder.GPSL1Data(state.data, C_ic = state.data.C_ic+1)
    )
    state = GNSSDecoder.confirm_data(state)
    @test state.data.C_ic == test_data.C_ic # erroneous data not accepted
    state = GNSSDecoder.GNSSDecoderState(
        state;
        raw_data = GNSSDecoder.GPSL1Data(state.data, C_ic = state.data.C_ic+1)
    )
    state = GNSSDecoder.confirm_data(state)
    @test state.data.C_ic == test_data.C_ic # erroneous data not accepted
    state = GNSSDecoder.GNSSDecoderState(
        state;
        raw_data = GNSSDecoder.GPSL1Data(state.data, C_ic = state.data.C_ic+1)
    )
    state = GNSSDecoder.confirm_data(state)
    @test state.data.C_ic == test_data.C_ic+1 # erroneous data accepted as it has been provided more often than true data
end

@testset "confirm_data branches" begin
    # Setup: create a base state with valid data
    decoder = GPSL1DecoderState(1)
    state = decode(decoder, GPSL1DATA, 1508)
    base_data = state.data

    # Branch 1: New IODC (not in cache) - should use data immediately
    @testset "new IODC - uses data immediately" begin
        # Start with empty cache
        state_empty_cache = GNSSDecoder.GNSSDecoderState(
            state;
            cache = GNSSDecoder.GPSL1Cache(),
            raw_data = base_data,
            data = GNSSDecoder.GPSL1Data(),
        )
        result = GNSSDecoder.confirm_data(state_empty_cache)
        @test result.data == base_data  # data is used
        @test length(result.cache.old_data) == 1  # added to cache
        @test result.cache.old_data[1].vote == 0  # initial vote is 0
        @test result.cache.old_data[1].data == base_data
    end

    # Branch 2: Same IODC but different data - adds to cache, doesn't use data yet
    @testset "same IODC different data - adds to cache, doesn't use" begin
        # Create cache with one entry
        existing_data = GNSSDecoder.GPSL1Data(base_data; C_ic = base_data.C_ic + 1.0)
        state_with_cache = GNSSDecoder.GNSSDecoderState(
            state;
            cache = GNSSDecoder.GPSL1Cache([GNSSDecoder.VotedGPSL1Data(5, existing_data)]),
            raw_data = base_data,  # different data, same IODC
            data = GNSSDecoder.GPSL1Data(),
        )
        result = GNSSDecoder.confirm_data(state_with_cache)
        @test result.data == GNSSDecoder.GPSL1Data()  # data NOT used (empty)
        @test result.raw_data == GNSSDecoder.GPSL1Data()  # raw_data reset
        @test length(result.cache.old_data) == 2  # new entry added
        @test result.cache.old_data[2].vote == 0  # new entry has vote 0
        @test result.cache.old_data[2].data == base_data
    end

    # Branch 3: Matching entry exists but another entry has higher score - rejects data
    @testset "matching entry with lower score than best - rejects data" begin
        # Create cache with two entries, one with higher score
        high_score_data = GNSSDecoder.GPSL1Data(base_data; C_ic = base_data.C_ic + 1.0)
        low_score_data = base_data
        state_with_competing = GNSSDecoder.GNSSDecoderState(
            state;
            cache = GNSSDecoder.GPSL1Cache([
                GNSSDecoder.VotedGPSL1Data(10, high_score_data),  # higher score
                GNSSDecoder.VotedGPSL1Data(2, low_score_data),    # lower score, matches raw_data
            ]),
            raw_data = base_data,  # matches low_score_data
            data = high_score_data,  # currently using high score data
        )
        result = GNSSDecoder.confirm_data(state_with_competing)
        @test result.data == high_score_data  # data NOT changed (still high score data)
        @test result.raw_data == GNSSDecoder.GPSL1Data()  # raw_data reset
        @test result.cache.old_data[2].vote == 3  # vote incremented from 2 to 3
    end

    # Branch 4: Matching entry has best (or tied best) score - uses data
    @testset "matching entry with best score - uses data" begin
        state_with_match = GNSSDecoder.GNSSDecoderState(
            state;
            cache = GNSSDecoder.GPSL1Cache([GNSSDecoder.VotedGPSL1Data(5, base_data)]),
            raw_data = base_data,
            data = GNSSDecoder.GPSL1Data(),
        )
        result = GNSSDecoder.confirm_data(state_with_match)
        @test result.data == base_data  # data is used
        @test result.cache.old_data[1].vote == 6  # vote incremented from 5 to 6
    end

    # Branch 5: Matching entry at max_vote - removes other entries
    @testset "max_vote reached - removes other entries" begin
        other_data = GNSSDecoder.GPSL1Data(base_data; C_ic = base_data.C_ic + 1.0)
        state_at_max = GNSSDecoder.GNSSDecoderState(
            state;
            cache = GNSSDecoder.GPSL1Cache([
                GNSSDecoder.VotedGPSL1Data(20, base_data),  # already at max_vote (20)
                GNSSDecoder.VotedGPSL1Data(5, other_data),  # should be removed
            ]),
            raw_data = base_data,
            data = GNSSDecoder.GPSL1Data(),
        )
        result = GNSSDecoder.confirm_data(state_at_max)
        @test result.data == base_data  # data is used
        @test length(result.cache.old_data) == 1  # only one entry remains
        @test result.cache.old_data[1].vote == 20  # stays at max_vote (capped)
        @test result.cache.old_data[1].data == base_data
    end

    # Branch 6: max_vote reached but only one entry - doesn't try to remove
    @testset "max_vote reached with single entry - keeps entry" begin
        state_single_max = GNSSDecoder.GNSSDecoderState(
            state;
            cache = GNSSDecoder.GPSL1Cache([GNSSDecoder.VotedGPSL1Data(19, base_data)]),
            raw_data = base_data,
            data = GNSSDecoder.GPSL1Data(),
        )
        result = GNSSDecoder.confirm_data(state_single_max)
        @test result.data == base_data
        @test length(result.cache.old_data) == 1
        @test result.cache.old_data[1].vote == 20
    end

    # Branch 7: Vote capped at max_vote
    @testset "vote capped at max_vote" begin
        state_at_max = GNSSDecoder.GNSSDecoderState(
            state;
            cache = GNSSDecoder.GPSL1Cache([GNSSDecoder.VotedGPSL1Data(20, base_data)]),
            raw_data = base_data,
            data = GNSSDecoder.GPSL1Data(),
        )
        result = GNSSDecoder.confirm_data(state_at_max)
        @test result.cache.old_data[1].vote == 20  # stays at max, doesn't overflow
    end

    # Edge case: Different IODC added to non-empty cache
    @testset "different IODC added to existing cache" begin
        existing_data = GNSSDecoder.GPSL1Data(base_data; IODC = "1111111111")
        new_data = base_data  # has different IODC
        state_diff_iodc = GNSSDecoder.GNSSDecoderState(
            state;
            cache = GNSSDecoder.GPSL1Cache([GNSSDecoder.VotedGPSL1Data(10, existing_data)]),
            raw_data = new_data,
            data = GNSSDecoder.GPSL1Data(),
        )
        result = GNSSDecoder.confirm_data(state_diff_iodc)
        @test result.data == new_data  # new IODC data is used immediately
        @test length(result.cache.old_data) == 2  # both entries kept
        @test result.cache.old_data[1].data.IODC == "1111111111"  # old entry preserved
        @test result.cache.old_data[2].data == new_data  # new entry added
    end
end

@testset "GPS L1 reset_decoder_state" begin
    decoder = GPSL1DecoderState(1)
    state = decode(decoder, GPSL1DATA, 1508)

    # test reset_decoder_state
    state = reset_decoder_state(state)
    @test state.raw_buffer == 0
    @test state.buffer == 0
    @test isnothing(state.raw_data.TOW)
    @test isnothing(state.data.TOW)
    @test state.num_bits_buffered == 0
    @test isnothing(state.num_bits_after_valid_syncro_sequence)
end