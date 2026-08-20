# BeiDou B3I decoder tests. The D1/D2 message core is shared with B1I (and
# fully exercised in `beidou_b1i.jl`); these tests pin the B3I signal layer:
# constructor mapping, signal identity, and one end-to-end D1 decode through
# the shared core.

@testset "BeiDou B3I" begin
    PI = GNSSDecoder.GNSS_PI

    @testset "Signal metadata and constructors" begin
        state = BeiDouB3IDecoderState(25)
        @test get_signal_type(state) == BeiDouB3I
        @test GNSSDecoderState(BeiDouB3I(), 25) == state
        @test get_data_frequency(state) == get_data_frequency(BeiDouB3I)
        @test get_time_system_name(state) == "BeiDou Time"
        # B1I and B3I share the data container but stay distinct signals.
        @test typeof(state.constants) != typeof(BeiDouB1IDecoderState(25).constants)
    end

    @testset "D1 subframes 1-3 decode through the shared core (PRN 25)" begin
        SOW0 = 86400
        d1_sf1 = (
            SOW = SOW0,
            SatH1 = 0,
            AODC = 2,
            URAI = 0,
            WN = 900,
            t_0c_raw = 900,
            T_GD1_raw = 11,
            T_GD2_raw = -22,
            α_raw = (1, 2, 3, 4),
            β_raw = (5, 6, 7, 8),
            a_2_raw = 1,
            a_0_raw = -100,
            a_1_raw = 200,
            AODE = 3,
        )
        d1_sf2 = (
            SOW = SOW0 + 6,
            Δn_raw = -111,
            C_uc_raw = 222,
            M_0_raw = -333,
            e_raw = 444,
            C_us_raw = -555,
            C_rc_raw = 666,
            C_rs_raw = -777,
            sqrt_A_raw = 888,
            t_0e_msb2 = 0,
        )
        d1_sf3 = (
            SOW = SOW0 + 12,
            t_0e_lsb15 = 900,
            i_0_raw = -999,
            C_ic_raw = 123,
            Ω_dot_raw = -456,
            C_is_raw = 789,
            i_dot_raw = -12,
            Ω_0_raw = 3456,
            ω_raw = -6789,
        )

        # Two broadcast cycles: the shared core publishes a dataset only after
        # the broadcast repeats it (see `beidou_b1i.jl`). Consecutive
        # subframes stay 6 s apart, so a three-subframe cycle repeats at 18 s.
        cycle(Δ) = (
            dnav_test_encode_subframe(
                dnav_test_d1_subframe1_content(; merge(d1_sf1, (SOW = SOW0 + Δ,))...),
            ),
            dnav_test_encode_subframe(
                dnav_test_d1_subframe2_content(; merge(d1_sf2, (SOW = SOW0 + Δ + 6,))...),
            ),
            dnav_test_encode_subframe(
                dnav_test_d1_subframe3_content(; merge(d1_sf3, (SOW = SOW0 + Δ + 12,))...),
            ),
        )
        state = BeiDouB3IDecoderState(25)
        symbols = dnav_test_soft_symbols(cycle(0)..., cycle(18)...)
        state = decode(state, symbols, length(symbols))
        @test is_decoding_completed_for_positioning(state)
        @test is_sat_healthy(state)
        @test state.data.WN == 900
        @test state.data.t_0c == 7200
        @test state.data.t_0e == 7200
        @test state.data.a_0 == -100 / 2.0^33
        @test state.data.Δn == -111 * PI / 2.0^43
        @test state.data.ω == -6789 * PI / 2.0^31
        @test state.data.ura == 2.0            # URAI = 0
        # The shared core's SOW epoch handling reaches B3I too: 311 symbols
        # (one subframe + the trailing preamble) past the last SOW.
        @test state.num_bits_after_valid_syncro_sequence == 311
        @test state.data.SOW == SOW0 + 30
    end
end
