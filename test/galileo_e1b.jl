const GALILEO_E1B_DATA = [
    uint4000"0x580ba78311d174abc085cca0ec39253e13c2ff0fae80ff9d381d939d486fb29609400f2e5cffc233a00002ae1ffbd94e002fda08ff43ce8c020bd5cfffc5a580f10b21cc68ef3ad013f196a7c8f6692a750cfca2db1a9885e2b89c03f6dd60fc00f1c5bffc283c400032e6ffbc504002f3b20ff4052ec0204574fffcea5804243d38ef0fac918402d459eff55640602d1b9fff7b74870786bbe9e7729609400f5f56ffc383a00000dfaffbd680402f583eff407608020294bfffc925808a00162e77ffef510001a0ae7ffefcca002e6affff6e9ae00067f87ff27960cc00f1543ffc327d400026e2ffbce4f002fc606ff43c6e40206752ffff1a5810010000bffbffcd000000adff3fff1807000f9fe3fffd0060015bfffffb160e000f4347ffc387c000046ffffbd404002f872fff42ba00020d74ffffdea5813fbb0b5e04fdc05e2d4f89b1959ee136ec032b3bd3af3ff92d31176a1fd960b400f4c51ffc22f8c0004de5ffbd0cfc02f0d0cff40c2a8020935ffffe5258143b6667ebc63c1971394c6b98cbf40af9d331fbf5befac40037a0cf1d57d60dc00f2e55ffc2178c0006decffbc0ca802f9d0dff43c220020cb65fffde2581145b160daf5a411a06a19992fc2f151cd5701b7c0af9a040aa5e72d1cb4560d000f1344ffc34be40006ff9ffbd281c02f923dff43863c020235efffe0a";
    uint4000"0x580c00000ddfffffe8bffffe8e0000078c00001bdfffffec3fffff9f000002960f400f0d56ffc30f980005af1ffbd303802ffd2cff42f6d40205276fffef6580c00000ddfffffe8bffffe8e0000078c00001bdfffffec3fffff9f000002960fc00f3854ffc283900003cf0ffbd5c3402fde2fff4396dc020df74fffdb2580fb52424f71e672cde7837ed304eca8bda017fcb6f61190dda584bf4f6c4d60b400f724affc32fc800029f1ffbdf84802f3630ff4192f4020467efffc0e58022ba308c6177fdc491b380f573f49da98ba788f87401c257a03096011ffd60a000f0b4fffc353b000053e3ffbc7c0802f1f38ff41e284020b968fffe125801a7c507d17ccc31d023f3c8b7b3a5911ade11dd0fcf38d97862ae1b2699d60bc00f5850ffc3c38000054feffbcfc7402f1222ff401e08020695ffffdba580c00000ddfffffe8bffffe8e0000078c00001bdfffffec3fffff9f000002960f000f2b57ffc3ab980007df2ffbdfc2402fa72bff43d6cc0201476fffe3e581145b178daf5a411a06a12192fc2d951cd5779b7c0ac5a040aa6672d1c8a560c000f4346ffc3d7fc00006ffffbc800c02fdd39ff42e61c0208e56fffc96580ba78311d174abc085cca0ec39253e13c2ff0fae80ff9d381d939d486fb29609400f2e5cffc233a00002ae1ffbd94e002fda08ff43ce8c020bd5cfffc5a";
    uint4000"0x580f10b21cc68ef3ad013f196a7c8f6692a750cfca2db1a9885e2b89c03f6dd60fc00f1c5bffc283c400032e6ffbc504002f3b20ff4052ec0204574fffcea5804243d38ef0fac918402d459eff55640602d1b9fff7b74870786bbe9e7729609c00f2e54ffc293800003ef4ffbdec2c02fbc3eff422e300208745fffd02581cc40114b68ffbad9600208eebfefbd4d00a18ae3fd7fc2900818f91ff7fd609c00f7b58ffc3578000036eaffbd60d402f410cff407e900200c5dfffcee581fa91e3e8d439401b3625bff9bb0b04e9f1dfdc5d34a0c1020c7f8f2f000960f800f4350ffc223bc00051f1ffbc607402ff033ff43fabc0205a72fffe165814d4bd80ce0353ad6516d389220ec7c8131046e8ba2f1b91879f3508446d960f800f5250ffc2cbc400046f4ffbce8d802ff50cff4382040201360ffff52581150a07ad638faa8bae8a37f8b522907653a06af456fe79ca997c172875cd60a800f0c5cffc28bbc00052efffbc249402fa61cff434e80020e147fffd1a581145b000daf5a0b1a06a74192fc35551cd57e9b7c0b61a040a86672d1db6560c800f3442ffc3e7f40003dfdffbd083c02f1c35ff40ce1c020b75affffe2580c00000ddfffffe8bffffe8e0000078c00001bdfffffec3fffff9f000002960e000f4550ffc36f9400020f7ffbcb42402fd42dff4216d4020dd70fffd62";
]

@testset "Galileo E1B constructor" begin
    galileo_e1b = GalileoE1B()

    @test GalileoE1BDecoderState(21) == GNSSDecoderState(galileo_e1b, 21)
end

@testset "Galileo E1B test data decoding" begin
    decoder = GalileoE1BDecoderState(21)

    test_data = GNSSDecoder.GalileoE1BData(;
        WN = 1082,
        TOW = 259235,
        SVID = 21,
        t_0e = 259200,
        M_0 = 2.583259511699057,
        e = 9.999994654208422e-5,
        sqrt_A = 5440.588203430176,
        Ω_0 = 1.4358151884944024,
        i_0 = 0.9885709926241616,
        ω = 3.1416165752262324e-5,
        i_dot = 3.142988060925546e-11,
        Ω_dot = 3.141595145762181e-8,
        Δn = 3.142988060925545e-10,
        C_uc = 1.000240445137024e-6,
        C_us = 2.000480890274048e-6,
        C_rc = 5,
        C_rs = 6,
        C_ic = 3.0007213354110718e-6,
        C_is = 3.999099135398865e-6,
        SISA_e1_e5b = 25,
        t_0c = 259200,
        a_f0 = 9.999785106629133e-5,
        a_f1 = 1.0000036354540498e-9,
        a_f2 = 1.734723475976807e-18,
        IOD_nav1 = 0x0000000000000030,
        IOD_nav2 = 0x0000000000000030,
        IOD_nav3 = 0x0000000000000030,
        IOD_nav4 = 0x0000000000000030,
        num_pages_after_last_TOW = 12,
        num_bits_after_valid_syncro_sequence_after_last_TOW = 2510,
        signal_health_e1b = GNSSDecoder.signal_ok,
        signal_health_e5b = GNSSDecoder.signal_ok,
        data_validity_status_e1b = GNSSDecoder.navigation_data_valid,
        data_validity_status_e5b = GNSSDecoder.navigation_data_valid,
        broadcast_group_delay_e1_e5a = -9.313225746154785e-10,
        broadcast_group_delay_e1_e5b = -1.1641532182693481e-9,
        a_i0 = 100,
        a_i1 = 1,
        a_i2 = 0.100006103515625,
        iono_storm_flag_region1 = true,
        iono_storm_flag_region2 = false,
        iono_storm_flag_region3 = true,
        iono_storm_flag_region4 = false,
        iono_storm_flag_region5 = false,
        A_0_utc = 1.000240445137024e-6,
        A_1_utc = 8.881784197001252e-16,
        Δt_LS = 18,
        t_0t = 259200,
        WN_0t = 58,
        WN_LSF = 56,
        DN = 2,
        Δt_LSF = 18,
        A_0G = -2.9103830456733704e-11,
        A_1G = -4.440892098500626e-16,
        t_0G = 918000,
        WN_0G = 63,
        almanacs = Dictionary{Int,GNSSDecoder.GalileoAlmanac}(
            [19, 20, 21],
            [
                GNSSDecoder.GalileoAlmanac(;
                    SVID = 19,
                    Δsqrt_A = 0,
                    e = 0,
                    ω = 0,
                    δi = 0,
                    Ω_0 = 2.815621736164101,
                    Ω_dot = 0,
                    M_0 = 0.9006384700873591,
                    a_f0 = 0,
                    a_f1 = 0,
                    signal_health_e5b = GNSSDecoder.signal_ok,
                    signal_health_e1b = GNSSDecoder.signal_ok,
                    IOD_a = 0,
                    WN_a = 2,
                    t_0a = 259200,
                ),
                GNSSDecoder.GalileoAlmanac(;
                    SVID = 20,
                    Δsqrt_A = 0,
                    e = 0,
                    ω = 0,
                    δi = 0,
                    Ω_0 = 2.815621736164101,
                    Ω_dot = 0,
                    M_0 = 1.6860366334848091,
                    a_f0 = 0,
                    a_f1 = 0,
                    signal_health_e5b = GNSSDecoder.signal_ok,
                    signal_health_e1b = GNSSDecoder.signal_ok,
                    IOD_a = 0,
                    WN_a = 2,
                    t_0a = 259200,
                ),
                GNSSDecoder.GalileoAlmanac(;
                    SVID = 21,
                    Δsqrt_A = 0,
                    e = 0.0001068115234375,
                    ω = 0,
                    δi = 0.011121360712170923,
                    Ω_0 = 1.4358060174609635,
                    Ω_dot = 3.1452738704244004e-8,
                    M_0 = 2.5832236467994254,
                    a_f0 = 9.918212890625e-5,
                    a_f1 = 1.000444171950221e-9,
                    signal_health_e5b = GNSSDecoder.signal_ok,
                    signal_health_e1b = GNSSDecoder.signal_ok,
                    IOD_a = 0,
                    WN_a = 2,
                    t_0a = 259200,
                ),
            ],
        ),
        reduced_ced = GNSSDecoder.GalileoReducedCED(
            0,
            9.989738464355469e-5,
            0,
            0.0111865249350938,
            1.4358165036577555,
            2.5888697147579616,
            0.00010004639625549316,
            9.022187441587448e-10,
        )
    )

    # Convert the hard-bit chunks to ±1.0f0 soft symbols at the test boundary.
    # Galileo's decoder still consumes hard bits internally in this slice; the
    # boundary conversion mirrors how Tracking.jl v2 callers will feed soft
    # prompts.
    state = reduce(
        (dec, data) -> decode(dec, to_soft_symbols(data, sizeof(data) * 8), sizeof(data) * 8),
        GALILEO_E1B_DATA;
        init = decoder,
    )
    @test state.data == test_data
    @test state.is_shifted_by_180_degrees == false
    @test state.num_bits_after_valid_syncro_sequence == 3500
    @test is_sat_healthy(state) == true

    decoder2 = GalileoE1BDecoderState(21)
    state = reduce(
        (dec, data) -> decode(dec, to_soft_symbols(~data, sizeof(data) * 8), sizeof(data) * 8),
        GALILEO_E1B_DATA;
        init = decoder2,
    )
    @test state.data == test_data
    @test state.is_shifted_by_180_degrees == true
    @test state.num_bits_after_valid_syncro_sequence == 3500
    @test is_sat_healthy(state) == true

    state = reset_decoder_state(state)
    @test state.cache.packed_buffer[] == 0
    @test state.cache.complemented_buffer[] == 0
    @test length(state.cache.soft_buffer) == 0
    @test isnothing(state.raw_data.TOW)
    @test isnothing(state.data.TOW)
    @test GNSSDecoder.num_bits_buffered(state) == 0
    @test isnothing(state.num_bits_after_valid_syncro_sequence)
end