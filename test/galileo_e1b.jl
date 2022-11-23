BitIntegers.@define_integers 7520
const GALILEO_E1B_DATA =
    uint7520"0xffea400000b780000261fffff953e562a175eec2f8024501fdf196b047f3c7e186a0787c8f92d3f4f6bca06174fe7ffffe04000003e80000223fffff467ffffea400000b780000261fffff953e47501d3ea73f9e259eff711af5040f3c380c4a28479ff9d09fd56cd0835f4fe7ffffe04000003e80000223fffff467ffffea400000b780000261fffff953e49ca15be1fa79e277cffad182e04f725981eea951fc87cb47f9d7477053f4fe7ffffe04000003e80000223fffff467ffffea400000b780000261fffff953e55001fde50c7aaa6005f5d1b7105173eda116ac037ee7c24bfd57631056f4fe56155760a21263eefb50a07414bf3e4e1a382663d7ec5143e4baa2af24553ebaac4bdd8480b7a8cf5ac12c29b8361163a02f3991246f8b70f14e1af95f4ffcc27e40658eac47e7efb103c05584467e4adcec50490f6ceb30af5ef54253e325965fe54cb87addcae730a7bb6d62f5140cc0aaf27e2933dbd7aa4baff4ff973c0786866e3bc15330d6c0b5cc4dd47406c584373cf274fb8a7e4d88053e863010bdd79c512f473995025f798f1f5edacbba52f17fcb3247448a53a74fe7ffffe04000003e80000223fffff467ffffea400000b780000261fffff953e44f9179e129de3272be2fb1c6ed51710a0d68a5978dd7ecb57e373e75d274fe7ffffe04000003e80000223fffff467ffffea400000b780000261fffff953e47fe1abe4487a12785ffe11bc885571d3a0bea7a17adfc33be8b79b005774feeb312644db4c5166d0ada45be4608d318479ab42bdfd52b2fd7e70a26b753e5c2e1c5dd5b7b9a367ffa9371406076f7c1eeb9877ea61cc1f8d013984cf4fe5ee27a26245b0bd75e5ff09280d4352de9e647845e9e42cf7320f4c2b66d3e988a143e86d7a72aa93f691f5105a6092814ae1e2fc4f4c7bec32d6601574fe00799f25fa0e59c81fa4349fc113ae00fa8a2ff815c9602f1cbeff0340453ef9ec101d800782a105df65047b870e11c4158c6727c2e2575f1d54c986ef4fc2fffc9e8e002884a8012023dff87edcfffb08880053df7001d8a0bffdef53ebaac121c7ccf91a92eff4d1ea30557b91c0228286f89fc4ffee7751605774fc6df8a9cea05d1c45185080230ff4eea30f3a6776db937b62e1b2ecda31653ef65217bd906fadaa919f2b0295062fc4ae09280107ba7c2a3f416ffa008f4fe7ffffe04000003e80000223fffff467ffffea400000b780000261fffff953e4cf0191e3b57a824ef1f0b1746060f2e6a0f2a7cb7857cbbbeb57919076f4fe7ffffe04000003e80000223fffff467ffffea4"

@testset "Galileo E1B constructor" begin
    galileo_e1b = GalileoE1B()

    @test GalileoE1BDecoderState(21) == GNSSDecoderState(galileo_e1b, 21)
end

@testset "Galileo E1B test data decoding" begin
    decoder = GalileoE1BDecoderState(21)

    test_data = GNSSDecoder.GalileoE1BData(;
        WN = 1170,
        TOW = 558041,
        t_0e = 556800.0,
        M_0 = 0.03588957467055676,
        e = 0.0002609505318105221,
        sqrt_A = 5440.625303268433,
        Ω_0 = 2.5016390578909675,
        i_0 = 0.9744856179803416,
        ω = 0.400715796701369,
        i_dot = -6.482412875658937e-10,
        Ω_dot = -5.1184274887640895e-9,
        Δn = 2.60617998642883e-9,
        C_uc = 5.239620804786682e-6,
        C_us = 1.1917203664779663e-5,
        C_rc = 90.4375,
        C_rs = 113.6875,
        C_ic = 7.450580596923828e-9,
        C_is = 0.0,
        t_0c = 556800.0,
        a_f0 = -0.0007140382076613605,
        a_f1 = -2.2311041902867146e-12,
        a_f2 = 0.0,
        IOD_nav1 = 0x0000000000000020,
        IOD_nav2 = 0x0000000000000020,
        IOD_nav3 = 0x0000000000000020,
        IOD_nav4 = 0x0000000000000020,
        num_pages_after_last_TOW = 1,
        signal_health_e1b = GNSSDecoder.signal_ok,
        signal_health_e5b = GNSSDecoder.signal_ok,
        data_validity_status_e1b = GNSSDecoder.navigation_data_valid,
        data_validity_status_e5b = GNSSDecoder.navigation_data_valid,
        broadcast_group_delay_e1_e5a = 2.7939677238464355e-9,
        broadcast_group_delay_e1_e5b = 3.026798367500305e-9,
    )

    state = decode(decoder, GALILEO_E1B_DATA, 7000)
    @test state.data == test_data
    @test state.is_shifted_by_180_degrees == true
    @test state.num_bits_after_valid_syncro_sequence == 665
    @test is_sat_healthy(state) == true

    state = decode(decoder, ~GALILEO_E1B_DATA, 7000)
    @test state.data == test_data
    @test state.is_shifted_by_180_degrees == false
    @test state.num_bits_after_valid_syncro_sequence == 665
    @test is_sat_healthy(state) == true
end