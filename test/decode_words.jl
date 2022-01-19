
test_integrity_status_flag = false
test_TOW = 34945
test_alert_flag = false
test_anti_spoof_flag = true
test_trans_week = 67
test_codeonl2 = 1
test_ura =  2.0
test_svhealth = "000000"
test_IODC = "0001001000"
test_l2pcode = false
test_T_GD = -1.0710209608078003e-8
test_t_oc =  216000
test_a_f2 =  0.0
test_a_f1 =  -4.774847184307873e-12
test_a_f0 =  -0.00018549291417002678
test_IODE_Sub_2 =  "01001000"
test_C_rs =  70.65625
test_Δn =  3.930878022562108e-9
test_M_0 =  2.4393048719362045
test_C_uc =  3.604218363761902e-6
test_e =  0.01144192845094949
test_C_us =  1.3023614883422852e-5
test_sqrt_A =  5153.7995529174805
test_t_oe =  216000
test_fit_interval =  false
test_AODO = 31
test_C_ic =  -1.73225998878479e-7
test_Ω_0 =  0.0600607702978756
test_C_is =  -2.2351741790771484e-7
test_i_0 =  0.9781895349147778
test_C_rc =  136.34375
test_ω =  0.635978551768012
test_Ω_dot =  -7.383521839035659e-9
test_IODE_Sub_3 =  "01001000"
test_IDOT =  -3.4465721349922174e-10

@testset "Test Decoding true data" begin
    dc = deepcopy(DECODER)
    dc.prev_30 = PREV_30
    dc.prev_29 = PREV_29
    buff = deepcopy(BUFFER_TEST_TRUE)
    for i = 1:5
        dc = GNSSDecoder.decode_words(dc, buff[i])
        dc.prev_30 = buff[i][10][30]
        dc.prev_29 = buff[i][10][29]
    end
    
    @test dc.data.integrity_status_flag == test_integrity_status_flag #1
    @test dc.data.TOW == test_TOW #2
    @test dc.data.alert_flag == test_alert_flag #3
    @test dc.data.anti_spoof_flag == test_anti_spoof_flag #4
    @test dc.data.trans_week == test_trans_week #5
    @test dc.data.codeonl2 == test_codeonl2 #6
    @test dc.data.ura == test_ura #7
    @test dc.data.svhealth == test_svhealth #8
    @test dc.data.IODC == test_IODC#9
    @test dc.data.l2pcode == test_l2pcode #10
    @test dc.data.T_GD == test_T_GD #11
    @test dc.data.t_oc == test_t_oc #12
    @test dc.data.a_f2 == test_a_f2 #13
    @test dc.data.a_f1 == test_a_f1 #14
    @test dc.data.a_f0 == test_a_f0 #15
    @test dc.data.IODE_Sub_2 == test_IODE_Sub_2 #16
    @test dc.data.C_rs == test_C_rs #17
    @test dc.data.Δn == test_Δn #18
    @test dc.data.M_0 == test_M_0 #19
    @test dc.data.C_uc == test_C_uc #20
    @test dc.data.e == test_e #21
    @test dc.data.C_us == test_C_us #22
    @test dc.data.sqrt_A == test_sqrt_A #23
    @test dc.data.t_oe == test_t_oe #24
    @test dc.data.fit_interval == test_fit_interval #25
    @test dc.data.AODO == test_AODO #26
    @test dc.data.C_ic == test_C_ic #27
    @test dc.data.Ω_0 == test_Ω_0 #28
    @test dc.data.C_is == test_C_is #29
    @test dc.data.i_0 == test_i_0 #30
    @test dc.data.C_rc == test_C_rc #31
    @test dc.data.ω == test_ω #32
    @test dc.data.Ω_dot == test_Ω_dot #33
    @test dc.data.IODE_Sub_3 == test_IODE_Sub_3 #34
    @test dc.data.IDOT == test_IDOT #35
end



@testset "Test parity errors" begin

    PREV_29_FALSE = 1
    PREV_30_FALSE = 1
    dc = deepcopy(DECODER)
    @suppress begin
    dc.prev_30 = PREV_30_FALSE
    dc.prev_29 = PREV_29_FALSE
    buff = deepcopy(BUFFER_TEST_TRUE)
    for i = 1:5
        dc = GNSSDecoder.decode_words(dc, buff[i])
        dc.prev_30 = PREV_30_FALSE
        dc.prev_29 = PREV_30_FALSE
    end
    end
    @test dc.data_integrity == false

    
    dc = deepcopy(DECODER)
    @suppress begin  
        BAD_LINE = BitArray([1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0, 1, 1])
        dc.prev_30 = PREV_30
        dc.prev_29 = PREV_29
        buff = deepcopy(BUFFER_TEST_TRUE)
        buff[5][2] = BAD_LINE
        for i = 1:5
            dc = GNSSDecoder.decode_words(dc, buff[i])
            dc.prev_30 = buff[i][10][30]
            dc.prev_29 = buff[i][10][29]
        end
    end

    @test dc.data_integrity == false
    @test dc.data.TOW != test_TOW
end


    

    


@testset "Test Decoding CEI-cutovers" begin


    IODC_CHANGE_LSB = BitArray([1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0, 1, 1])
    dc = deepcopy(DECODER)
    dc.prev_30 = PREV_30
    dc.prev_29 = PREV_29
    buff = deepcopy(BUFFER_TEST_TRUE)
    buff[1][8] = IODC_CHANGE_LSB
    out = true
    @suppress begin
        TLM_HOW_data,sub_1_data = GNSSDecoder.decode_subframe_1(buff[1])
        TLM_HOW_data,sub_2_data = GNSSDecoder.decode_subframe_2(buff[2])
        TLM_HOW_data,sub_3_data = GNSSDecoder.decode_subframe_3(buff[3])
    
        out = GNSSDecoder.control_data(TLM_HOW_data, sub_1_data, sub_2_data, sub_3_data)
    end

    @test out == false

    out = true

    IODE_SUB2_CHANGE = BitArray([1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 1, 0])
    dc = deepcopy(DECODER)
    dc.prev_30 = PREV_30
    dc.prev_29 = PREV_29
    buff = deepcopy(BUFFER_TEST_TRUE)
    buff[2][3] = IODE_SUB2_CHANGE
    @suppress begin
        TLM_HOW_data,sub_1_data = GNSSDecoder.decode_subframe_1(buff[1])
        TLM_HOW_data,sub_2_data = GNSSDecoder.decode_subframe_2(buff[2])
        TLM_HOW_data,sub_3_data = GNSSDecoder.decode_subframe_3(buff[3])
    
    out = GNSSDecoder.control_data(TLM_HOW_data, sub_1_data, sub_2_data, sub_3_data)
    end
    @test out == false




    out = true
    IODE_SUB3_CHANGE = BitArray([1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0])
    dc = deepcopy(DECODER)
    dc.prev_30 = PREV_30
    dc.prev_29 = PREV_29
    buff = deepcopy(BUFFER_TEST_TRUE)
    buff[3][10] = IODE_SUB3_CHANGE
    
    @suppress begin
        TLM_HOW_data,sub_1_data = GNSSDecoder.decode_subframe_1(buff[1])
        TLM_HOW_data,sub_2_data = GNSSDecoder.decode_subframe_2(buff[2])
        TLM_HOW_data,sub_3_data = GNSSDecoder.decode_subframe_3(buff[3])
        out = GNSSDecoder.control_data(TLM_HOW_data, sub_1_data, sub_2_data, sub_3_data)
    end
    @test out == false
    
end