@testset "Test Parity Buffer Control" begin
    
    dc = GNSSDecoderState(PRN = 1)
    dc.prev_30 = PREV_30
    dc.prev_29 = PREV_29
    buff = deepcopy(BUFFER_TEST_TRUE)
    dc = GNSSDecoder.check_buffer(dc, buff[1])
    @test dc.data_integrity == true
    
    dc = GNSSDecoderState(PRN = 1)
    @suppress begin
        dc.prev_30 = 1
        dc = GNSSDecoder.check_buffer(dc, BUFFER_TEST_TRUE[1])
    end
    @test dc.data_integrity == false

    dc = GNSSDecoderState(PRN = 1)
    @suppress begin
        BUFFER_TEST_FALSE = deepcopy(BUFFER_TEST_TRUE)
        BUFFER_TEST_FALSE[2][3][11] = !BUFFER_TEST_FALSE[2][3][11]
        dc = GNSSDecoder.check_buffer(dc, BUFFER_TEST_FALSE[2])
    end
    @test dc.data_integrity == false
end