@testset "Test Parity Buffer Control" begin
    
    dc = deepcopy(DECODER)
    @suppress begin
        buffer = deepcopy(BUFFER_TEST_TRUE)
        dc = GNSSDecoder.buffer_control(dc, buffer[1])
    end
    @test dc.data_integrity == true
    
    dc = deepcopy(DECODER)
    @suppress begin
        dc = deepcopy(DECODER)
        buffer = deepcopy(BUFFER_TEST_TRUE)
        dc.prev_30 = 1
        dc = GNSSDecoder.buffer_control(dc, buffer[1])
    end
    @test dc.data_integrity == false

    dc = deepcopy(DECODER)
    @suppress begin
        dc = deepcopy(DECODER)
        BUFFER_TEST_FALSE = deepcopy(BUFFER_TEST_TRUE)
        BUFFER_TEST_FALSE[2][3][11] = !BUFFER_TEST_FALSE[2][3][11]
        dc = GNSSDecoder.buffer_control(dc, BUFFER_TEST_FALSE[2])
    end
    @test dc.data_integrity == false
end