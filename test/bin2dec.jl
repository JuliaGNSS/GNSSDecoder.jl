@testset "Bin 2 Dec" begin
    test_signal_pos_1 = BitArray([1;0;0;0;1;0;0;0;1;1;0;1;0;0])
    test_signal_pos_2 = BitArray([1;0;1;0;0;1;1;0;0;1])

    @test @inferred(GNSSDecoder.bin2dec(test_signal_pos_1)) == 8756
    @test @inferred(GNSSDecoder.bin2dec(test_signal_pos_2)) == 665
end


@testset "Bin 2 Dec two's complement" begin
    test_signal_neg_8_bits = BitArray([1;1;0;1;1;1;1;1])
    test_signal_neg_40 =  BitArray([1;1;1;1;1;1;1;0;1;1;1;1;1;0;0;1;0;0;0;0;0;1;0;1;0;1;1;1;1;0;1;0;1;0;1;1;1;0;1;1])
    test_signal_pos_16_bits = BitArray([0;0;0;0;0;0;1;1;1;1;0;0;1;0;1;1])


    @test @inferred(GNSSDecoder.bin2dec_twoscomp(test_signal_neg_8_bits)) == -33
    @test @inferred(GNSSDecoder.bin2dec_twoscomp(test_signal_neg_40)) == -4412048709
    @test @inferred(GNSSDecoder.bin2dec_twoscomp(test_signal_pos_16_bits)) == 971

end
