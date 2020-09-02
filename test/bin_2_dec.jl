@testset "Binary to Decimal" begin
    
    arr = BitArray([1,0,0,0,0])
    out = GNSSDecoder.bin2dec(arr)
    @test out == 16

    arr = BitArray([1,0,1,0,0])
    out = GNSSDecoder.bin2dec(arr)
    @test out == 20

    arr = BitArray([1,0,0,0,0,0,0,1])
    out = GNSSDecoder.bin2dec(arr)
    @test out == 129
end


@testset "Binary to Decimal Two's complement" begin

    arr = BitArray([1,0,0,0,0])
    out = GNSSDecoder.bin2dec_twoscomp(arr)
    @test out == -16

    arr = BitArray([1,0,1,0,0])
    out = GNSSDecoder.bin2dec_twoscomp(arr)
    @test out == -12

    arr = BitArray([1,0,0,0,0,0,0,1])
    out = GNSSDecoder.bin2dec_twoscomp(arr)
    @test out == -127
end
