@testset "Bit fiddling" begin
    
    @test GNSSDecoder.get_bits(0b111, 4, 2, 3) == 0b111
    @test GNSSDecoder.get_bits(0b111, 5, 3, 3) == 0b111
    @test GNSSDecoder.get_bits(0b111, 5, 3, 2) == 0b11
    @test GNSSDecoder.get_bits(0b111, 4, 2, 2) == 0b11
    @test GNSSDecoder.get_bits(0b101, 4, 2, 2) == 0b10
    @test GNSSDecoder.get_bits(0b101, 3, 1, 2) == 0b10
    @test GNSSDecoder.get_bits(0b101, 3, 1, 3) == 0b101

    @test GNSSDecoder.get_bit(0b101, 3, 1) == true
    @test GNSSDecoder.get_bit(0b101, 3, 2) == false
    @test GNSSDecoder.get_bit(0b101, 3, 3) == true
    @test GNSSDecoder.get_bit(0b101, 4, 2) == true
    @test GNSSDecoder.get_bit(0b101, 4, 3) == false
    @test GNSSDecoder.get_bit(0b101, 4, 4) == true

    @test GNSSDecoder.get_two_complement_num(0b000, 3, 1, 3) == 0
    @test GNSSDecoder.get_two_complement_num(0b001, 3, 1, 3) == 1
    @test GNSSDecoder.get_two_complement_num(0b011, 3, 1, 3) == 3
    @test GNSSDecoder.get_two_complement_num(0b111, 3, 1, 3) == -1
    @test GNSSDecoder.get_two_complement_num(0b100, 3, 1, 3) == -4
    @test GNSSDecoder.get_two_complement_num(0b101, 3, 1, 3) == -3
end
