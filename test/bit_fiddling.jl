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

    @test GNSSDecoder.get_twos_complement_num(0b000, 3, 1, 3) == 0
    @test GNSSDecoder.get_twos_complement_num(0b001, 3, 1, 3) == 1
    @test GNSSDecoder.get_twos_complement_num(0b011, 3, 1, 3) == 3
    @test GNSSDecoder.get_twos_complement_num(0b111, 3, 1, 3) == -1
    @test GNSSDecoder.get_twos_complement_num(0b100, 3, 1, 3) == -4
    @test GNSSDecoder.get_twos_complement_num(0b101, 3, 1, 3) == -3

    # Exemplary bits from the Galileo specification, run through the same
    # soft-symbol FEC chain the E1B decoder uses: ±1 Float32 LLRs ('0' ⇒ +1,
    # '1' ⇒ -1; positive ⇒ bit 0) → 30×8 block deinterleave → negate every
    # second symbol → AFF3CT K=7 NSC Viterbi.
    to_soft(s) = Float32[c == '1' ? -1.0f0 : 1.0f0 for c in s]
    to_bits(v) = join(ifelse.(v .< 0, '1', '0'))

    ex_encoded_bits = "101000000101111110001100111100000101111010100000000110011110001110101000010100000100111101010111011110001000011011111010111001111001100000011111111000100000100111110110000010011100011101100000100101110100100011000110110110010000011100111010"
    ex_deinterleaved_bits = "100011000001101010101010011100110011000101011010011011110101100101111000100101010101010110001100110011101010010110010000101001101000010000000010000000100001101110011001111010110101110000011000101110111111110111111101111001000101001100100010"
    deinterleaved = GNSSDecoder.deinterleave(to_soft(ex_encoded_bits), 30, 8)
    @test to_bits(deinterleaved) == ex_deinterleaved_bits

    ex_inv_deinterleaved_bits = "110110010100111111111111001001100110010000001111001110100000110000101101110000000000000011011001100110111111000011000101111100111101000101010111010101110100111011001100101111100000100101001101111011101010100010101000101100010000011001110111"
    # Inverting a hard bit is negating the soft symbol.
    inverted = copy(deinterleaved)
    inverted[2:2:end] .= -inverted[2:2:end]
    @test to_bits(inverted) == ex_inv_deinterleaved_bits

    ex_true_bits = "111111111111000011001100101010100000000000001111001100110101010111100011111011001101111110001010000111000001001101000000"
    # The K=7 NSC Viterbi step returns the 114 information bits (the 6 tail
    # bits are absorbed by trellis termination).
    dec = Aff3ct.ConvViterbiDecoder(114, 240, [0o171, 0o133])
    decoded = join(string.(Int.(Aff3ct.decode(dec, inverted))))
    @test decoded == ex_true_bits[1:114]
end
