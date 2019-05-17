##All words from subframe 2
@testset "word1sub2" begin

    decode = GNSSDecoder.init_decode()
    _interm = GNSSDecoder.GPSData_interm()
    decode.parameters1.prev_30 = false

    #words
    word1 = [true,false,false,false,true,false,true,true,false,false,false,false,true,true,true,true,true,false,false,false,false,false,false,false,true,false,false,true,true,false]
    word1 = word1[30:-1:1]

    #first word
    decode.parameters1.word_count = decode.parameters1.word_count + 1
    GNSSDecoder.decodeword(word1,decode.data,decode.parameters1,_interm)

    @test decode.parameters1.word_count == 1
    @test decode.parameters1.prev_29 == true
    @test decode.parameters1.prev_30 == false


end

@testset "word2sub2" begin

    decode = GNSSDecoder.init_decode()
    _interm = GNSSDecoder.GPSData_interm()
    decode.parameters1.prev_30 = false

    #words
    word1 = [true,false,false,false,true,false,true,true,false,false,false,false,true,true,true,true,true,false,false,false,false,false,false,false,true,false,false,true,true,false]
    word1 = word1[30:-1:1]

    #first word
    decode.parameters1.word_count = decode.parameters1.word_count + 1
    GNSSDecoder.decodeword(word1,decode.data,decode.parameters1,_interm)

    #second word
    word2 = [true,false,false,true,true,false,true,true,true,true,true,false,true,false,false,true,true,false,true,false,true,false,false,false,true,false,true,false,false,false]
    word2 = word2[30:-1:1]

    decode.parameters1.word_count = decode.parameters1.word_count + 1
    GNSSDecoder.decodeword(word2,decode.data,decode.parameters1,_interm)

    @test decode.parameters1.word_count == 2
    @test decode.parameters1.prev_29 == false
    @test decode.parameters1.prev_30 == false
    @test decode.parameters1.subframe_count == 2

end

@testset "word3sub2" begin
    #initialization
    decode = GNSSDecoder.init_decode()
    _interm = GNSSDecoder.GPSData_interm()
    decode.parameters1.prev_30 = false

    #first word
    word1 = [true,false,false,false,true,false,true,true,false,false,false,false,true,true,true,true,true,false,false,false,false,false,false,false,true,false,false,true,true,false]
    word1 = word1[30:-1:1]

    decode.parameters1.word_count = decode.parameters1.word_count + 1
    GNSSDecoder.decodeword(word1,decode.data,decode.parameters1,_interm)

    #second word
    word2 = [true,false,false,true,true,false,true,true,true,true,true,false,true,false,false,true,true,false,true,false,true,false,false,false,true,false,true,false,false,false]
    word2 = word2[30:-1:1]

    decode.parameters1.word_count = decode.parameters1.word_count + 1
    GNSSDecoder.decodeword(word2,decode.data,decode.parameters1,_interm)

    #third word
    word3 = [false,true,true,true,true,false,true,true,false,false,false,false,false,false,true,true,true,true,false,false,true,false,true,true,false,false,true,true,true,false]
    word3 = word3[30:-1:1]

    decode.parameters1.word_count = decode.parameters1.word_count + 1
    GNSSDecoder.decodeword(word3,decode.data,decode.parameters1,_interm)

    @test decode.parameters1.word_count == 3
    @test decode.parameters1.prev_29 == true
    @test decode.parameters1.prev_30 == false

    @test decode.data.IODE == "01111011"
    @test decode.data.C_rs ≈ 30.34375
end

@testset "word4sub2" begin
    #initialization
    decode = GNSSDecoder.init_decode()
    _interm = GNSSDecoder.GPSData_interm()
    decode.parameters1.prev_30 = false

    #first word
    word1 = [true,false,false,false,true,false,true,true,false,false,false,false,true,true,true,true,true,false,false,false,false,false,false,false,true,false,false,true,true,false]
    word1 = word1[30:-1:1]

    decode.parameters1.word_count = decode.parameters1.word_count + 1
    GNSSDecoder.decodeword(word1,decode.data,decode.parameters1,_interm)

    #second word
    word2 = [true,false,false,true,true,false,true,true,true,true,true,false,true,false,false,true,true,false,true,false,true,false,false,false,true,false,true,false,false,false]
    word2 = word2[30:-1:1]

    decode.parameters1.word_count = decode.parameters1.word_count + 1
    GNSSDecoder.decodeword(word2,decode.data,decode.parameters1,_interm)

    #third word
    word3 = [false,true,true,true,true,false,true,true,false,false,false,false,false,false,true,true,true,true,false,false,true,false,true,true,false,false,true,true,true,false]
    word3 = word3[30:-1:1]

    decode.parameters1.word_count = decode.parameters1.word_count + 1
    GNSSDecoder.decodeword(word3,decode.data,decode.parameters1,_interm)

    #fourth word
    word4 = [false,false,true,true,false,false,true,false,true,true,true,true,true,false,false,false,false,false,true,false,false,true,true,true,true,false,false,false,false,false]
    word4 = word4[30:-1:1]

    decode.parameters1.word_count = decode.parameters1.word_count + 1
    GNSSDecoder.decodeword(word4,decode.data,decode.parameters1,_interm)

    @test decode.parameters1.word_count == 4
    @test decode.parameters1.prev_29 == false
    @test decode.parameters1.prev_30 == false

    @test decode.data.Δn ≈ 4.660194116*10^-9

end

@testset "word5sub2" begin
    #initialization
    decode = GNSSDecoder.init_decode()
    _interm = GNSSDecoder.GPSData_interm()
    decode.parameters1.prev_30 = false

    #first word
    word1 = [true,false,false,false,true,false,true,true,false,false,false,false,true,true,true,true,true,false,false,false,false,false,false,false,true,false,false,true,true,false]
    word1 = word1[30:-1:1]

    decode.parameters1.word_count = decode.parameters1.word_count + 1
    GNSSDecoder.decodeword(word1,decode.data,decode.parameters1,_interm)

    #second word
    word2 = [true,false,false,true,true,false,true,true,true,true,true,false,true,false,false,true,true,false,true,false,true,false,false,false,true,false,true,false,false,false]
    word2 = word2[30:-1:1]

    decode.parameters1.word_count = decode.parameters1.word_count + 1
    GNSSDecoder.decodeword(word2,decode.data,decode.parameters1,_interm)

    #third word
    word3 = [false,true,true,true,true,false,true,true,false,false,false,false,false,false,true,true,true,true,false,false,true,false,true,true,false,false,true,true,true,false]
    word3 = word3[30:-1:1]

    decode.parameters1.word_count = decode.parameters1.word_count + 1
    GNSSDecoder.decodeword(word3,decode.data,decode.parameters1,_interm)

    #fourth word
    word4 = [false,false,true,true,false,false,true,false,true,true,true,true,true,false,false,false,false,false,true,false,false,true,true,true,true,false,false,false,false,false]
    word4 = word4[30:-1:1]

    decode.parameters1.word_count = decode.parameters1.word_count + 1
    GNSSDecoder.decodeword(word4,decode.data,decode.parameters1,_interm)

    #fifth word
    word5 = [true,true,true,false,false,true,false,false,false,true,true,false,true,true,false,true,true,true,false,false,false,false,false,true,false,false,true,false,true,false]
    word5 = word5[30:-1:1]

    decode.parameters1.word_count = decode.parameters1.word_count + 1
    GNSSDecoder.decodeword(word5,decode.data,decode.parameters1,_interm)

    @test decode.parameters1.word_count == 5
    @test decode.parameters1.prev_29 == true
    @test decode.parameters1.prev_30 == false

    @test decode.data.M0 ≈ 0.9791043415

end

@testset "word6sub2" begin
    #initialization
    decode = GNSSDecoder.init_decode()
    _interm = GNSSDecoder.GPSData_interm()
    decode.parameters1.prev_30 = false

    #first word
    word1 = [true,false,false,false,true,false,true,true,false,false,false,false,true,true,true,true,true,false,false,false,false,false,false,false,true,false,false,true,true,false]
    word1 = word1[30:-1:1]

    decode.parameters1.word_count = decode.parameters1.word_count + 1
    GNSSDecoder.decodeword(word1,decode.data,decode.parameters1,_interm)

    #second word
    word2 = [true,false,false,true,true,false,true,true,true,true,true,false,true,false,false,true,true,false,true,false,true,false,false,false,true,false,true,false,false,false]
    word2 = word2[30:-1:1]

    decode.parameters1.word_count = decode.parameters1.word_count + 1
    GNSSDecoder.decodeword(word2,decode.data,decode.parameters1,_interm)

    #third word
    word3 = [false,true,true,true,true,false,true,true,false,false,false,false,false,false,true,true,true,true,false,false,true,false,true,true,false,false,true,true,true,false]
    word3 = word3[30:-1:1]

    decode.parameters1.word_count = decode.parameters1.word_count + 1
    GNSSDecoder.decodeword(word3,decode.data,decode.parameters1,_interm)

    #fourth word
    word4 = [false,false,true,true,false,false,true,false,true,true,true,true,true,false,false,false,false,false,true,false,false,true,true,true,true,false,false,false,false,false]
    word4 = word4[30:-1:1]

    decode.parameters1.word_count = decode.parameters1.word_count + 1
    GNSSDecoder.decodeword(word4,decode.data,decode.parameters1,_interm)

    #fifth word
    word5 = [true,true,true,false,false,true,false,false,false,true,true,false,true,true,false,true,true,true,false,false,false,false,false,true,false,false,true,false,true,false]
    word5 = word5[30:-1:1]

    decode.parameters1.word_count = decode.parameters1.word_count + 1
    GNSSDecoder.decodeword(word5,decode.data,decode.parameters1,_interm)

    #sixth word
    word6 = [false,false,false,false,false,false,true,true,true,false,true,true,true,false,false,false,false,false,false,false,false,true,false,false,false,false,false,true,true,false]
    word6 = word6[30:-1:1]

    decode.parameters1.word_count = decode.parameters1.word_count + 1
    GNSSDecoder.decodeword(word6,decode.data,decode.parameters1,_interm)

    @test decode.parameters1.word_count == 6
    @test decode.parameters1.prev_29 == true
    @test decode.parameters1.prev_30 == false

    @test decode.data.C_uc ≈ 1.773238182*10^-6

end

@testset "word7sub2" begin
    #initialization
    decode = GNSSDecoder.init_decode()
    _interm = GNSSDecoder.GPSData_interm()
    decode.parameters1.prev_30 = false

    #first word
    word1 = [true,false,false,false,true,false,true,true,false,false,false,false,true,true,true,true,true,false,false,false,false,false,false,false,true,false,false,true,true,false]
    word1 = word1[30:-1:1]

    decode.parameters1.word_count = decode.parameters1.word_count + 1
    GNSSDecoder.decodeword(word1,decode.data,decode.parameters1,_interm)

    #second word
    word2 = [true,false,false,true,true,false,true,true,true,true,true,false,true,false,false,true,true,false,true,false,true,false,false,false,true,false,true,false,false,false]
    word2 = word2[30:-1:1]

    decode.parameters1.word_count = decode.parameters1.word_count + 1
    GNSSDecoder.decodeword(word2,decode.data,decode.parameters1,_interm)

    #third word
    word3 = [false,true,true,true,true,false,true,true,false,false,false,false,false,false,true,true,true,true,false,false,true,false,true,true,false,false,true,true,true,false]
    word3 = word3[30:-1:1]

    decode.parameters1.word_count = decode.parameters1.word_count + 1
    GNSSDecoder.decodeword(word3,decode.data,decode.parameters1,_interm)

    #fourth word
    word4 = [false,false,true,true,false,false,true,false,true,true,true,true,true,false,false,false,false,false,true,false,false,true,true,true,true,false,false,false,false,false]
    word4 = word4[30:-1:1]

    decode.parameters1.word_count = decode.parameters1.word_count + 1
    GNSSDecoder.decodeword(word4,decode.data,decode.parameters1,_interm)

    #fifth word
    word5 = [true,true,true,false,false,true,false,false,false,true,true,false,true,true,false,true,true,true,false,false,false,false,false,true,false,false,true,false,true,false]
    word5 = word5[30:-1:1]

    decode.parameters1.word_count = decode.parameters1.word_count + 1
    GNSSDecoder.decodeword(word5,decode.data,decode.parameters1,_interm)

    #sixth word
    word6 = [false,false,false,false,false,false,true,true,true,false,true,true,true,false,false,false,false,false,false,false,false,true,false,false,false,false,false,true,true,false]
    word6 = word6[30:-1:1]

    decode.parameters1.word_count = decode.parameters1.word_count + 1
    GNSSDecoder.decodeword(word6,decode.data,decode.parameters1,_interm)

    #seventh word
    word7 = [false,false,true,false,false,true,false,true,true,true,true,true,true,true,true,true,false,false,false,false,true,false,false,false,false,false,false,false,false,true]
    word7 = word7[30:-1:1]

    decode.parameters1.word_count = decode.parameters1.word_count + 1
    GNSSDecoder.decodeword(word7,decode.data,decode.parameters1,_interm)

    @test decode.parameters1.word_count == 7
    @test decode.parameters1.prev_29 == false
    @test decode.parameters1.prev_30 == true

    @test decode.data.e ≈ 8.102388121*10^-3

end

@testset "word8sub2" begin
    #initialization
    decode = GNSSDecoder.init_decode()
    _interm = GNSSDecoder.GPSData_interm()
    decode.parameters1.prev_30 = false

    #first word
    word1 = [true,false,false,false,true,false,true,true,false,false,false,false,true,true,true,true,true,false,false,false,false,false,false,false,true,false,false,true,true,false]
    word1 = word1[30:-1:1]

    decode.parameters1.word_count = decode.parameters1.word_count + 1
    GNSSDecoder.decodeword(word1,decode.data,decode.parameters1,_interm)

    #second word
    word2 = [true,false,false,true,true,false,true,true,true,true,true,false,true,false,false,true,true,false,true,false,true,false,false,false,true,false,true,false,false,false]
    word2 = word2[30:-1:1]

    decode.parameters1.word_count = decode.parameters1.word_count + 1
    GNSSDecoder.decodeword(word2,decode.data,decode.parameters1,_interm)

    #third word
    word3 = [false,true,true,true,true,false,true,true,false,false,false,false,false,false,true,true,true,true,false,false,true,false,true,true,false,false,true,true,true,false]
    word3 = word3[30:-1:1]

    decode.parameters1.word_count = decode.parameters1.word_count + 1
    GNSSDecoder.decodeword(word3,decode.data,decode.parameters1,_interm)

    #fourth word
    word4 = [false,false,true,true,false,false,true,false,true,true,true,true,true,false,false,false,false,false,true,false,false,true,true,true,true,false,false,false,false,false]
    word4 = word4[30:-1:1]

    decode.parameters1.word_count = decode.parameters1.word_count + 1
    GNSSDecoder.decodeword(word4,decode.data,decode.parameters1,_interm)

    #fifth word
    word5 = [true,true,true,false,false,true,false,false,false,true,true,false,true,true,false,true,true,true,false,false,false,false,false,true,false,false,true,false,true,false]
    word5 = word5[30:-1:1]

    decode.parameters1.word_count = decode.parameters1.word_count + 1
    GNSSDecoder.decodeword(word5,decode.data,decode.parameters1,_interm)

    #sixth word
    word6 = [false,false,false,false,false,false,true,true,true,false,true,true,true,false,false,false,false,false,false,false,false,true,false,false,false,false,false,true,true,false]
    word6 = word6[30:-1:1]

    decode.parameters1.word_count = decode.parameters1.word_count + 1
    GNSSDecoder.decodeword(word6,decode.data,decode.parameters1,_interm)

    #seventh word
    word7 = [false,false,true,false,false,true,false,true,true,true,true,true,true,true,true,true,false,false,false,false,true,false,false,false,false,false,false,false,false,true]
    word7 = word7[30:-1:1]

    decode.parameters1.word_count = decode.parameters1.word_count + 1
    GNSSDecoder.decodeword(word7,decode.data,decode.parameters1,_interm)

    #eighth word
    word8 = [true,true,true,true,false,false,true,false,true,true,true,false,false,false,false,true,false,true,false,true,true,true,true,false,false,false,true,false,true,true]
    word8 = word8[30:-1:1]

    decode.parameters1.word_count = decode.parameters1.word_count + 1
    GNSSDecoder.decodeword(word8,decode.data,decode.parameters1,_interm)

    @test decode.parameters1.word_count == 8
    @test decode.parameters1.prev_29 == true
    @test decode.parameters1.prev_30 == true

    @test decode.data.C_us ≈ 6.254762411*10^-6

end

@testset "word9sub2" begin
    #initialization
    decode = GNSSDecoder.init_decode()
    _interm = GNSSDecoder.GPSData_interm()
    decode.parameters1.prev_30 = false

    #first word
    word1 = [true,false,false,false,true,false,true,true,false,false,false,false,true,true,true,true,true,false,false,false,false,false,false,false,true,false,false,true,true,false]
    word1 = word1[30:-1:1]

    decode.parameters1.word_count = decode.parameters1.word_count + 1
    GNSSDecoder.decodeword(word1,decode.data,decode.parameters1,_interm)

    #second word
    word2 = [true,false,false,true,true,false,true,true,true,true,true,false,true,false,false,true,true,false,true,false,true,false,false,false,true,false,true,false,false,false]
    word2 = word2[30:-1:1]

    decode.parameters1.word_count = decode.parameters1.word_count + 1
    GNSSDecoder.decodeword(word2,decode.data,decode.parameters1,_interm)

    #third word
    word3 = [false,true,true,true,true,false,true,true,false,false,false,false,false,false,true,true,true,true,false,false,true,false,true,true,false,false,true,true,true,false]
    word3 = word3[30:-1:1]

    decode.parameters1.word_count = decode.parameters1.word_count + 1
    GNSSDecoder.decodeword(word3,decode.data,decode.parameters1,_interm)

    #fourth word
    word4 = [false,false,true,true,false,false,true,false,true,true,true,true,true,false,false,false,false,false,true,false,false,true,true,true,true,false,false,false,false,false]
    word4 = word4[30:-1:1]

    decode.parameters1.word_count = decode.parameters1.word_count + 1
    GNSSDecoder.decodeword(word4,decode.data,decode.parameters1,_interm)

    #fifth word
    word5 = [true,true,true,false,false,true,false,false,false,true,true,false,true,true,false,true,true,true,false,false,false,false,false,true,false,false,true,false,true,false]
    word5 = word5[30:-1:1]

    decode.parameters1.word_count = decode.parameters1.word_count + 1
    GNSSDecoder.decodeword(word5,decode.data,decode.parameters1,_interm)

    #sixth word
    word6 = [false,false,false,false,false,false,true,true,true,false,true,true,true,false,false,false,false,false,false,false,false,true,false,false,false,false,false,true,true,false]
    word6 = word6[30:-1:1]

    decode.parameters1.word_count = decode.parameters1.word_count + 1
    GNSSDecoder.decodeword(word6,decode.data,decode.parameters1,_interm)

    #seventh word
    word7 = [false,false,true,false,false,true,false,true,true,true,true,true,true,true,true,true,false,false,false,false,true,false,false,false,false,false,false,false,false,true]
    word7 = word7[30:-1:1]

    decode.parameters1.word_count = decode.parameters1.word_count + 1
    GNSSDecoder.decodeword(word7,decode.data,decode.parameters1,_interm)

    #eighth word
    word8 = [true,true,true,true,false,false,true,false,true,true,true,false,false,false,false,true,false,true,false,true,true,true,true,false,false,false,true,false,true,true]
    word8 = word8[30:-1:1]

    decode.parameters1.word_count = decode.parameters1.word_count + 1
    GNSSDecoder.decodeword(word8,decode.data,decode.parameters1,_interm)

    #ninth word
    word9 = [true,true,true,true,false,false,true,false,true,false,true,true,false,false,true,false,true,true,true,true,true,true,false,true,true,false,false,false,false,true]
    word9 = word9[30:-1:1]

    decode.parameters1.word_count = decode.parameters1.word_count + 1
    GNSSDecoder.decodeword(word9,decode.data,decode.parameters1,_interm)

    @test decode.parameters1.word_count == 9
    @test decode.parameters1.prev_29 == false
    @test decode.parameters1.prev_30 == true

    @test decode.data.sqrt_A ≈ 5153.662601

end

@testset "word10sub2" begin
    #initialization
    decode = GNSSDecoder.init_decode()
    _interm = GNSSDecoder.GPSData_interm()
    decode.parameters1.prev_30 = false

    #first word
    word1 = [true,false,false,false,true,false,true,true,false,false,false,false,true,true,true,true,true,false,false,false,false,false,false,false,true,false,false,true,true,false]
    word1 = word1[30:-1:1]

    decode.parameters1.word_count = decode.parameters1.word_count + 1
    GNSSDecoder.decodeword(word1,decode.data,decode.parameters1,_interm)

    #second word
    word2 = [true,false,false,true,true,false,true,true,true,true,true,false,true,false,false,true,true,false,true,false,true,false,false,false,true,false,true,false,false,false]
    word2 = word2[30:-1:1]

    decode.parameters1.word_count = decode.parameters1.word_count + 1
    GNSSDecoder.decodeword(word2,decode.data,decode.parameters1,_interm)

    #third word
    word3 = [false,true,true,true,true,false,true,true,false,false,false,false,false,false,true,true,true,true,false,false,true,false,true,true,false,false,true,true,true,false]
    word3 = word3[30:-1:1]

    decode.parameters1.word_count = decode.parameters1.word_count + 1
    GNSSDecoder.decodeword(word3,decode.data,decode.parameters1,_interm)

    #fourth word
    word4 = [false,false,true,true,false,false,true,false,true,true,true,true,true,false,false,false,false,false,true,false,false,true,true,true,true,false,false,false,false,false]
    word4 = word4[30:-1:1]

    decode.parameters1.word_count = decode.parameters1.word_count + 1
    GNSSDecoder.decodeword(word4,decode.data,decode.parameters1,_interm)

    #fifth word
    word5 = [true,true,true,false,false,true,false,false,false,true,true,false,true,true,false,true,true,true,false,false,false,false,false,true,false,false,true,false,true,false]
    word5 = word5[30:-1:1]

    decode.parameters1.word_count = decode.parameters1.word_count + 1
    GNSSDecoder.decodeword(word5,decode.data,decode.parameters1,_interm)

    #sixth word
    word6 = [false,false,false,false,false,false,true,true,true,false,true,true,true,false,false,false,false,false,false,false,false,true,false,false,false,false,false,true,true,false]
    word6 = word6[30:-1:1]

    decode.parameters1.word_count = decode.parameters1.word_count + 1
    GNSSDecoder.decodeword(word6,decode.data,decode.parameters1,_interm)

    #seventh word
    word7 = [false,false,true,false,false,true,false,true,true,true,true,true,true,true,true,true,false,false,false,false,true,false,false,false,false,false,false,false,false,true]
    word7 = word7[30:-1:1]

    decode.parameters1.word_count = decode.parameters1.word_count + 1
    GNSSDecoder.decodeword(word7,decode.data,decode.parameters1,_interm)

    #eighth word
    word8 = [true,true,true,true,false,false,true,false,true,true,true,false,false,false,false,true,false,true,false,true,true,true,true,false,false,false,true,false,true,true]
    word8 = word8[30:-1:1]

    decode.parameters1.word_count = decode.parameters1.word_count + 1
    GNSSDecoder.decodeword(word8,decode.data,decode.parameters1,_interm)

    #ninth word
    word9 = [true,true,true,true,false,false,true,false,true,false,true,true,false,false,true,false,true,true,true,true,true,true,false,true,true,false,false,false,false,true]
    word9 = word9[30:-1:1]

    decode.parameters1.word_count = decode.parameters1.word_count + 1
    GNSSDecoder.decodeword(word9,decode.data,decode.parameters1,_interm)

    #tenth word
    word10 = [true,false,false,false,true,false,true,false,false,false,true,true,true,false,false,true,true,false,false,false,false,false,true,true,true,false,false,false,false,false]
    word10 = word10[30:-1:1]

    decode.parameters1.word_count = decode.parameters1.word_count + 1
    GNSSDecoder.decodeword(word10,decode.data,decode.parameters1,_interm)

    @test decode.parameters1.word_count == 0
    @test decode.parameters1.prev_29 == false
    @test decode.parameters1.prev_30 == false

    @test decode.data.t_oe ≈ 482400

end
