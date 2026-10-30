@testset "decode! symbol count API" begin
    messages = load_l2c_fixture_messages()
    stream = fec_encode_soft(CNAVTestEncoder(), reduce(vcat, messages))

    @testset "default count consumes the full input" begin
        defaulted = decode!(GPSL2CMDecoderState(9), stream)
        explicit = decode!(GPSL2CMDecoderState(9), stream, length(stream))
        @test defaulted == explicit
        @test GNSSDecoder.has_positioning_data(defaulted)
    end

    @testset "explicit count consumes only the leading symbols" begin
        prefix_length = length(stream) - 100
        prefix = view(stream, 1:prefix_length)
        with_tail = vcat(stream, fill(0.0f0, 100))
        defaulted = decode!(GPSL2CMDecoderState(9), prefix)
        explicit = decode!(GPSL2CMDecoderState(9), with_tail, prefix_length)
        @test explicit == defaulted
    end

    @testset "empty input and count validation" begin
        state = GPSL2CMDecoderState(9)
        @test decode!(state, Float32[]) == decode!(state, Float32[], 0)
        @test_throws ArgumentError decode!(state, Float32[], -1)
        @test_throws ArgumentError decode!(state, Float32[], 1)
    end
end
