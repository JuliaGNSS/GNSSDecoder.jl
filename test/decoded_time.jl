@testset "Decoded time at the last consumed symbol" begin
    @test isnothing(get_decoded_time_of_week(GPSL1CADecoderState(1)))
    @test isnothing(get_decoded_time_of_week(GalileoE6BDecoderState(1)))

    gps = complete_state(GPSL1CADecoderState(1), GPSLNAVData; TOW = Int64(604_799))
    # A decoded stamp alone is insufficient after the time anchor was reset.
    @test isnothing(get_decoded_time_of_week(gps))
    gps = GNSSDecoderState(gps; num_bits_since_time_reference = 100)
    @test get_decoded_time_of_week(gps) == 604_801.0
    # The counter advances even when the new symbols contain no complete unit.
    gps = decode!(gps, zeros(Float32, 25))
    @test get_decoded_time_of_week(gps) == 604_801.5

    # Decoding rate, not the ranging pilot's rate; two signals share CNAV data.
    for (state, expected) in
        ((GPSL5IDecoderState(1), 100.5), (GPSL2CMDecoderState(1), 101.0))
        state = complete_state(state, GPSCNAVData; TOW = Int64(100))
        state = GNSSDecoderState(state; num_bits_since_time_reference = 50)
        @test get_decoded_time_of_week(state) == expected
    end

    # The same message format has a different rate on GEO and non-GEO BeiDou.
    for (prn, expected) in ((1, 100.1), (20, 101.0))
        state = complete_state(BeiDouB1IDecoderState(prn), BeiDouDNAVData; SOW = Int64(100))
        state = GNSSDecoderState(state; num_bits_since_time_reference = 50)
        @test get_decoded_time_of_week(state) == expected
    end

    # A vector-loop snapshot can place the last complete bit just before the
    # reference while keeping its fractional code phase on the other side.
    shifted = GNSSDecoderState(gps; num_bits_since_time_reference = -1)
    @test get_decoded_time_of_week(shifted) == 604_798.98

    allocations(state) = @allocated get_decoded_time_of_week(state)
    allocations(gps)
    @test allocations(gps) == 0 skip = !CHECK_ALLOCATIONS
end
