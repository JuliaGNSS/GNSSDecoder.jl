@testset "Navigation-message data frequency" begin
    # `get_data_frequency` on a decoder state forwards to the GNSSSignals data
    # rate of the signal it demodulates, so a decoder state can report its own
    # nav-message symbol rate without the caller re-deriving which signal it is.
    signal_of_state = [
        (GPSL1CADecoderState(1), GPSL1CA),
        (GPSL1C_DDecoderState(1), GPSL1C_D),
        (GPSL5IDecoderState(1), GPSL5I),
        (GPSL2CMDecoderState(1), GPSL2CM),
        (GalileoE1BDecoderState(1), GalileoE1B),
        (GalileoE5aDecoderState(1), GalileoE5aI),
    ]
    for (state, signal) in signal_of_state
        @test get_data_frequency(state) == get_data_frequency(signal)
    end

    # GPS L5-I and L2C-M share `GPSCNAVData` but run at different symbol rates;
    # dispatch on the constants type must keep them distinct.
    @test get_data_frequency(GPSL5IDecoderState(1)) !=
          get_data_frequency(GPSL2CMDecoderState(1))

    # The E1B BOC(1,1) approximation decodes the identical I/NAV stream, so its
    # decoder reports the same rate as full E1B.
    @test get_data_frequency(GNSSDecoder.GNSSDecoderState(GalileoE1B_BOC11(), 1)) ==
          get_data_frequency(GalileoE1B)
end
