@testset "Signal accessors" begin
    # Every decoder state maps to exactly one GNSSSignals signal type, and each
    # forwarded accessor must agree with what GNSSSignals reports for that
    # signal. The pairs below are the whole mapping; the data-bearing component
    # is the one named, since that is what carries the navigation message
    # (`GPSL2CM` not `GPSL2CL`, `GalileoE5aI` not `GalileoE5aQ`, `GPSL1C_D` not
    # `GPSL1C_P`).
    signal_of_state = [
        (GPSL1CADecoderState(1), GPSL1CA),
        (GPSL1C_DDecoderState(1), GPSL1C_D),
        (GPSL5IDecoderState(1), GPSL5I),
        (GPSL2CMDecoderState(1), GPSL2CM),
        (GalileoE1BDecoderState(1), GalileoE1B),
        (GalileoE5aDecoderState(1), GalileoE5aI),
    ]

    # The accessors GNSSDecoder forwards from a decoder state to its signal.
    forwarded_accessors = [
        get_signal_id,
        get_signal_name,
        get_constellation_id,
        get_constellation_name,
        get_band,
        get_band_id,
        get_band_name,
        get_data_frequency,
        get_time_system,
        get_time_system_id,
        get_time_system_name,
        get_system_start_time,
        get_tai_offset,
    ]

    @testset "forwarded to the signal" begin
        for (state, signal) in signal_of_state
            @test get_signal_type(state) === signal
            for accessor in forwarded_accessors
                @test accessor(state) == accessor(signal)
            end
        end
    end

    @testset "identity" begin
        # A decoder state can name itself without the caller carrying the signal
        # alongside it — the ids for keying, the names for log output.
        state_of_id = Dict(get_signal_id(state) => state for (state, _) in signal_of_state)
        @test length(state_of_id) == length(signal_of_state)   # ids are distinct
        @test get_signal_id(state_of_id[:GPSL1CA]) == :GPSL1CA
        @test get_signal_name(state_of_id[:GPSL1CA]) == "GPS L1 C/A"

        for (state, _) in signal_of_state
            # Constellation membership agrees with the data supertype this
            # package sorts its decoded data into (see `gnss_supertype.jl`).
            expected = state.data isa GNSSDecoder.AbstractGPSData ? :GPS : :Galileo
            @test get_constellation_id(state) == expected
            @test get_constellation_name(state) == String(expected)
            @test get_band_name(state) == String(get_band_id(state))
        end

        @test get_band_id(GPSL1CADecoderState(1)) == :L1
        @test get_band_id(GPSL2CMDecoderState(1)) == :L2
        @test get_band_id(GalileoE5aDecoderState(1)) == :L5
        # Bands are identified by RF frequency, not by ICD label: Galileo E1B
        # shares GPS L1's carrier and so reports `:L1`, not `:E1`.
        @test get_band_id(GalileoE1BDecoderState(1)) == :L1
    end

    @testset "GPS L5-I vs L2C-M" begin
        # The two share `GPSCNAVData` and every line of decoding; only the
        # constants type tells them apart, so the signal mapping — and with it
        # every forwarded accessor — must stay keyed on the constants.
        l5 = GPSL5IDecoderState(1)
        l2c = GPSL2CMDecoderState(1)
        @test get_signal_type(l5) !== get_signal_type(l2c)
        @test get_data_frequency(l5) != get_data_frequency(l2c)
        @test get_band_id(l5) == :L5
        @test get_band_id(l2c) == :L2
        # Same constellation, same time scale, though.
        @test get_constellation_id(l5) == get_constellation_id(l2c) == :GPS
        @test get_time_system(l5) == get_time_system(l2c)
    end

    @testset "time system" begin
        # The decoded week numbers and times of week are counted in the
        # constellation's time scale; a consumer turning them into an absolute
        # instant needs its epoch and TAI offset, and gets both from the state.
        for (state, signal) in signal_of_state
            is_gps = state.data isa GNSSDecoder.AbstractGPSData
            @test get_time_system(state) == (is_gps ? GPST() : GST())
            @test get_time_system_id(state) == (is_gps ? :GPST : :GST)
            @test get_system_start_time(state) == get_system_start_time(signal)
        end
        @test get_time_system_name(GPSL1CADecoderState(1)) == "GPS Time"
        @test get_time_system_name(GalileoE1BDecoderState(1)) == "Galileo System Time"
        # The Galileo epoch is 13 s shy of a UTC midnight, so it is worth
        # sourcing rather than restating — the two epochs genuinely differ.
        @test get_system_start_time(GPSL1CADecoderState(1)) !=
              get_system_start_time(GalileoE1BDecoderState(1))
        # Both currently modelled scales run 19 s behind TAI.
        @test get_tai_offset(GPSL1CADecoderState(1)) ==
              get_tai_offset(GalileoE1BDecoderState(1))
    end

    @testset "E1B BOC(1,1) approximation" begin
        # The BOC(1,1) approximation decodes the identical I/NAV stream, so its
        # decoder *is* an E1B decoder and reports E1B throughout — the
        # approximation is a tracking/acquisition concern and leaves no trace in
        # the symbol domain.
        boc11 = GNSSDecoder.GNSSDecoderState(GalileoE1B_BOC11(), 1)
        @test get_signal_type(boc11) === GalileoE1B
        for accessor in forwarded_accessors
            @test accessor(boc11) == accessor(GalileoE1B)
        end
    end

    @testset "inferred from the type" begin
        # The mapping is keyed on the constants type, so the signal — and every
        # accessor forwarded through it — is known to the compiler and costs
        # nothing at run time.
        state = GPSL1CADecoderState(1)
        @test @inferred(get_signal_type(state)) === GPSL1CA
        @test @inferred(get_signal_id(state)) == :GPSL1CA
        @test @inferred(get_band_id(state)) == :L1
        @test @inferred(get_data_frequency(state)) == get_data_frequency(GPSL1CA)
        @test @inferred(get_time_system(state)) == GPST()
    end
end
