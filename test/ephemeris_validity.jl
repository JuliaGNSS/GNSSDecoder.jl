@testset "Broadcast ephemeris validity metadata" begin
    @testset "week placement uses the broadcast stamp" begin
        data = complete_data(
            GPSLNAVData;
            WN = Int64(0),
            TOW = Int64(10),
            t_0e = Int64(604_784),
            t_0c = Int64(12),
        )
        state = GNSSDecoderState(GPSL1CADecoderState(1); data)
        validity = get_broadcast_ephemeris_validity(state; approximate_year = 2020)
        week = get_broadcast_week(state; approximate_year = 2020)
        week_start = get_week_start_epoch(state, week)

        @test validity isa BroadcastEphemerisValidity
        @test validity.ephemeris_reference == week_start - 16
        @test validity.clock_reference == week_start + 12
        @test validity.ephemeris_reference isa GNSSSignals.TAITime
        @test validity.clock_reference isa GNSSSignals.TAITime
    end

    @testset "references can fall in the following week" begin
        data = complete_data(
            GPSLNAVData;
            WN = Int64(0),
            TOW = Int64(604_790),
            t_0e = Int64(8),
            t_0c = Int64(604_780),
        )
        state = GNSSDecoderState(GPSL1CADecoderState(1); data)
        validity = get_broadcast_ephemeris_validity(state; approximate_year = 2020)
        week = get_broadcast_week(state; approximate_year = 2020)
        week_start = get_week_start_epoch(state, week)

        @test validity.ephemeris_reference ==
              week_start + GNSSSignals.SECONDS_PER_WEEK + 8
        @test validity.clock_reference == week_start + 604_780
    end

    @testset "incomplete data and ephemeris-free messages" begin
        for state in (
            GPSL1CADecoderState(1),
            GalileoE6BDecoderState(1),
        )
            @test isnothing(
                get_broadcast_ephemeris_validity(state; approximate_year = 2020),
            )
        end

        missing_stamp = complete_state(
            GPSL1CADecoderState(1), GPSLNAVData;
            TOW = nothing,
        )
        @test isnothing(
            get_broadcast_ephemeris_validity(missing_stamp; approximate_year = 2020),
        )

        missing_week = complete_state(
            GPSL1CADecoderState(1), GPSLNAVData;
            WN = nothing,
        )
        @test isnothing(
            get_broadcast_ephemeris_validity(missing_week; approximate_year = 2020),
        )
    end

    @testset "LNAV fit duration follows flag and IODC" begin
        function fit_seconds(flag, iodc)
            state = complete_state(
                GPSL1CADecoderState(1), GPSLNAVData;
                fit_interval = flag,
                IODC = Int64(iodc),
            )
            validity = get_broadcast_ephemeris_validity(
                state;
                approximate_year = 2020,
            )
            validity.fit_interval_seconds
        end

        @test fit_seconds(false, 497) == 4 * 3600
        @test fit_seconds(true, 12) == 6 * 3600
        @test fit_seconds(true, 300) == 6 * 3600
        @test fit_seconds(true, 244) == 8 * 3600
        @test fit_seconds(true, 496) == 14 * 3600
        @test fit_seconds(true, 500) == 26 * 3600
        @test fit_seconds(true, 1022) == 26 * 3600
        for reserved_iodc in (504, 511, 752, 767, 1008, 1020)
            @test isnothing(get_broadcast_ephemeris_validity(
                complete_state(
                    GPSL1CADecoderState(1), GPSLNAVData;
                    fit_interval = true,
                    IODC = Int64(reserved_iodc),
                );
                approximate_year = 2020,
            ).fit_interval_seconds)
        end

        galileo = complete_state(GalileoE1BDecoderState(1), GalileoINAVData)
        validity = get_broadcast_ephemeris_validity(galileo; approximate_year = 2020)
        @test isnothing(validity.fit_interval_seconds)

        function allocations(state)
            get_broadcast_ephemeris_validity(state; approximate_year = 2020)
            @allocated get_broadcast_ephemeris_validity(
                state;
                approximate_year = 2020,
            )
        end
        complete = complete_state(GPSL1CADecoderState(1), GPSLNAVData)
        @test allocations(complete) == 0
        @test allocations(galileo) == 0
    end
end
