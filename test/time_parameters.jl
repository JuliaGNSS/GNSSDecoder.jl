# The broadcast UTC parameters, leap seconds and Galileo GGTO, with their truncated
# reference weeks resolved onto the decoder's own week count. Every week is checked
# on both sides of its field's modulus wrap, and across a week crossover, where the
# decoded week still lags the reference week of the next one.

using GNSSDecoder: GPSCNAV2Data, BeiDouCNAV2Data

@testset "Broadcast time parameters" begin
    utc(state; year = 2026) = get_utc_parameters(state; approximate_year = year)
    leap(state; year = 2026) = get_leap_second_parameters(state; approximate_year = year)

    @testset "GPS LNAV: 8-bit weeks against a 10-bit week" begin
        # GPS week 2304 (March 2024) is 256 as broadcast, and a multiple of 256.
        lnav(WN, WN_0t, WN_LSF = WN_0t) = GNSSDecoderState(
            GPSL1CADecoderState(1);
            data = GPSLNAVData(;
                WN = Int64(WN),
                A_0UTC = 1.0e-9,
                A_1UTC = 2.0e-15,
                t_0t = Int64(405504),
                WN_0t = Int64(WN_0t),
                Δt_LS = Int64(18),
                Δt_LSF = Int64(18),
                WN_LSF = Int64(WN_LSF),
                DN = Int64(7),
            ),
        )
        p = utc(lnav(256, 0); year = 2024)
        @test p == UTCParameters(1.0e-9, 2.0e-15, 0.0, 405504, 2304)
        # Either side of the 8-bit wrap.
        @test utc(lnav(256, 255); year = 2024).WN_0t == 2303
        @test utc(lnav(256, 1); year = 2024).WN_0t == 2305
        @test leap(lnav(256, 0, 255); year = 2024) == LeapSecondParameters(18, 18, 2303, 7)
        @test leap(lnav(256, 0, 1); year = 2024).WN_LSF == 2305
        # Week crossover: the decoded week is still 2303 while the reference week is
        # 2304, which is 0 as broadcast — the next week, not 255 weeks back.
        @test utc(lnav(255, 0); year = 2024).WN_0t == 2304
        # The same across the 1024-week rollover of 2019-04-07: week 2047 is 1023
        # as broadcast, and the reference week 2048 is 0.
        @test utc(lnav(1023, 0); year = 2019).WN_0t == 2048
        @test leap(lnav(1023, 0, 0); year = 2019).WN_LSF == 2048
        # LNAV broadcasts no `A_2`.
        @test p.A_2 == 0.0
    end

    @testset "GPS CNAV and CNAV-2: 13-bit weeks" begin
        fields(WN, WN_0t) = (;
            WN = Int64(WN),
            A_0UTC = 1.0e-9,
            A_1UTC = 2.0e-15,
            A_2UTC = 3.0e-20,
            t_0t = Int64(345600),
            WN_0t = Int64(WN_0t),
            Δt_LS = Int64(18),
            Δt_LSF = Int64(19),
            WN_LSF = Int64(WN_0t),
            DN = Int64(1),
        )
        for (decoder, Data) in (
            (GPSL5IDecoderState(1), GPSCNAVData),
            (GPSL2CMDecoderState(1), GPSCNAVData),
            (GPSL1C_DDecoderState(1), GPSCNAV2Data),
        )
            state(WN, WN_0t) =
                GNSSDecoderState(decoder; data = Data(; fields(WN, WN_0t)...))
            @test utc(state(2348, 2348)) ==
                  UTCParameters(1.0e-9, 2.0e-15, 3.0e-20, 345600, 2348)
            @test leap(state(2348, 2348)) == LeapSecondParameters(18, 19, 2348, 1)
            # The field is as wide as the week: only its own wrap, at week 8192, moves it.
            @test utc(state(8191, 8191)).WN_0t == 8191
            @test utc(state(8191, 0)).WN_0t == 8192
            @test leap(state(8191, 0)).WN_LSF == 8192
        end
    end

    @testset "Galileo I/NAV and F/NAV: 8-bit UTC and 6-bit GGTO weeks" begin
        # GST week 1280 is a multiple of both 256 and 64.
        fields(WN, WN_0t, WN_0G) = (;
            WN = Int64(WN),
            A_0UTC = 4.0e-9,
            A_1UTC = 5.0e-15,
            t_0t = 3600 * 24,
            WN_0t = WN_0t,
            Δt_LS = 18,
            Δt_LSF = 18,
            WN_LSF = WN_0t,
            DN = 7,
            A_0G = 6.0e-9,
            A_1G = 7.0e-15,
            t_0G = 3600 * 48,
            WN_0G = WN_0G,
        )
        for (decoder, Data) in (
            (GalileoE1BDecoderState(1), GalileoINAVData),
            (GalileoE5bDecoderState(1), GalileoINAVData),
            (GalileoE5aDecoderState(1), GalileoFNAVData),
        )
            state(WN, WN_0t, WN_0G) =
                GNSSDecoderState(decoder; data = Data(; fields(WN, WN_0t, WN_0G)...))
            s = state(1280, 0, 0)
            @test utc(s) == UTCParameters(4.0e-9, 5.0e-15, 0.0, 86400, 1280)
            @test leap(s) == LeapSecondParameters(18, 18, 1280, 7)
            @test get_ggto_parameters(s) == GGTOParameters(6.0e-9, 7.0e-15, 172800, 1280)
            # Either side of each wrap.
            @test utc(state(1280, 255, 63)).WN_0t == 1279
            @test leap(state(1280, 255, 63)).WN_LSF == 1279
            @test get_ggto_parameters(state(1280, 255, 63)).WN_0G == 1279
            @test get_ggto_parameters(state(1280, 1, 1)).WN_0G == 1281
            # Week crossover: decoded week 1279, reference week 1280.
            @test utc(state(1279, 0, 0)).WN_0t == 1280
            @test get_ggto_parameters(state(1279, 0, 0)).WN_0G == 1280
            # The GGTO epoch is the one the time-offset accessors evaluate from.
            g = get_ggto_parameters(s)
            @test get_week_start_epoch(GST(), g.WN_0G) + g.t_0G ==
                  get_broadcast_time_offset_epoch(s, GPST())
            # Galileo weeks are full: the year does not enter.
            @test utc(s; year = 1990) == utc(s)
        end
    end

    @testset "BeiDou D1/D2: no UTC reference epoch, an 8-bit leap-second week" begin
        # BDT week 1024 is a multiple of 256.
        dnav(WN, WN_LSF) = GNSSDecoderState(
            BeiDouB1IDecoderState(20);
            data = BeiDouDNAVData(;
                WN = Int64(WN),
                A_0UTC = 8.0e-9,
                A_1UTC = 9.0e-15,
                Δt_LS = Int64(4),
                Δt_LSF = Int64(4),
                WN_LSF = Int64(WN_LSF),
                DN = Int64(6),
            ),
        )
        p = utc(dnav(1024, 0))
        @test p == UTCParameters(8.0e-9, 9.0e-15, 0.0, nothing, nothing)
        @test isnothing(p.t_0t) && isnothing(p.WN_0t)
        @test leap(dnav(1024, 0)) == LeapSecondParameters(4, 4, 1024, 6)
        @test leap(dnav(1024, 255)).WN_LSF == 1023
        @test leap(dnav(1023, 0)).WN_LSF == 1024
        # B3I carries the same message.
        @test utc(GNSSDecoderState(BeiDouB3IDecoderState(20); data = dnav(1024, 0).data)) ==
              p
    end

    @testset "BeiDou B-CNAV: 13-bit weeks" begin
        b2a(WN, WN_0t) = GNSSDecoderState(
            BeiDouB2aDecoderState(20);
            data = BeiDouCNAV2Data(;
                WN = Int64(WN),
                A_0UTC = 1.0e-9,
                A_1UTC = 2.0e-15,
                A_2UTC = 3.0e-20,
                t_0t = Int64(345600),
                WN_0t = Int64(WN_0t),
                Δt_LS = Int64(4),
                Δt_LSF = Int64(5),
                WN_LSF = Int64(WN_0t),
                DN = Int64(0),
            ),
        )
        @test utc(b2a(1024, 1024)) == UTCParameters(1.0e-9, 2.0e-15, 3.0e-20, 345600, 1024)
        @test leap(b2a(1024, 1024)) == LeapSecondParameters(4, 5, 1024, 0)
        @test utc(b2a(8191, 0)).WN_0t == 8192
    end

    @testset "nothing until decoded, and where not broadcast" begin
        for decoder in (
            GPSL1CADecoderState(1),
            GPSL5IDecoderState(1),
            GPSL1C_DDecoderState(1),
            GalileoE1BDecoderState(1),
            GalileoE5aDecoderState(1),
            GalileoE6BDecoderState(1),
            BeiDouB1IDecoderState(20),
            BeiDouB1CDecoderState(20),
            BeiDouB2bDecoderState(20),
        )
            @test isnothing(utc(decoder))
            @test isnothing(leap(decoder))
            @test isnothing(get_ggto_parameters(decoder))
            # No default year, as for `get_broadcast_week`.
            @test_throws UndefKeywordError get_utc_parameters(decoder)
            @test_throws UndefKeywordError get_leap_second_parameters(decoder)
        end
        # The fields alone are not enough: the week to resolve against must be there.
        no_week = GNSSDecoderState(
            GalileoE1BDecoderState(1);
            data = GalileoINAVData(;
                A_0UTC = 4.0e-9,
                A_1UTC = 5.0e-15,
                t_0t = 0,
                WN_0t = 0,
                A_0G = 6.0e-9,
                A_1G = 7.0e-15,
                t_0G = 0,
                WN_0G = 0,
            ),
        )
        @test isnothing(utc(no_week))
        @test isnothing(get_ggto_parameters(no_week))
        # GPS and BeiDou broadcast no Galileo GGTO.
        @test isnothing(
            get_ggto_parameters(
                GNSSDecoderState(
                    GPSL5IDecoderState(1);
                    data = GPSCNAVData(; WN = Int64(2348)),
                ),
            ),
        )
    end

    @testset "allocation-free" begin
        lnav = GNSSDecoderState(
            GPSL1CADecoderState(1);
            data = GPSLNAVData(;
                WN = Int64(256),
                A_0UTC = 1.0e-9,
                A_1UTC = 2.0e-15,
                t_0t = Int64(0),
                WN_0t = Int64(0),
                Δt_LS = Int64(18),
                Δt_LSF = Int64(18),
                WN_LSF = Int64(0),
                DN = Int64(7),
            ),
        )
        dnav = GNSSDecoderState(
            BeiDouB1IDecoderState(20);
            data = BeiDouDNAVData(; WN = Int64(1024), A_0UTC = 8.0e-9, A_1UTC = 9.0e-15),
        )
        galileo = GNSSDecoderState(
            GalileoE1BDecoderState(1);
            data = GalileoINAVData(;
                WN = Int64(1280),
                A_0G = 6.0e-9,
                A_1G = 7.0e-15,
                t_0G = 0,
                WN_0G = 0,
            ),
        )
        # Called directly, not through the keyword closures above, which Julia 1.10
        # does not inline and so would allocate on their own.
        function allocations(state)
            get_utc_parameters(state; approximate_year = 2026)
            get_leap_second_parameters(state; approximate_year = 2026)
            get_ggto_parameters(state)
            (@allocated get_utc_parameters(state; approximate_year = 2026)) +
            (@allocated get_leap_second_parameters(state; approximate_year = 2026)) +
            (@allocated get_ggto_parameters(state))
        end
        @test allocations(lnav) == 0
        @test allocations(dnav) == 0
        @test allocations(galileo) == 0
    end
end
