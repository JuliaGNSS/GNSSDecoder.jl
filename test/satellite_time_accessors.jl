using GNSSSignals: GPST, GST, BDT, Hz
using GNSSDecoder: geostationary_orbit, inclined_geosynchronous_orbit, medium_earth_orbit

@testset "Satellite and time accessors" begin
    # Every decoder state this package can build. The three accessors under test
    # are cross-signal, so "every signal answers" is part of the contract: a new
    # decoder that forgets to state its orbit class or time of week must fail
    # here rather than throw a MethodError at a consumer.
    all_states() = [
        GPSL1CADecoderState(1),
        GPSL1C_DDecoderState(1),
        GPSL5IDecoderState(1),
        GPSL2CMDecoderState(1),
        GalileoE1BDecoderState(1),
        GalileoE5bDecoderState(1),
        GalileoE5aDecoderState(1),
        GalileoE6BDecoderState(1),
        BeiDouB1IDecoderState(20),
        BeiDouB3IDecoderState(20),
        BeiDouB1CDecoderState(20),
        BeiDouB2aDecoderState(20),
        BeiDouB2bDecoderState(20),
    ]

    @testset "every signal answers all three" begin
        for state in all_states()
            @test get_orbit_class(state) isa Union{Nothing,GNSSDecoder.OrbitClass}
            @test isnothing(get_time_of_week(state))   # nothing decoded yet
            for target in (GPST(), GST(), BDT())
                @test isnothing(get_time_offset(state, target))
            end
        end
    end

    @testset "orbit class" begin
        # GPS and Galileo are MEO-only constellations, so the answer needs no
        # decoded data.
        for state in (
            GPSL1CADecoderState(1),
            GPSL1C_DDecoderState(1),
            GPSL5IDecoderState(1),
            GPSL2CMDecoderState(1),
            GalileoE1BDecoderState(1),
            GalileoE5aDecoderState(1),
            GalileoE6BDecoderState(1),
        )
            @test get_orbit_class(state) === medium_earth_orbit
        end

        # B1I/B3I: the ICD's GEO PRN partition, and nothing finer. This is the
        # fact `PositionVelocityTime.jl` used to infer from the 500 bps D2
        # symbol rate — the two must agree, since the rate is chosen from the
        # same partition.
        for prn in (1, 5, 59, 63)
            state = BeiDouB1IDecoderState(prn)
            @test get_orbit_class(state) === geostationary_orbit
            @test get_data_frequency(state) == 500Hz
        end
        for prn in (6, 20, 58)
            state = BeiDouB3IDecoderState(prn)
            @test isnothing(get_orbit_class(state))   # IGSO and MEO are not distinguished
            @test get_data_frequency(state) == 50Hz
        end

        # B-CNAV: from the satellite's own broadcast `sat_type`, reserved code 0
        # and an undecoded field both reading `nothing`.
        expected = Dict(
            0 => nothing,
            1 => geostationary_orbit,
            2 => inclined_geosynchronous_orbit,
            3 => medium_earth_orbit,
        )
        for (sat_type, class) in expected
            b1c = GNSSDecoderState(
                BeiDouB1CDecoderState(20);
                data = BeiDouB1CData(; sat_type = Int64(sat_type)),
            )
            b2a = GNSSDecoderState(
                BeiDouB2aDecoderState(20);
                data = BeiDouB2aData(; sat_type = Int64(sat_type)),
            )
            b2b = GNSSDecoderState(
                BeiDouB2bDecoderState(20);
                data = BeiDouB2bData(; sat_type = Int64(sat_type)),
            )
            @test get_orbit_class(b1c) === class
            @test get_orbit_class(b2a) === class
            @test get_orbit_class(b2b) === class
        end
    end

    @testset "time of week" begin
        # The signals that broadcast seconds outright pass them through.
        @test get_time_of_week(GPSL1CAData(; TOW = Int64(432_006))) == 432_006
        @test get_time_of_week(GalileoINAVData(; TOW = Int64(123_456))) == 123_456
        @test get_time_of_week(GalileoE5aData(; TOW = Int64(123_456))) == 123_456
        @test get_time_of_week(BeiDouDNAVData(; SOW = Int64(99))) == 99
        @test get_time_of_week(BeiDouB2aData(; SOW = Int64(300_003))) == 300_003
        @test get_time_of_week(BeiDouB2bData(; SOW = Int64(300_001))) == 300_001

        # GPS L1C-D: ITOW counts two-hour intervals, toi 18-second frames.
        @test get_time_of_week(GPSL1C_DData(; ITOW = Int64(0), toi = 0)) == 0
        @test get_time_of_week(GPSL1C_DData(; ITOW = Int64(3), toi = 7)) ==
              3 * 7200 + 7 * 18
        # The last frame of the week: interval 83, frame 399.
        @test get_time_of_week(GPSL1C_DData(; ITOW = Int64(83), toi = 399)) == 604_800 - 18
        # Either field alone is not a time of week.
        @test isnothing(get_time_of_week(GPSL1C_DData(; ITOW = Int64(3))))
        @test isnothing(get_time_of_week(GPSL1C_DData(; toi = 7)))

        # BeiDou B1C: HOW counts hours, soh 18-second frames.
        @test get_time_of_week(BeiDouB1CData(; HOW = Int64(0), soh = 0)) == 0
        @test get_time_of_week(BeiDouB1CData(; HOW = Int64(5), soh = 9)) ==
              5 * 3600 + 9 * 18
        # The last frame of the week: hour 167, frame 199.
        @test get_time_of_week(BeiDouB1CData(; HOW = Int64(167), soh = 199)) == 604_800 - 18
        @test isnothing(get_time_of_week(BeiDouB1CData(; HOW = Int64(5))))
        @test isnothing(get_time_of_week(BeiDouB1CData(; soh = 9)))

        # E6-B stamps HAS blocks with a time of hour, not a time of week.
        @test isnothing(get_time_of_week(GalileoE6BData()))
    end

    @testset "inter-system time offset" begin
        # Galileo: GGTO is GPS-only, two-term, and absent as a set of all-ones.
        scaled = GNSSDecoder.galileo_ggto(0x1234, 0x0abc, 0x10, 0x05)
        for data in (GalileoINAVData(; scaled...), GalileoE5aData(; scaled...))
            state = GNSSDecoderState(
                data isa GalileoINAVData ? GalileoE1BDecoderState(1) :
                GalileoE5aDecoderState(1);
                data,
            )
            offset = get_time_offset(state, GPST())
            @test offset.target === GPST()
            @test offset.A_0 == scaled.A_0G
            @test offset.A_1 == scaled.A_1G
            @test offset.A_2 == 0.0        # GGTO has no quadratic term
            @test offset.t_0 == scaled.t_0G
            @test offset.WN_0 == scaled.WN_0G
            # Galileo broadcasts no offset to any other scale.
            @test isnothing(get_time_offset(state, GST()))
            @test isnothing(get_time_offset(state, BDT()))
        end
        # All four fields all ones: "GGTO not valid" (OS SIS ICD 5.1.8).
        absent = GNSSDecoder.galileo_ggto(0xffff, 0x0fff, 0xff, 0x3f)
        @test isnothing(
            get_time_offset(
                GNSSDecoderState(
                    GalileoE1BDecoderState(1);
                    data = GalileoINAVData(; absent...),
                ),
                GPST(),
            ),
        )

        # BeiDou B2a/B2b: one flat set at a time, tagged by GNSS_ID.
        bgto = (;
            WN_0BGTO = Int64(900),
            t_0BGTO = Int64(432_000),
            A_0BGTO = 1.0e-9,
            A_1BGTO = 2.0e-14,
            A_2BGTO = 3.0e-20,
        )
        for (GNSS_ID, target) in ((1, GPST()), (2, GST()))
            for (state, data) in (
                (
                    BeiDouB2aDecoderState(20),
                    BeiDouB2aData(; GNSS_ID = Int64(GNSS_ID), bgto...),
                ),
                (
                    BeiDouB2bDecoderState(20),
                    BeiDouB2bData(; GNSS_ID = Int64(GNSS_ID), bgto...),
                ),
            )
                s = GNSSDecoderState(state; data)
                offset = get_time_offset(s, target)
                @test offset.target === target
                @test offset.A_0 == bgto.A_0BGTO
                @test offset.A_1 == bgto.A_1BGTO
                @test offset.A_2 == bgto.A_2BGTO
                @test offset.t_0 == bgto.t_0BGTO
                @test offset.WN_0 == bgto.WN_0BGTO
                # Only the tagged system is answerable.
                for other in (GPST(), GST(), BDT())
                    other === target || @test isnothing(get_time_offset(s, other))
                end
            end
        end
        # GNSS_ID 0 means "not available" and must not be published as an offset,
        # and GNSS_ID 3 (GLONASS) has no TimeSystem to be asked for.
        for GNSS_ID in (0, 3)
            s = GNSSDecoderState(
                BeiDouB2aDecoderState(20);
                data = BeiDouB2aData(; GNSS_ID = Int64(GNSS_ID), bgto...),
            )
            for target in (GPST(), GST(), BDT())
                @test isnothing(get_time_offset(s, target))
            end
        end

        # BeiDou B1C: keyed by GNSS ID, so more than one system stays answerable.
        bgtos = Dictionary(
            [1, 2],
            [
                BeiDouB1CBGTO(;
                    GNSS_ID = 1,
                    WN_0BGTO = 900,
                    t_0BGTO = 432_000,
                    A_0BGTO = 1.0e-9,
                    A_1BGTO = 2.0e-14,
                    A_2BGTO = 3.0e-20,
                ),
                BeiDouB1CBGTO(;
                    GNSS_ID = 2,
                    WN_0BGTO = 901,
                    t_0BGTO = 86_400,
                    A_0BGTO = 4.0e-9,
                    A_1BGTO = 5.0e-14,
                    A_2BGTO = 6.0e-20,
                ),
            ],
        )
        s = GNSSDecoderState(BeiDouB1CDecoderState(20); data = BeiDouB1CData(; bgtos))
        @test get_time_offset(s, GPST()).A_0 == 1.0e-9
        @test get_time_offset(s, GST()).A_0 == 4.0e-9
        @test get_time_offset(s, GST()).WN_0 == 901
        @test isnothing(get_time_offset(s, BDT()))

        # BeiDou D1/D2: named per system, two-term, and with no reference epoch.
        s = GNSSDecoderState(
            BeiDouB1IDecoderState(20);
            data = BeiDouDNAVData(;
                A_0GPS = 1.5e-9,
                A_1GPS = 2.5e-14,
                A_0Gal = 3.5e-9,
                A_1Gal = 4.5e-14,
            ),
        )
        gps = get_time_offset(s, GPST())
        @test (gps.A_0, gps.A_1, gps.A_2) == (1.5e-9, 2.5e-14, 0.0)
        @test isnothing(gps.t_0) && isnothing(gps.WN_0)
        gal = get_time_offset(s, GST())
        @test (gal.A_0, gal.A_1) == (3.5e-9, 4.5e-14)
        @test isnothing(get_time_offset(s, BDT()))
        # Unlike the tagged shapes, both are answerable at once.
        @test !isnothing(get_time_offset(s, GPST())) &&
              !isnothing(get_time_offset(s, GST()))

        # GPS: the single GGTO set is tagged, and the two ICDs assign code 3
        # differently — BeiDou on CNAV, reserved on CNAV-2.
        ggto = (;
            A_0GGTO = 1.0e-9,
            A_1GGTO = 2.0e-14,
            A_2GGTO = 3.0e-20,
            t_GGTO = Int64(432_000),
            WN_GGTO = Int64(2300),
        )
        cnav = GNSSDecoderState(
            GPSL5IDecoderState(1);
            data = GPSCNAVData(; GNSS_ID = Int64(3), ggto...),
        )
        @test get_time_offset(cnav, BDT()).A_2 == 3.0e-20
        @test isnothing(get_time_offset(cnav, GST()))
        l1cd = GNSSDecoderState(
            GPSL1C_DDecoderState(1);
            data = GPSL1C_DData(; GGTO_ID = Int64(3), ggto...),
        )
        for target in (GPST(), GST(), BDT())
            @test isnothing(get_time_offset(l1cd, target))   # code 3 is reserved here
        end
        l1cd_gal = GNSSDecoderState(
            GPSL1C_DDecoderState(1);
            data = GPSL1C_DData(; GGTO_ID = Int64(1), ggto...),
        )
        @test get_time_offset(l1cd_gal, GST()).A_0 == 1.0e-9

        # GPS L1 C/A carries UTC parameters but no inter-GNSS offset.
        for target in (GPST(), GST(), BDT())
            @test isnothing(get_time_offset(GPSL1CADecoderState(1), target))
        end
    end
end
