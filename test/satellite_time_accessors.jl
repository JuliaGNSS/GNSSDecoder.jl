using GNSSSignals: GPST, GST, BDT, Hz, get_tai_offset
# `GNSSSignals.s` / `.ustrip` qualified, not imported: this file uses `s` as a
# local for decoder states, which would shadow the second unit.
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
        # `WN` must be present: `WN_0G` is 6 bits and has to be lifted into the
        # 12-bit week numbering before `WN - WN_0` means anything.
        WN_now = Int64(1234)
        for data in (
            GalileoINAVData(; WN = WN_now, scaled...),
            GalileoE5aData(; WN = WN_now, scaled...),
        )
            state = GNSSDecoderState(
                data isa GalileoINAVData ? GalileoE1BDecoderState(1) :
                GalileoE5aDecoderState(1);
                data,
            )
            offset = get_time_offset(state, GPST())
            @test offset.target === GPST()
            # GST and GPST are both TAI - 19 s, so the defined offset is zero
            # and `A_0` is the broadcast coefficient unchanged. This is the case
            # that hides the BeiDou bug, so it is asserted as an equality on
            # purpose rather than left implicit.
            @test GNSSSignals.ustrip(
                GNSSSignals.s,
                get_tai_offset(GPST()) - get_tai_offset(GST()),
            ) == 0
            @test offset.A_0 == scaled.A_0G
            @test offset.A_1 == scaled.A_1G
            @test offset.A_2 == 0.0        # GGTO has no quadratic term
            @test offset.t_0 == scaled.t_0G
            # Resolved, not as broadcast: the raw field is 5, and 5 is not a
            # week number a 12-bit `WN` of 1234 can be differenced against.
            @test scaled.WN_0G == 5
            @test offset.WN_0 != scaled.WN_0G
            @test mod(offset.WN_0, 64) == scaled.WN_0G
            @test abs(WN_now - offset.WN_0) <= 32
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
                    BeiDouB2aData(; WN = Int64(901), GNSS_ID = Int64(GNSS_ID), bgto...),
                ),
                (
                    BeiDouB2bDecoderState(20),
                    BeiDouB2bData(; WN = Int64(901), GNSS_ID = Int64(GNSS_ID), bgto...),
                ),
            )
                s = GNSSDecoderState(state; data)
                offset = get_time_offset(s, target)
                @test offset.target === target
                # BDT is TAI - 33 s and both targets are TAI - 19 s, so the
                # defined part of the offset is -14 s and the broadcast
                # coefficient is only the steering residual on top of it. A
                # 16-bit field at 2^-35 s spans +/- 0.95 us and could not carry
                # 14 s even if the ICD meant it to.
                @test offset.A_0 == bgto.A_0BGTO - 14.0
                @test offset.A_0 ≈ -14.0 atol = 1e-6
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
                data = BeiDouB2aData(; WN = Int64(901), GNSS_ID = Int64(GNSS_ID), bgto...),
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
        s = GNSSDecoderState(
            BeiDouB1CDecoderState(20);
            data = BeiDouB1CData(; WN = Int64(901), bgtos),
        )
        @test get_time_offset(s, GPST()).A_0 == 1.0e-9 - 14.0
        @test get_time_offset(s, GST()).A_0 == 4.0e-9 - 14.0
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
        @test (gps.A_0, gps.A_1, gps.A_2) == (1.5e-9 - 14.0, 2.5e-14, 0.0)
        @test isnothing(gps.t_0) && isnothing(gps.WN_0)
        gal = get_time_offset(s, GST())
        @test (gal.A_0, gal.A_1) == (3.5e-9 - 14.0, 4.5e-14)
        @test isnothing(get_time_offset(s, BDT()))
        # Unlike the tagged shapes, both are answerable at once.
        @test !isnothing(get_time_offset(s, GPST())) &&
              !isnothing(get_time_offset(s, GST()))

        # GPS: the single GGTO set is tagged, and CNAV, CNAV-2 and the L2C
        # carriage of message type 35 all use one table — 0 no data available,
        # 1 Galileo, 2 GLONASS, 3-7 reserved (IS-GPS-705J §20.3.3.8.1,
        # IS-GPS-200N §30.3.3.8.1, IS-GPS-800J §3.5.4.2.1). Galileo is the only
        # target that can be answered: GLONASS has no `TimeSystem`, and the ICDs
        # direct that a reserved code be read as presently unusable rather than
        # guessed at. There is no BeiDou code on any GPS signal.
        ggto = (;
            A_0GGTO = 1.0e-9,
            A_1GGTO = 2.0e-14,
            A_2GGTO = 3.0e-20,
            t_GGTO = Int64(432_000),
            WN_GGTO = Int64(2300),
        )
        for (state, mk) in (
            (
                GPSL5IDecoderState(1),
                id -> GPSCNAVData(; WN = Int64(2300), GNSS_ID = Int64(id), ggto...),
            ),
            (
                GPSL2CMDecoderState(1),
                id -> GPSCNAVData(; WN = Int64(2300), GNSS_ID = Int64(id), ggto...),
            ),
            (
                GPSL1C_DDecoderState(1),
                id -> GPSL1C_DData(; WN = Int64(2300), GGTO_ID = Int64(id), ggto...),
            ),
        )
            galileo = GNSSDecoderState(state; data = mk(1))
            offset = get_time_offset(galileo, GST())
            @test offset.A_0 == 1.0e-9     # GPST and GST are both TAI - 19 s
            @test offset.A_2 == 3.0e-20     # GPS GGTO is the one three-term shape
            @test isnothing(get_time_offset(galileo, GPST()))
            @test isnothing(get_time_offset(galileo, BDT()))
            # 0 no data, 2 GLONASS (no TimeSystem), 3-7 reserved: none answerable.
            for id in (0, 2, 3, 4, 7)
                unusable = GNSSDecoderState(state; data = mk(id))
                for target in (GPST(), GST(), BDT())
                    @test isnothing(get_time_offset(unusable, target))
                end
            end
        end

        @testset "truncated reference week is resolved" begin
            resolve = GNSSDecoder.resolve_reference_week
            # Galileo: 6-bit field against a 12-bit week. Same week, one back,
            # one ahead — all must land next to `WN`, never 63 weeks away.
            @test resolve(mod(1234, 64), 1234, 64) == 1234
            @test resolve(mod(1233, 64), 1234, 64) == 1233
            @test resolve(mod(1235, 64), 1234, 64) == 1235
            # Across the 6-bit field's own wrap, in both directions.
            @test resolve(mod(1216, 64), 1217, 64) == 1216
            @test resolve(mod(1215, 64), 1216, 64) == 1215
            # Every representative resolves back to the field it came from and
            # stays inside half a modulus of the current week.
            for WN = 1200:1300, offset = -31:31
                raw = mod(WN + offset, 64)
                r = resolve(raw, WN, 64)
                @test mod(r, 64) == raw
                @test abs(r - WN) <= 32
            end
            # 13-bit fields: unchanged in the ordinary case, and rescued at a
            # week-number rollover, which is the only time they are wrong.
            @test resolve(900, 901, 8192) == 900
            @test resolve(8190, 2, 8192) == -2      # 4 weeks back, not 8188 forward
            @test 2 - resolve(8190, 2, 8192) == 4
            # Unknown either side is unresolvable rather than guessed.
            @test isnothing(resolve(nothing, 1234, 64))
            @test isnothing(resolve(5, nothing, 64))
        end

        @testset "an unresolvable reference week yields no offset" begin
            # A Galileo satellite can broadcast word type 10 (GGTO) before word
            # type 0 or 5 (WN). Until the week arrives the reference epoch
            # cannot be placed, and `Δτ` could not be formed either, so there is
            # no usable offset rather than one with a 6-bit week in it.
            scaled = GNSSDecoder.galileo_ggto(0x1234, 0x0abc, 0x10, 0x05)
            without = GNSSDecoderState(
                GalileoE1BDecoderState(1);
                data = GalileoINAVData(; scaled...),
            )
            @test isnothing(without.data.WN)
            @test isnothing(get_time_offset(without, GPST()))
            with = GNSSDecoderState(
                GalileoE1BDecoderState(1);
                data = GalileoINAVData(; WN = Int64(1234), scaled...),
            )
            @test !isnothing(get_time_offset(with, GPST()))
            # BeiDou D1/D2 broadcasts no reference epoch at all, so it is
            # unaffected by the week being unknown.
            dnav = GNSSDecoderState(
                BeiDouB1IDecoderState(20);
                data = BeiDouDNAVData(; A_0GPS = 1.5e-9, A_1GPS = 2.5e-14),
            )
            offset = get_time_offset(dnav, GPST())
            @test !isnothing(offset)
            @test isnothing(offset.WN_0) && isnothing(offset.t_0)
        end

        # The contract itself, on every decoder that answers, rather than the
        # arithmetic of any one of them: `t_target = t_own - Δt` must hold for
        # the *defined* part of the offset, which is the part no message
        # carries. Broadcasting zero coefficients isolates it.
        @testset "defined scale offset is included" begin
            zeroed = (;
                A_0GGTO = 0.0,
                A_1GGTO = 0.0,
                A_2GGTO = 0.0,
                t_GGTO = Int64(0),
                WN_GGTO = Int64(0),
            )
            zbgto = (;
                GNSS_ID = Int64(1),
                WN_0BGTO = Int64(0),
                t_0BGTO = Int64(0),
                A_0BGTO = 0.0,
                A_1BGTO = 0.0,
                A_2BGTO = 0.0,
            )
            zggto = GNSSDecoder.galileo_ggto(0x0000, 0x0000, 0x00, 0x00)
            cases = [
                (
                    GNSSDecoderState(
                        GPSL5IDecoderState(1);
                        data = GPSCNAVData(;
                            WN = Int64(2300),
                            GNSS_ID = Int64(1),
                            zeroed...,
                        ),
                    ),
                    GST(),
                ),
                (
                    GNSSDecoderState(
                        GPSL1C_DDecoderState(1);
                        data = GPSL1C_DData(;
                            WN = Int64(2300),
                            GGTO_ID = Int64(1),
                            zeroed...,
                        ),
                    ),
                    GST(),
                ),
                (
                    GNSSDecoderState(
                        GalileoE1BDecoderState(1);
                        data = GalileoINAVData(; WN = Int64(1234), zggto...),
                    ),
                    GPST(),
                ),
                (
                    GNSSDecoderState(
                        GalileoE5aDecoderState(1);
                        data = GalileoE5aData(; WN = Int64(1234), zggto...),
                    ),
                    GPST(),
                ),
                (
                    GNSSDecoderState(
                        BeiDouB2aDecoderState(20);
                        data = BeiDouB2aData(; WN = Int64(901), zbgto...),
                    ),
                    GPST(),
                ),
                (
                    GNSSDecoderState(
                        BeiDouB2bDecoderState(20);
                        data = BeiDouB2bData(; WN = Int64(901), zbgto...),
                    ),
                    GPST(),
                ),
                (
                    GNSSDecoderState(
                        BeiDouB1CDecoderState(20);
                        data = BeiDouB1CData(;
                            WN = Int64(901),
                            bgtos = Dictionary([1], [BeiDouB1CBGTO(; zbgto...)]),
                        ),
                    ),
                    GPST(),
                ),
                (
                    GNSSDecoderState(
                        BeiDouB1IDecoderState(20);
                        data = BeiDouDNAVData(; A_0GPS = 0.0, A_1GPS = 0.0),
                    ),
                    GPST(),
                ),
            ]
            for (state, target) in cases
                offset = get_time_offset(state, target)
                own = get_time_system(state)
                # A zero broadcast bias must still leave the defined offset,
                # which is what `t_target = t_own - Δt` needs to be true.
                expected = GNSSSignals.ustrip(
                    GNSSSignals.s,
                    get_tai_offset(target) - get_tai_offset(own),
                )
                @test offset.A_0 == expected
                # And that is 0 for every GPS/Galileo pair but -14 s for BeiDou,
                # so a suite that only exercised Galileo would pass either way.
                @test expected == (own === BDT() ? -14.0 : 0.0)
            end
        end

        # GPS L1 C/A carries UTC parameters but no inter-GNSS offset.
        for target in (GPST(), GST(), BDT())
            @test isnothing(get_time_offset(GPSL1CADecoderState(1), target))
        end
    end
end
