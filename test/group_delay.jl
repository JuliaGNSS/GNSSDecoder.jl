using GNSSSignals: get_center_frequency

@testset "Group delay" begin
    T_GD, ISC, ISC_pilot = 1.0e-8, 2.0e-9, 3.0e-9

    @testset "nothing until the positioning set is complete" begin
        for (state, signal) in (
            (GPSL1CADecoderState(1), GPSL1CA()),
            (GPSL5IDecoderState(1), GPSL5I()),
            (GPSL2CMDecoderState(1), GPSL2CM()),
            (GPSL1C_DDecoderState(1), GPSL1C_D()),
            (GalileoE1BDecoderState(1), GalileoE1B()),
            (GalileoE5bDecoderState(1), GalileoE5bI()),
            (GalileoE5aDecoderState(1), GalileoE5aI()),
            (BeiDouB1IDecoderState(20), BeiDouB1I()),
            (BeiDouB3IDecoderState(20), BeiDouB3I()),
            (BeiDouB1CDecoderState(20), BeiDouB1C_D()),
            (BeiDouB2aDecoderState(20), BeiDouB2aI()),
            (BeiDouB2bDecoderState(20), BeiDouB2bI()),
        )
            @test isnothing(get_group_delay(state, signal))
        end
    end

    @testset "GPS" begin
        l1ca = complete_state(GPSL1CADecoderState(1), GPSL1CAData; T_GD)
        @test is_decoding_completed_for_positioning(l1ca)
        @test get_group_delay(l1ca, GPSL1CA()) == T_GD

        # CNAV: T_GD less the ranging signal's ISC; the pilot has its own, and
        # L2 CL shares L2 CM's.
        cnav = complete_state(
            GPSL5IDecoderState(1),
            GPSCNAVData;
            T_GD,
            ISC_L5I5 = ISC,
            ISC_L5Q5 = ISC_pilot,
            ISC_L2C = 4.0e-9,
        )
        @test get_group_delay(cnav, GPSL5I()) == T_GD - ISC
        @test get_group_delay(cnav, GPSL5Q()) == T_GD - ISC_pilot
        l2c = GNSSDecoderState(GPSL2CMDecoderState(1); data = cnav.data)
        @test get_group_delay(l2c, GPSL2CM()) ==
              get_group_delay(l2c, GPSL2CL()) ==
              T_GD - 4.0e-9
        # Neither term is gated: the fix may go ahead without them, but there is
        # no group delay until both are decoded, and the caller decides what to
        # use meanwhile.
        partial = complete_state(
            GPSL5IDecoderState(1),
            GPSCNAVData;
            T_GD,
            ISC_L5I5 = nothing,
            ISC_L5Q5 = nothing,
        )
        @test is_decoding_completed_for_positioning(partial)
        @test isnothing(get_group_delay(partial, GPSL5I()))
        @test something(get_group_delay(partial, GPSL5I()), 0.0) == 0.0
        no_T_GD = complete_state(GPSL5IDecoderState(1), GPSCNAVData; T_GD = nothing)
        @test is_decoding_completed_for_positioning(no_T_GD)
        @test isnothing(get_group_delay(no_T_GD, GPSL5I()))

        # CNAV-2 carries the L1 set: L1C-D, L1C-P and C/A.
        cnav2 = complete_state(
            GPSL1C_DDecoderState(1),
            GPSL1C_DData;
            T_GD,
            ISC_L1CD = ISC,
            ISC_L1CP = ISC_pilot,
            ISC_L1CA = 4.0e-9,
        )
        @test get_group_delay(cnav2, GPSL1C_D()) == T_GD - ISC
        @test get_group_delay(cnav2, GPSL1C_P()) == T_GD - ISC_pilot
        @test get_group_delay(cnav2, GPSL1CA()) == T_GD - 4.0e-9
        @test isnothing(
            get_group_delay(
                complete_state(
                    GPSL1C_DDecoderState(1),
                    GPSL1C_DData;
                    T_GD,
                    ISC_L1CD = nothing,
                ),
                GPSL1C_D(),
            ),
        )
    end

    @testset "Galileo" begin
        BGD_E1_E5a, BGD_E1_E5b = 2.0e-9, 3.0e-9
        ratio(signal) =
            (get_center_frequency(GalileoE1B()) / get_center_frequency(signal))^2
        # I/NAV: E1 unscaled, E5b scaled by (f_E1/f_E5b)² ≈ 1.70.
        inav = complete_state(
            GalileoE1BDecoderState(1),
            GalileoINAVData;
            BGD_E1_E5a,
            BGD_E1_E5b,
        )
        for signal in (GalileoE1B(), GalileoE1C(), GalileoE1B_BOC11(), GalileoE1C_BOC11())
            @test get_group_delay(inav, signal) == BGD_E1_E5b
        end
        e5b = GNSSDecoderState(GalileoE5bDecoderState(1); data = inav.data)
        for signal in (GalileoE5bI(), GalileoE5bQ())
            @test get_group_delay(e5b, signal) == ratio(GalileoE5bI()) * BGD_E1_E5b
            @test get_group_delay(e5b, signal) ≈ 1.70 * BGD_E1_E5b rtol = 1e-2
        end
        # An E5a range carries an F/NAV decoder, so I/NAV refuses it.
        @test_throws MethodError get_group_delay(inav, GalileoE5aI())

        # F/NAV: E5a scaled by (f_E1/f_E5a)² ≈ 1.79.
        fnav = complete_state(GalileoE5aDecoderState(1), GalileoE5aData; BGD_E1_E5a)
        for signal in (GalileoE5aI(), GalileoE5aQ(), GalileoE5aQP())
            @test get_group_delay(fnav, signal) ≈ 1.79 * BGD_E1_E5a rtol = 1e-2
        end
        @test_throws MethodError get_group_delay(fnav, GalileoE1B())

        # E6-B broadcasts no group delay.
        @test_throws MethodError get_group_delay(GalileoE6BDecoderState(1), GalileoE6B())
    end

    @testset "BeiDou" begin
        # D1/D2: T_GD1 on B1I; the clock is B3I's own.
        dnav = complete_state(BeiDouB1IDecoderState(20), BeiDouDNAVData; T_GD1 = T_GD)
        @test get_group_delay(dnav, BeiDouB1I()) == T_GD
        b3i = GNSSDecoderState(BeiDouB3IDecoderState(20); data = dnav.data)
        @test get_group_delay(b3i, BeiDouB3I()) == 0.0

        # B-CNAV1: the pilot's delay, plus the data component's ISC.
        b1c = complete_state(
            BeiDouB1CDecoderState(20),
            BeiDouB1CData;
            T_GD_B1Cp = T_GD,
            ISC_B1Cd = ISC,
        )
        @test get_group_delay(b1c, BeiDouB1C_D()) == T_GD + ISC
        @test get_group_delay(b1c, BeiDouB1C_P()) == T_GD

        # B-CNAV2 mirrors it, with its terms ungated.
        b2a = complete_state(
            BeiDouB2aDecoderState(20),
            BeiDouB2aData;
            SOW_mt11 = Int64(4),   # broadcast in the frame after MT10's
            T_GD_B2ap = T_GD,
            ISC_B2ad = ISC,
        )
        @test is_decoding_completed_for_positioning(b2a)
        @test get_group_delay(b2a, BeiDouB2aI()) == T_GD + ISC
        @test get_group_delay(b2a, BeiDouB2aQ()) == T_GD
        b2a_partial = complete_state(
            BeiDouB2aDecoderState(20),
            BeiDouB2aData;
            SOW_mt11 = Int64(4),
            T_GD_B2ap = T_GD,
            ISC_B2ad = nothing,
        )
        @test is_decoding_completed_for_positioning(b2a_partial)
        # The data component needs both terms, the pilot only its own.
        @test isnothing(get_group_delay(b2a_partial, BeiDouB2aI()))
        @test get_group_delay(b2a_partial, BeiDouB2aQ()) == T_GD

        b2b = complete_state(BeiDouB2bDecoderState(20), BeiDouB2bData; T_GD_B2bI = T_GD)
        @test get_group_delay(b2b, BeiDouB2bI()) == T_GD
    end

    @testset "the ranging signal is an instance" begin
        l1ca = complete_state(GPSL1CADecoderState(1), GPSL1CAData; T_GD)
        @test get_group_delay(l1ca, GPSL1CA()) == T_GD
        # A type is refused rather than dispatched on dynamically.
        @test_throws MethodError get_group_delay(l1ca, GPSL1CA)
    end

    @testset "allocation-free" begin
        # Read per satellite per epoch, as the clock correction is, with the
        # signal instance the tracking loop holds.
        function allocations(state, signal)
            get_group_delay(state, signal)
            @allocated get_group_delay(state, signal)
        end
        @test allocations(complete_state(GPSL1CADecoderState(1), GPSL1CAData), GPSL1CA()) ==
              0
        @test allocations(complete_state(GPSL5IDecoderState(1), GPSCNAVData), GPSL5I()) == 0
        @test allocations(
            complete_state(GalileoE5aDecoderState(1), GalileoE5aData),
            GalileoE5aQ(),
        ) == 0
    end
end
