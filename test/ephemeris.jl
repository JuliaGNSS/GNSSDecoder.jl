@testset "Broadcast ephemeris" begin
    # Distinct values for the elements every message names alike, so a field
    # read from the wrong place shows.
    shared = (;
        t_0e = Int64(7200),
        M_0 = 0.1,
        e = 0.01,
        ω = 0.2,
        Ω_0 = 0.3,
        i_0 = 0.9,
        i_dot = 1.0e-10,
        C_rs = 10.0,
        C_rc = 200.0,
        C_us = 1.0e-6,
        C_uc = 2.0e-6,
        C_is = 3.0e-8,
        C_ic = 4.0e-8,
        t_0c = Int64(7000),
        a_f0 = 1.0e-4,
        a_f1 = 1.0e-11,
        a_f2 = 1.0e-18,
    )

    function test_shared(ephemeris, state)
        for name in keys(shared)
            @test getfield(ephemeris, name) == getfield(shared, name)
        end
        # The constants each ICD fixes for its ephemeris algorithm: WGS-84's for GPS
        # (IS-GPS-200 Table 20-IV), Galileo's own (OS SIS ICD Table 68), BDCS's for
        # BeiDou (BDS-SIS-ICD-B1I-3.0 §5.2.4.12), whose rotation rate differs.
        μ, F, Ω_dot_e =
            state.data isa GNSSDecoder.AbstractGPSData ?
            (3.986005e14, -4.442807633e-10, 7.2921151467e-5) :
            state.data isa GNSSDecoder.AbstractGalileoData ?
            (3.986004418e14, -4.442807309e-10, 7.2921151467e-5) :
            (3.986004418e14, -4.442807309e-10, 7.2921150e-5)
        @test ephemeris.μ == μ
        @test ephemeris.F == F
        @test ephemeris.Ω_dot_e == Ω_dot_e
    end

    @testset "nothing until the positioning set is complete" begin
        for state in (
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
        )
            @test isnothing(get_broadcast_ephemeris(state))
        end
    end

    @testset "Keplerian messages" begin
        keplerian = (; shared..., sqrt_A = 5153.6, Δn = 4.0e-9, Ω_dot = -8.0e-9)
        for (state, D) in (
            (GPSL1CADecoderState(1), GPSL1CAData),
            (GalileoE1BDecoderState(1), GalileoINAVData),
            (GalileoE5aDecoderState(1), GalileoE5aData),
            (BeiDouB1IDecoderState(20), BeiDouDNAVData),
        )
            s = complete_state(state, D; keplerian...)
            ephemeris = get_broadcast_ephemeris(s)
            @test ephemeris isa BroadcastEphemeris
            test_shared(ephemeris, s)
            @test ephemeris.A == 5153.6^2
            @test ephemeris.Δn_0 == 4.0e-9
            @test ephemeris.Ω_dot == -8.0e-9
            @test ephemeris.A_dot == ephemeris.Δn_0_dot == 0.0
            @test !ephemeris.is_geostationary
        end
        # BeiDou D1/D2 tells its GEOs apart by PRN.
        geo = complete_state(BeiDouB1IDecoderState(1), BeiDouDNAVData; keplerian...)
        @test get_broadcast_ephemeris(geo).is_geostationary
    end

    @testset "GPS CNAV and CNAV-2" begin
        quasi = (;
            shared...,
            ΔA = 120.0,
            A_dot = 0.01,
            Δn_0 = 4.0e-9,
            Δn_0_dot = 1.0e-14,
            ΔΩ_dot = 1.0e-10,
        )
        for (state, D) in (
            (GPSL5IDecoderState(1), GPSCNAVData),
            (GPSL2CMDecoderState(1), GPSCNAVData),
            (GPSL1C_DDecoderState(1), GPSL1C_DData),
        )
            s = complete_state(state, D; quasi...)
            ephemeris = get_broadcast_ephemeris(s)
            test_shared(ephemeris, s)
            # Off A_REF and Ω̇_REF (IS-GPS-200N §30.3.3.1.3).
            @test ephemeris.A == 26_559_710.0 + 120.0
            @test ephemeris.Ω_dot == -2.6e-9 * GNSSDecoder.GNSS_PI + 1.0e-10
            @test (ephemeris.A_dot, ephemeris.Δn_0, ephemeris.Δn_0_dot) ==
                  (0.01, 4.0e-9, 1.0e-14)
            @test !ephemeris.is_geostationary
        end
    end

    @testset "BeiDou B-CNAV" begin
        quasi = (;
            shared...,
            ΔA = -120.0,
            A_dot = 0.01,
            Δn_0 = 4.0e-9,
            Δn_0_dot = 1.0e-14,
            Ω_dot = -7.0e-9,
        )
        for (state, D, overrides) in (
            (BeiDouB1CDecoderState(20), BeiDouB1CData, (;)),
            (BeiDouB2aDecoderState(20), BeiDouB2aData, (; SOW_mt11 = Int64(4))),
            (BeiDouB2bDecoderState(20), BeiDouB2bData, (;)),
        )
            # The reference semi-major axis follows the orbit class: MEO (3),
            # IGSO (2) and GEO (1).
            for (sat_type, A_REF) in
                ((3, 27_906_100.0), (2, 42_162_200.0), (1, 42_162_200.0))
                s = complete_state(
                    state,
                    D;
                    quasi...,
                    overrides...,
                    sat_type = Int64(sat_type),
                )
                ephemeris = get_broadcast_ephemeris(s)
                test_shared(ephemeris, s)
                @test ephemeris.A == A_REF - 120.0
                # Broadcast outright, not off a reference.
                @test ephemeris.Ω_dot == -7.0e-9
                @test ephemeris.is_geostationary == (sat_type == 1)
            end
        end
        # The reserved orbit type has no reference axis; the gate keeps it out.
        @test_throws ArgumentError GNSSDecoder.beidou_cnav_reference_semi_major_axis(0)
        @test !is_decoding_completed_for_positioning(
            complete_state(BeiDouB1CDecoderState(20), BeiDouB1CData; sat_type = Int64(0)),
        )
    end

    @testset "a bits type, compared by value, built without allocating" begin
        @test isbitstype(BroadcastEphemeris)
        a = complete_state(GPSL1CADecoderState(1), GPSL1CAData)
        b = complete_state(GPSL1CADecoderState(1), GPSL1CAData)
        @test get_broadcast_ephemeris(a) == get_broadcast_ephemeris(b)
        c = complete_state(GPSL1CADecoderState(1), GPSL1CAData; M_0 = 0.5)
        @test get_broadcast_ephemeris(a) != get_broadcast_ephemeris(c)

        function allocations(state)
            get_broadcast_ephemeris(state)
            @allocated get_broadcast_ephemeris(state)
        end
        @test allocations(a) == 0
        @test allocations(complete_state(GPSL5IDecoderState(1), GPSCNAVData)) == 0
        @test allocations(
            complete_state(BeiDouB2aDecoderState(20), BeiDouB2aData; SOW_mt11 = Int64(4)),
        ) == 0
    end
end
