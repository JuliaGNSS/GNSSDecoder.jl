@testset "Ionospheric parameters" begin
    klobuchar = (;
        α_0 = 1.0e-8,
        α_1 = 2.0e-8,
        α_2 = -3.0e-8,
        α_3 = 4.0e-8,
        β_0 = 1.0e5,
        β_1 = 2.0e5,
        β_2 = -3.0e5,
        β_3 = 4.0e5,
    )
    expected_klobuchar = values(klobuchar)
    nequick_g = (; a_i0 = 60.0, a_i1 = 0.3, a_i2 = 0.002)
    bdgim = (;
        α_bdgim_1 = 1.0,
        α_bdgim_2 = 2.0,
        α_bdgim_3 = 3.0,
        α_bdgim_4 = 4.0,
        α_bdgim_5 = -5.0,
        α_bdgim_6 = 6.0,
        α_bdgim_7 = 7.0,
        α_bdgim_8 = 8.0,
        α_bdgim_9 = 9.0,
    )
    fields(x) = ntuple(i -> getfield(x, i), fieldcount(typeof(x)))

    every_state = (
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

    @testset "every signal answers every model" begin
        for state in every_state
            @test isnothing(get_klobuchar_parameters(state))
            @test isnothing(get_nequick_g_parameters(state))
            @test isnothing(get_bdgim_parameters(state))
        end
    end

    @testset "Klobuchar on GPS and BeiDou D1/D2" begin
        for (state, D) in (
            (GPSL1CADecoderState(1), GPSLNAVData),
            (GPSL5IDecoderState(1), GPSCNAVData),
            (GPSL1C_DDecoderState(1), GPSCNAV2Data),
        )
            s = GNSSDecoderState(state; data = D(; klobuchar...))
            parameters = get_klobuchar_parameters(s)
            @test parameters isa KlobucharParameters
            @test fields(parameters) == expected_klobuchar
            # A partly decoded set is no set.
            partial = GNSSDecoderState(state; data = D(; klobuchar..., β_3 = nothing))
            @test isnothing(get_klobuchar_parameters(partial))
        end
        dnav = GNSSDecoderState(
            BeiDouB1IDecoderState(20);
            data = BeiDouDNAVData(; klobuchar...),
        )
        parameters = get_klobuchar_parameters(dnav)
        @test parameters isa BeiDouKlobucharParameters
        @test fields(parameters) == expected_klobuchar
    end

    @testset "NeQuick-G on Galileo I/NAV and F/NAV" begin
        for (state, D) in (
            (GalileoE1BDecoderState(1), GalileoINAVData),
            (GalileoE5bDecoderState(1), GalileoINAVData),
            (GalileoE5aDecoderState(1), GalileoFNAVData),
        )
            parameters =
                get_nequick_g_parameters(GNSSDecoderState(state; data = D(; nequick_g...)))
            @test parameters == NeQuickGParameters(60.0, 0.3, 0.002)
            # The set carries no reference epoch, so it does not wait for the
            # week number, which comes in another word type.
            @test isnothing(GNSSDecoderState(state; data = D(; nequick_g...)).data.WN)
            # A partly decoded set is no set.
            partial = GNSSDecoderState(state; data = D(; nequick_g..., a_i2 = nothing))
            @test isnothing(get_nequick_g_parameters(partial))
        end
    end

    @testset "BDGIM on BeiDou B-CNAV" begin
        for (state, D) in (
            (BeiDouB1CDecoderState(20), BeiDouCNAV1Data),
            (BeiDouB2aDecoderState(20), BeiDouCNAV2Data),
            (BeiDouB2bDecoderState(20), BeiDouCNAV3Data),
        )
            parameters = get_bdgim_parameters(GNSSDecoderState(state; data = D(; bdgim...)))
            @test parameters ==
                  BDGIMParameters(1.0, 2.0, 3.0, 4.0, -5.0, 6.0, 7.0, 8.0, 9.0)
            # No week number needed, as for NeQuick-G.
            @test isnothing(GNSSDecoderState(state; data = D(; bdgim...)).data.WN)
            partial = GNSSDecoderState(state; data = D(; bdgim..., α_bdgim_9 = nothing))
            @test isnothing(get_bdgim_parameters(partial))
        end
    end

    @testset "allocation-free" begin
        function allocations(read, state)
            read(state)
            @allocated read(state)
        end
        gps = GNSSDecoderState(GPSL1CADecoderState(1); data = GPSLNAVData(; klobuchar...))
        @test allocations(get_klobuchar_parameters, gps) == 0
        galileo = GNSSDecoderState(
            GalileoE1BDecoderState(1);
            data = GalileoINAVData(; nequick_g...),
        )
        @test allocations(get_nequick_g_parameters, galileo) == 0
        b1c =
            GNSSDecoderState(BeiDouB1CDecoderState(20); data = BeiDouCNAV1Data(; bdgim...))
        @test allocations(get_bdgim_parameters, b1c) == 0
    end
end
