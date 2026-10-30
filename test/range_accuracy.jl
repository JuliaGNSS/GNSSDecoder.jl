# The broadcast range accuracy in metres, against each ICD's table: the nominal
# value of every index, `Inf` where no prediction is broadcast, and `nothing` until
# decoded and where no single figure is broadcast.

using GNSSDecoder: GPSCNAV2Data

@testset "Broadcast range accuracy" begin
    lnav(N) =
        GNSSDecoderState(GPSL1CADecoderState(1); data = GPSLNAVData(; URA_index = Int64(N)))
    inav(index) = GNSSDecoderState(
        GalileoE1BDecoderState(1);
        data = GalileoINAVData(; SISA_E1_E5b = index),
    )
    fnav(index) = GNSSDecoderState(
        GalileoE5aDecoderState(1);
        data = GalileoFNAVData(; SISA_E1_E5a = index),
    )
    dnav(N) = GNSSDecoderState(
        BeiDouB1IDecoderState(20);
        data = BeiDouDNAVData(; URAI = Int64(N)),
    )
    # The upper bounds of the URA ranges, shared by IS-GPS-200 Table 20-I and
    # BDS-SIS-ICD-B1I-3.0 Table 5-4.
    upper = [
        2.4,
        3.4,
        4.85,
        6.85,
        9.65,
        13.65,
        24.0,
        48.0,
        96.0,
        192.0,
        384.0,
        768.0,
        1536.0,
        3072.0,
        6144.0,
    ]

    @testset "GPS LNAV URA: IS-GPS-200N §20.3.3.3.1.3, nominal values" begin
        # The nominal value of each index, the roundings IS-GPS-200 asks for included.
        nominal = [
            2.0,
            2.8,
            4.0,
            5.7,
            8.0,
            11.3,
            16.0,
            32.0,
            64.0,
            128.0,
            256.0,
            512.0,
            1024.0,
            2048.0,
            4096.0,
        ]
        @test [get_range_accuracy(lnav(N)) for N = 0:14] == nominal
        # Each lies within its index's range of Table 20-I.
        @test all(nominal .<= upper)
        @test get_range_accuracy(lnav(15)) == Inf
    end

    @testset "Galileo SISA: OS SIS ICD Table 91" begin
        for state in (inav, fnav)
            @test get_range_accuracy(state(0)) == 0.0
            @test get_range_accuracy(state(49)) ≈ 0.49
            @test get_range_accuracy(state(50)) ≈ 0.5
            @test get_range_accuracy(state(74)) ≈ 0.98
            @test get_range_accuracy(state(75)) ≈ 1.0
            @test get_range_accuracy(state(99)) ≈ 1.96
            @test get_range_accuracy(state(100)) ≈ 2.0
            @test get_range_accuracy(state(125)) ≈ 6.0
            @test get_range_accuracy(state(126)) == Inf   # spare
            @test get_range_accuracy(state(255)) == Inf   # NAPA
        end
        # The decoder's own pair: E5b decodes I/NAV, E1-E5b.
        @test get_range_accuracy(
            GNSSDecoderState(GalileoE5bDecoderState(1); data = inav(100).data),
        ) ≈ 2.0
    end

    @testset "BeiDou D1/D2 URA: BDS-SIS-ICD-B1I-3.0 Table 5-4" begin
        values = [get_range_accuracy(dnav(N)) for N = 0:14]
        @test values == [N < 6 ? 2.0^(N / 2 + 1) : 2.0^(N - 2) for N = 0:14]
        @test all(values .<= upper)
        @test get_range_accuracy(dnav(15)) == Inf
        @test get_range_accuracy(
            GNSSDecoderState(BeiDouB3IDecoderState(20); data = dnav(3).data),
        ) == get_range_accuracy(dnav(3))
    end

    @testset "nothing until decoded, and where no single figure is broadcast" begin
        for decoder in (
            GPSL1CADecoderState(1),
            GPSL5IDecoderState(1),
            GPSL2CMDecoderState(1),
            GPSL1C_DDecoderState(1),
            GalileoE1BDecoderState(1),
            GalileoE5bDecoderState(1),
            GalileoE5aDecoderState(1),
            GalileoE6BDecoderState(1),
            BeiDouB1IDecoderState(20),
            BeiDouB1CDecoderState(20),
            BeiDouB2aDecoderState(20),
            BeiDouB2bDecoderState(20),
        )
            @test isnothing(get_range_accuracy(decoder))
        end
        # CNAV's and CNAV-2's URA combines an elevation- and a time-dependent
        # component, so decoded indices still give no single figure.
        @test isnothing(
            get_range_accuracy(
                GNSSDecoderState(
                    GPSL5IDecoderState(1);
                    data = GPSCNAVData(;
                        ura_ed_index = Int64(0),
                        ura_ned0_index = Int64(0),
                    ),
                ),
            ),
        )
        @test isnothing(
            get_range_accuracy(
                GNSSDecoderState(
                    GPSL1C_DDecoderState(1);
                    data = GPSCNAV2Data(;
                        ura_ed_index = Int64(0),
                        ura_ned0_index = Int64(0),
                    ),
                ),
            ),
        )
    end

    @testset "allocation-free" begin
        function allocations(state)
            get_range_accuracy(state)
            @allocated get_range_accuracy(state)
        end
        @test allocations(lnav(5)) == 0
        @test allocations(inav(100)) == 0
        @test allocations(dnav(15)) == 0
    end
end
