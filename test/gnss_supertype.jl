@testset "GNSS data supertypes" begin
    # The per-constellation abstract supertypes sit between `AbstractGNSSData`
    # and the concrete per-signal data types.
    @test GNSSDecoder.AbstractGPSData <: GNSSDecoder.AbstractGNSSData
    @test GNSSDecoder.AbstractGalileoData <: GNSSDecoder.AbstractGNSSData

    gps_data = [GNSSDecoder.GPSL1CAData, GNSSDecoder.GPSL1C_DData, GNSSDecoder.GPSCNAVData]
    galileo_data = [GNSSDecoder.GalileoE1BData, GNSSDecoder.GalileoE5aData]

    # Constellation membership is encoded at each struct's definition site.
    for D in gps_data
        @test D <: GNSSDecoder.AbstractGPSData
        @test !(D <: GNSSDecoder.AbstractGalileoData)
    end
    for D in galileo_data
        @test D <: GNSSDecoder.AbstractGalileoData
        @test !(D <: GNSSDecoder.AbstractGPSData)
    end

    # `is_ephemeris_decoded` / `is_clock_correction_decoded` are stated once per
    # constellation: every Galileo data type dispatches to the single method on
    # `AbstractGalileoData` rather than a per-signal copy.
    for D in galileo_data
        @test which(GNSSDecoder.is_ephemeris_decoded, (D,)).sig == Tuple{
            typeof(GNSSDecoder.is_ephemeris_decoded),
            GNSSDecoder.AbstractGalileoData,
        }
        @test which(GNSSDecoder.is_clock_correction_decoded, (D,)).sig == Tuple{
            typeof(GNSSDecoder.is_clock_correction_decoded),
            GNSSDecoder.AbstractGalileoData,
        }
    end

    # The collapsed methods still behave: all fields present ⇒ decoded, a
    # missing field ⇒ not decoded. Exercised for both Galileo signals.
    for D in galileo_data
        empty = D()
        @test !GNSSDecoder.is_ephemeris_decoded(empty)
        @test !GNSSDecoder.is_clock_correction_decoded(empty)

        full_eph = D(;
            t_0e = 0,
            M_0 = 0.0,
            e = 0.0,
            sqrt_A = 0.0,
            Ω_0 = 0.0,
            i_0 = 0.0,
            ω = 0.0,
            i_dot = 0.0,
            Ω_dot = 0.0,
            Δn = 0.0,
            C_uc = 0.0,
            C_us = 0.0,
            C_rc = 0.0,
            C_rs = 0.0,
            C_ic = 0.0,
            C_is = 0.0,
        )
        @test GNSSDecoder.is_ephemeris_decoded(full_eph)

        full_clock = D(; t_0c = 0, a_f0 = 0.0, a_f1 = 0.0, a_f2 = 0.0)
        @test GNSSDecoder.is_clock_correction_decoded(full_clock)
    end
end
