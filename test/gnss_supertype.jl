@testset "GNSS data supertypes" begin
    # The per-constellation abstract supertypes sit between `AbstractGNSSData`
    # and the concrete per-signal data types.
    @test GNSSDecoder.AbstractGPSData <: GNSSDecoder.AbstractGNSSData
    @test GNSSDecoder.AbstractGalileoData <: GNSSDecoder.AbstractGNSSData
    # ...and one level further down for the Galileo signals that carry an
    # ephemeris, so the shared completeness checks dispatch on a type that has
    # the orbital fields.
    @test GNSSDecoder.AbstractGalileoEphemerisData <: GNSSDecoder.AbstractGalileoData

    gps_data = [GNSSDecoder.GPSL1CAData, GNSSDecoder.GPSL1C_DData, GNSSDecoder.GPSCNAVData]
    # The Galileo data types that carry an ephemeris. `GalileoINAVData` is shared
    # by E1-B and E5b; `GalileoE6BData` deliberately sits outside this list (see
    # below).
    galileo_data = [GNSSDecoder.GalileoINAVData, GNSSDecoder.GalileoE5aData]

    # Constellation membership is encoded at each struct's definition site.
    for D in gps_data
        @test D <: GNSSDecoder.AbstractGPSData
        @test !(D <: GNSSDecoder.AbstractGalileoData)
    end
    for D in galileo_data
        @test D <: GNSSDecoder.AbstractGalileoData
        @test !(D <: GNSSDecoder.AbstractGPSData)
    end

    # `is_ephemeris_decoded` / `is_clock_correction_decoded` are stated once for
    # the ephemeris-bearing Galileo signals: both dispatch to the single method on
    # `AbstractGalileoEphemerisData` rather than a per-signal copy.
    for D in galileo_data
        @test D <: GNSSDecoder.AbstractGalileoEphemerisData
        @test which(GNSSDecoder.is_ephemeris_decoded, (D,)).sig == Tuple{
            typeof(GNSSDecoder.is_ephemeris_decoded),
            GNSSDecoder.AbstractGalileoEphemerisData,
        }
        @test which(GNSSDecoder.is_clock_correction_decoded, (D,)).sig == Tuple{
            typeof(GNSSDecoder.is_clock_correction_decoded),
            GNSSDecoder.AbstractGalileoEphemerisData,
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

    # Galileo E6-B is a Galileo signal but carries no ephemeris of its own — its
    # C/NAV message corrects *another* signal's — so it joins the constellation
    # supertype without joining the ephemeris-bearing one, and states its own
    # positioning-readiness answer. It must not pick up the completeness checks:
    # they field-access `t_0e` / `M_0` / `t_0c`, which it does not have.
    @test GalileoE6BData <: GNSSDecoder.AbstractGalileoData
    @test !(GalileoE6BData <: GNSSDecoder.AbstractGPSData)
    @test !(GalileoE6BData <: GNSSDecoder.AbstractGalileoEphemerisData)
    @test !hasmethod(GNSSDecoder.is_ephemeris_decoded, (GalileoE6BData,))
    @test !hasmethod(GNSSDecoder.is_clock_correction_decoded, (GalileoE6BData,))
    @test !is_decoding_completed_for_positioning(GalileoE6BData())
    @test which(is_decoding_completed_for_positioning, (GalileoE6BData,)).sig ==
          Tuple{typeof(is_decoding_completed_for_positioning),GalileoE6BData}

    # Containers shared by two signals were renamed away from the first signal
    # that happened to implement them — `GPSL5I*` → `GPSCNAV*` when L2C joined,
    # `GalileoE1B*` → `GalileoINAV*` when E5b-I did. Both sets of old names are
    # gone outright rather than deprecated: this is the breaking release, so the
    # migration happens once.
    for old in (
        :GalileoE1BData,
        :GalileoE1BCache,
        :GPSL5IData,
        :GPSL5IReducedAlmanac,
        :GPSL5IMidiAlmanac,
        :GPSL5IClockDifferentialCorrection,
        :GPSL5IEphemerisDifferentialCorrection,
        :GPSL5IIntegritySupportMessage,
    )
        @test !isdefined(GNSSDecoder, old)
    end

    # Every Galileo data container is exported, as the GPS CNAV and BeiDou ones
    # are — a caller naming `GNSSDecoderState{GalileoINAVData}` should not have to
    # qualify it.
    for D in (:GalileoINAVData, :GalileoE5aData, :GalileoE6BData)
        @test D in names(GNSSDecoder)
    end
end
