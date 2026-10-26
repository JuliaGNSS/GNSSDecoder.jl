# A cache type with no `clear_satellite_cache` method, as a new signal might add.
struct CacheWithoutMethod <: GNSSDecoder.AbstractGNSSCache end

@testset "Reset for another satellite" begin
    @testset "every decoder" begin
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
            D = typeof(state.data)
            # A decoder in the middle of decoding one satellite.
            used = GNSSDecoderState(
                state;
                raw_data = complete_data(D),
                data = complete_data(D),
                num_symbols_since_sync = 5,
                is_shifted_by_180_degrees = true,
            )
            push!(used.cache.soft_buffer, 1.0f0)
            moved = reset_decoder_state!(used, 9)
            @test moved.prn == 9
            @test moved.raw_data == D()
            @test moved.data == D()
            @test isnothing(moved.num_symbols_since_sync)
            @test !moved.is_shifted_by_180_degrees
            @test isempty(moved.cache.soft_buffer)
            # The buffers are kept, not allocated anew.
            @test moved.cache.soft_buffer === used.cache.soft_buffer
            @test moved.cache.storage === used.cache.storage
            @test typeof(moved) === typeof(used)
        end
    end

    @testset "the GPS L1 C/A vote tally is cleared" begin
        state = GPSL1CADecoderState(1)
        push!(state.cache.old_data, GNSSDecoder.VotedGPSL1CAData(1, GPSL1CAData()))
        # Kept for the same satellite, for a fast recovery ...
        @test !isempty(reset_decoder_state!(copy(state)).cache.old_data)
        # ... and cleared for another one.
        @test isempty(reset_decoder_state!(state, 2).cache.old_data)
    end

    @testset "a cache without a method is refused" begin
        # No fallback: a new cache must state what of it belongs to the satellite.
        @test_throws MethodError GNSSDecoder.clear_satellite_cache(CacheWithoutMethod())
    end

    @testset "message assembly in flight is cleared" begin
        inav = GalileoE1BDecoderState(1)
        inav = GNSSDecoderState(
            inav;
            cache = GNSSDecoder.GalileoINAVCache(
                inav.cache;
                even_page_part_bits = UInt128(5),
                almanac_chain_pos1 = GNSSDecoder.GalileoAlmanac(; SVID = 3),
            ),
        )
        moved = reset_decoder_state!(inav, 2)
        @test isnothing(moved.cache.even_page_part_bits)
        @test moved.cache.almanac_chain_pos1 == GNSSDecoder.GalileoAlmanac()

        fnav = GalileoE5aDecoderState(1)
        fnav = GNSSDecoderState(
            fnav;
            cache = GNSSDecoder.GalileoE5aCache(
                fnav.cache;
                almanac_chain_partial = GNSSDecoder.GalileoAlmanac(; SVID = 3),
                almanac_chain_omega0_msb = 7,
            ),
        )
        moved = reset_decoder_state!(fnav, 2)
        @test moved.cache.almanac_chain_partial == GNSSDecoder.GalileoAlmanac()
        @test isnothing(moved.cache.almanac_chain_omega0_msb)
    end

    @testset "allocation-free" begin
        function allocations(state, prn)
            reset_decoder_state!(copy(state), prn)
            state = copy(state)
            @allocated reset_decoder_state!(state, prn)
        end
        for state in (
            GPSL1CADecoderState(1),
            GalileoE1BDecoderState(1),
            GalileoE5aDecoderState(1),
            BeiDouB1IDecoderState(20),
            BeiDouB1CDecoderState(20),
        )
            @test allocations(state, 9) == 0 skip = !CHECK_ALLOCATIONS
        end
    end
end
