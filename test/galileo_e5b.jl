# Galileo E5b I/NAV decoder tests.
#
# E5b-I carries the *same* I/NAV message as E1-B — the OS SIS ICD (Issue 2.2
# §4.3.1) says so outright: "the same page layout … only page sequencing is
# different" — so the decoder is the shared core in `src/galileo/inav.jl` and
# `e5b.jl` is a thin signal layer. These tests therefore check two distinct
# things:
#
#  1. That the signal layer is right: the E5b-I signal is reported (not E1-B),
#     its metadata comes out of GNSSSignals, and `is_sat_healthy` reads the *E5b*
#     facet of word type 5 rather than the E1-B/C one. That health facet is the
#     only decode-level difference between the two decoders.
#  2. That the shared core really is layout-compatible with an E5b page pair.
#     The page parts differ from E1-B in fields this decoder never reads — a
#     64-bit "Reserved 1" where E1-B has OSNMA(40) + SAR(22) + Spare(2), and an
#     8-bit "Reserved 2" where E1-B has the SSP (Table 38) — so the tests below
#     drive synthetic page pairs with random content in both, and pin down the
#     CRC scope by tampering: a flipped Reserved 1 bit must break the page (it is
#     CRC-protected), a flipped Reserved 2 bit must not (it is not).
#
# For (1) the golden E1-B symbol fixture from `galileo_e1b.jl` is replayed
# through an E5b decoder: identical symbols must yield identical navigation data,
# which is the claim "one decoder, two signals" stated as a test.

"""
I/NAV page synchronisation pattern `0101100000` (OS SIS ICD §4.3.2.1).
"""
const GALILEO_INAV_SYNC = (0, 1, 0, 1, 1, 0, 0, 0, 0, 0)

"""
    inav_word_type_5(; kwargs...) -> Vector{Bool}

The 128 data bits of an I/NAV word type 5 (OS SIS ICD Table 40 / §5.1.6-5.1.10):
ionospheric coefficients, storm flags, both broadcast group delays, the E5b and
E1-B/C signal-health and data-validity fields, and WN/TOW. Chosen for these
tests because it is the one word carrying the health facets the two signal layers
disagree about.
"""
function inav_word_type_5(;
    a_i0_raw = 300,
    a_i1_raw = 12,
    a_i2_raw = -7,
    storm_flags = (false, false, false, false, false),
    bgd_e1_e5a_raw = 24,
    bgd_e1_e5b_raw = -31,
    signal_health_e5b = 0,
    signal_health_e1b = 0,
    data_validity_e5b = 0,
    data_validity_e1b = 0,
    WN = 1234,
    TOW = 345_600,
)
    bits = Bool[]
    galileo_push_field!(bits, 5, 6)                        # word type
    galileo_push_field!(bits, a_i0_raw, 11)
    galileo_push_field!(bits, a_i1_raw & 0x7ff, 11)
    galileo_push_field!(bits, a_i2_raw & 0x3fff, 14)
    for flag in storm_flags
        push!(bits, flag)
    end
    galileo_push_field!(bits, bgd_e1_e5a_raw & 0x3ff, 10)
    galileo_push_field!(bits, bgd_e1_e5b_raw & 0x3ff, 10)
    galileo_push_field!(bits, signal_health_e5b, 2)
    galileo_push_field!(bits, signal_health_e1b, 2)
    galileo_push_field!(bits, data_validity_e5b, 1)
    galileo_push_field!(bits, data_validity_e1b, 1)
    galileo_push_field!(bits, WN, 12)
    galileo_push_field!(bits, TOW, 20)
    galileo_push_field!(bits, 0, 23)                       # spare
    @assert length(bits) == 128
    bits
end

"""
    inav_page_pair(word; reserved_1, reserved_2) -> (even_bits, odd_bits)

Split one 128-bit I/NAV word into the 114 information bits of an even and an odd
nominal page part in the *E5b-I* layout (OS SIS ICD Table 38):

    even: Even/Odd=0 | Page Type=0 | Data(1/2) 112
    odd:  Even/Odd=1 | Page Type=0 | Data(2/2) 16 | Reserved 1 (64) | CRC (24) | Reserved 2 (8)

The CRC-24Q covers the even part's 114 bits plus the odd part's leading 82 —
everything up to and including Reserved 1, and nothing of Reserved 2 (§4.3.2.3).
"""
function inav_page_pair(
    word::Vector{Bool};
    reserved_1::Vector{Bool},
    reserved_2::Vector{Bool},
)
    @assert length(word) == 128
    @assert length(reserved_1) == 64
    @assert length(reserved_2) == 8
    even = vcat([false, false], word[1:112])
    @assert length(even) == 114
    odd_prefix = vcat([true, false], word[113:128], reserved_1)
    @assert length(odd_prefix) == 82
    crc = Bool[]
    galileo_push_field!(crc, crc24q(vcat(even, odd_prefix)), 24)
    odd = vcat(odd_prefix, crc, reserved_2)
    @assert length(odd) == 114
    return even, odd
end

"""
Concatenate I/NAV page parts (250 on-air symbols each: 10-symbol sync + 240
interleaved) into a soft-symbol stream, closed by a trailing sync.
"""
inav_symbol_stream(parts::Vector{Vector{Bool}}) =
    galileo_symbol_stream(parts, GALILEO_INAV_SYNC, 30)

@testset "Galileo E5b constructor" begin
    decoder = GalileoE5bDecoderState(21)
    @test decoder.prn == 21
    @test decoder.data isa GNSSDecoder.GalileoINAVData
    @test decoder.constants isa GNSSDecoder.GalileoE5bConstants
    # Framing is I/NAV's, identical to E1-B (OS SIS ICD Table 37).
    @test decoder.constants.syncro_sequence_length == 250
    @test decoder.constants.preamble_length == 10
    @test decoder.constants.preamble == 0b0101100000
    # Same field values as E1-B's constants — only the type tag differs.
    let e1b_constants = GalileoE1BDecoderState(21).constants
        @test all(
            getfield(decoder.constants, f) == getfield(e1b_constants, f) for
            f in fieldnames(GNSSDecoder.GalileoINAVConstants)
        )
    end
    # ...but the signal is E5b-I, the data-bearing component (E5b-Q is the pilot).
    @test get_signal_type(decoder) == GalileoE5bI
    @test get_band(decoder) == GNSSSignals.E5b()
    @test get_signal_name(decoder) == "Galileo E5b-I"
    @test get_data_frequency(decoder) == 250GNSSSignals.Hz
    @test get_time_system(decoder) == get_time_system(GalileoE1BDecoderState(21))
    @test hasmethod(GNSSDecoderState, Tuple{GalileoE5bI,Int})
    @test !hasmethod(GNSSDecoderState, Tuple{GalileoE5bQ,Int})
    @test !is_decoding_completed_for_positioning(decoder)
    @test !is_sat_healthy(decoder)
end

@testset "Galileo E5b decodes the E1B golden fixture identically" begin
    # The same I/NAV symbols must decode the same way on either component: the
    # page layout, FEC, interleaver and CRC scope are shared, so a divergence
    # here would mean the signal layer had leaked into the message decoding.
    e1b = reduce(
        (dec, data) ->
            decode(dec, to_soft_symbols(data, sizeof(data) * 8), sizeof(data) * 8),
        GALILEO_E1B_DATA;
        init = GalileoE1BDecoderState(21),
    )
    e5b = reduce(
        (dec, data) ->
            decode(dec, to_soft_symbols(data, sizeof(data) * 8), sizeof(data) * 8),
        GALILEO_E1B_DATA;
        init = GalileoE5bDecoderState(21),
    )
    @test e5b.data == e1b.data
    @test e5b.raw_data == e1b.raw_data
    @test e5b.num_bits_after_valid_syncro_sequence ==
          e1b.num_bits_after_valid_syncro_sequence
    @test e5b.is_shifted_by_180_degrees == e1b.is_shifted_by_180_degrees
    # Spot-check absolute values too, so this is not merely "equal to whatever
    # E1-B produced": these are the fixture's known TOW/WN and health.
    @test e5b.data.TOW == 259_235
    @test e5b.data.WN == 1082
    @test e5b.data.signal_health_e5b == GNSSDecoder.signal_ok
    @test e5b.data.data_validity_status_e5b == GNSSDecoder.navigation_data_valid
    @test is_sat_healthy(e5b)
    # Word type 16 (Reduced CED) is broadcast on E1-B only (Table 40); the
    # fixture is an E1-B capture, so it is present here. The parser is shared and
    # deliberately not gated on the component — the word is CRC-protected either
    # way, and the ICD warns its dissemination sequence is only indicative.
    @test e5b.data.reduced_ced == e1b.data.reduced_ced
end

@testset "Galileo E5b page-part layout and CRC scope" begin
    rng = MersenneTwister(20260821)
    reserved_1 = rand(rng, Bool, 64)
    reserved_2 = rand(rng, Bool, 8)
    word = inav_word_type_5(;
        signal_health_e5b = 0,
        signal_health_e1b = 1,
        data_validity_e5b = 0,
        data_validity_e1b = 1,
        WN = 1234,
        TOW = 345_600,
    )
    even, odd = inav_page_pair(word; reserved_1, reserved_2)

    decoder = GalileoE5bDecoderState(21)
    decoder = decode(
        decoder,
        inav_symbol_stream([even, odd]),
        length(inav_symbol_stream([even, odd])),
    )
    # Arbitrary Reserved 1 / Reserved 2 content must not disturb the decode.
    @test decoder.raw_data.WN == 1234
    @test decoder.raw_data.TOW == 345_600
    @test decoder.raw_data.signal_health_e5b == GNSSDecoder.signal_ok
    @test decoder.raw_data.signal_health_e1b == GNSSDecoder.signal_out_of_service
    @test decoder.raw_data.data_validity_status_e5b == GNSSDecoder.navigation_data_valid
    @test decoder.raw_data.data_validity_status_e1b == GNSSDecoder.working_without_guarantee
    @test decoder.raw_data.a_i0 == 300 / 4
    @test decoder.raw_data.broadcast_group_delay_e1_e5a == 24 / 2^32
    @test decoder.raw_data.broadcast_group_delay_e1_e5b == -31 / 2^32

    @testset "Reserved 2 is not CRC-protected" begin
        # Flipping a Reserved 2 bit leaves the checksum valid, so the word must
        # still be decoded (OS SIS ICD §4.3.2.3 note).
        tampered = copy(odd)
        tampered[107] = !tampered[107]
        state = GalileoE5bDecoderState(21)
        stream = inav_symbol_stream([even, tampered])
        state = decode(state, stream, length(stream))
        @test state.raw_data.TOW == 345_600
    end

    @testset "Reserved 1 is CRC-protected" begin
        # Flipping a Reserved 1 bit invalidates the checksum, so the whole page
        # pair must be rejected and nothing decoded from it.
        tampered = copy(odd)
        tampered[20] = !tampered[20]
        state = GalileoE5bDecoderState(21)
        stream = inav_symbol_stream([even, tampered])
        state = decode(state, stream, length(stream))
        @test isnothing(state.raw_data.TOW)
        @test isnothing(state.raw_data.WN)
    end

    @testset "the same page pair decodes on E1-B" begin
        # The 64 + 8 reserved bits are where the components differ, so decoding
        # an E5b-shaped pair with the E1-B decoder must give the same word: the
        # CRC-protected prefix is 82 bits either way.
        state = GalileoE1BDecoderState(21)
        stream = inav_symbol_stream([even, odd])
        state = decode(state, stream, length(stream))
        @test state.raw_data.TOW == 345_600
        @test state.raw_data.WN == 1234
    end
end

@testset "Galileo E5b health uses the E5b facet of word type 5" begin
    # Word type 5 carries both facets, and the two signal layers differ *only* in
    # which one they report (OS SIS ICD Tables 81 and 84).
    rng = MersenneTwister(1)
    reserved_1 = rand(rng, Bool, 64)
    reserved_2 = rand(rng, Bool, 8)
    # E5b healthy, E1-B/C out of service and working without guarantee.
    word = inav_word_type_5(;
        signal_health_e5b = 0,
        signal_health_e1b = 1,
        data_validity_e5b = 0,
        data_validity_e1b = 1,
    )
    even, odd = inav_page_pair(word; reserved_1, reserved_2)
    stream = inav_symbol_stream([even, odd])
    e5b = decode(GalileoE5bDecoderState(21), stream, length(stream))
    e1b = decode(GalileoE1BDecoderState(21), stream, length(stream))
    # `is_sat_healthy` reads the *validated* data, which needs a full positioning
    # set, so drive it directly off a state whose `data` is this word.
    e5b = GNSSDecoder.GNSSDecoderState(e5b; data = e5b.raw_data)
    e1b = GNSSDecoder.GNSSDecoderState(e1b; data = e1b.raw_data)
    @test is_sat_healthy(e5b)
    @test !is_sat_healthy(e1b)

    # And the other way round.
    word = inav_word_type_5(;
        signal_health_e5b = 3,
        signal_health_e1b = 0,
        data_validity_e5b = 1,
        data_validity_e1b = 0,
    )
    even, odd = inav_page_pair(word; reserved_1, reserved_2)
    stream = inav_symbol_stream([even, odd])
    e5b = decode(GalileoE5bDecoderState(21), stream, length(stream))
    e1b = decode(GalileoE1BDecoderState(21), stream, length(stream))
    e5b = GNSSDecoder.GNSSDecoderState(e5b; data = e5b.raw_data)
    e1b = GNSSDecoder.GNSSDecoderState(e1b; data = e1b.raw_data)
    @test !is_sat_healthy(e5b)
    @test is_sat_healthy(e1b)
    @test e5b.raw_data.signal_health_e5b == GNSSDecoder.signal_component_currently_in_test
end

@testset "Galileo E5b reset_decoder_state" begin
    state = reduce(
        (dec, data) ->
            decode(dec, to_soft_symbols(data, sizeof(data) * 8), sizeof(data) * 8),
        GALILEO_E1B_DATA;
        init = GalileoE5bDecoderState(21),
    )
    @test !isnothing(state.data.TOW)
    state = reset_decoder_state(state)
    @test GNSSDecoder.num_bits_buffered(state) == 0
    @test isnothing(state.raw_data.TOW)
    @test isnothing(state.data.TOW)
    @test isnothing(state.num_bits_after_valid_syncro_sequence)
    # Ephemeris survives the reset, as on E1-B.
    @test !isnothing(state.raw_data.sqrt_A)
end
