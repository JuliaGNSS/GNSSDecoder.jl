# Galileo E6-B C/NAV (High Accuracy Service) decoder tests.
#
# Ground truth is the Galileo HAS SIS ICD's own worked examples, transcribed in
# `has_test_vectors.jl`: Annex C's sample C/NAV page and Annex D's two complete
# HPVRS decoding examples, including the field values the ICD states they decode
# to. The Reed-Solomon layer underneath is tested separately in
# `reed_solomon.jl` (against the ICD's Annex B generator matrix and Table 42
# generator polynomial).
#
# The ICD gives decoded pages, not E6-B symbols, so — as in the E5a tests — each
# page is re-encoded through the transmit chain (rate-1/2 K=7 NSC convolutional
# code with G2 inverted, 123×8 block interleave, 16-symbol sync prefix, shared
# helpers in `galileo_test_utils.jl`) to synthesise the on-air soft-symbol
# stream, and the full `decode` path is driven from there. Because the *decoded
# field values* are compared against the ICD's independent statement of them, a
# shared encode/decode error cannot mask a field-layout mistake; and for Annex C
# the 486 page bits come from the ICD verbatim, CRC included, so the CRC scope is
# checked against a checksum this package did not compute.

using GNSSDecoder: crc24q

"""
C/NAV synchronisation pattern `1011011101110000` (HAS SIS ICD §2.3.1).
"""
const GALILEO_E6B_SYNC = (1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0)

"""
    e6b_page_bits(header, encoded_page; reserved) -> Vector{Bool}

Build the 486 information bits of one C/NAV page from a 24-bit HAS Page Header
and 53 HAS Encoded Page octets, appending the CRC-24Q over the leading 462 bits
(HAS SIS ICD Table 5, §2.3.3). The 14 reserved bits default to all ones, as they
are in the ICD's Annex C sample page.
"""
function e6b_page_bits(
    header::Integer,
    encoded_page::AbstractVector{UInt8};
    reserved = 0x3fff,
)
    bits = Bool[]
    galileo_push_field!(bits, reserved, 14)
    galileo_push_field!(bits, header, 24)
    for octet in encoded_page
        galileo_push_field!(bits, octet, 8)
    end
    @assert length(bits) == 462
    galileo_push_field!(bits, crc24q(bits), 24)
    @assert length(bits) == 486
    bits
end

"""
Compose a 24-bit HAS Page Header (HAS SIS ICD Table 7: HASS, reserved, MT, MID, MS, PID).
"""
function e6b_header(;
    has_status = 1,
    message_type = 1,
    message_id = 3,
    message_size = 4,
    page_id = 40,
)
    bits = Bool[]
    galileo_push_field!(bits, has_status, 2)
    galileo_push_field!(bits, 0, 2)                    # reserved
    galileo_push_field!(bits, message_type, 2)
    galileo_push_field!(bits, message_id, 5)
    galileo_push_field!(bits, message_size - 1, 5)     # "0" = 1 page (Table 8)
    galileo_push_field!(bits, page_id, 8)
    value = UInt32(0)
    for b in bits
        value = (value << 1) | UInt32(b)
    end
    value
end

"""
Concatenate C/NAV pages (1000 on-air symbols each, from 486 information bits)
into a soft-symbol stream, closed by a trailing sync pattern so the final page's
window is bounded by the next page's sync.
"""
e6b_symbol_stream(pages::Vector{Vector{Bool}}) =
    galileo_symbol_stream(pages, GALILEO_E6B_SYNC, 123)

"""
The 486 bits of the ICD's Annex C sample page, taken verbatim from its hex.
"""
function e6b_annex_c_page_bits()
    value = parse(BigInt, HAS_ANNEX_C_PAGE_HEX; base = 16)
    Bool[(value >> (512 - i)) & 1 == 1 for i = 1:486]
end

"""
Pages of one Annex D example, wrapped in HAS Page Headers. `has_status` defaults
to operational so the decoder reports the satellite as usable.
"""
function e6b_example_pages(example; has_status = 1, message_id = example.message_id)
    map(enumerate(example.page_ids)) do (index, page_id)
        octets = hex2bytes(example.encoded_pages_hex[index])
        header = e6b_header(;
            has_status,
            message_type = example.message_type,
            message_id,
            message_size = example.message_size,
            page_id,
        )
        e6b_page_bits(header, octets)
    end
end

"""
Drive `decode` over a whole symbol stream in one call.
"""
e6b_decode_stream(state, stream) = decode(state, stream, length(stream))

@testset "Galileo E6-B constructor" begin
    decoder = GalileoE6BDecoderState(1)
    @test decoder.prn == 1
    @test decoder.data isa GalileoE6BData
    @test decoder.constants.syncro_sequence_length == 1000
    @test decoder.constants.preamble_length == 16
    @test decoder.constants.preamble == 0b1011011101110000
    @test decoder.constants.preamble == 0xB770
    @test get_signal_type(decoder) == GalileoE6B
    @test get_band(decoder) == GNSSSignals.E6()
    @test get_data_frequency(decoder) == 1000GNSSSignals.Hz
    # C/NAV rides on E6-B; E6-C is the dataless pilot, so only E6-B maps to a
    # decoder. Checked by method existence rather than by constructing the
    # signal, which would build its spreading-code table.
    @test hasmethod(GNSSDecoderState, Tuple{GalileoE6B,Int})
    @test !hasmethod(GNSSDecoderState, Tuple{GalileoE6C,Int})
    # No ephemeris of its own — corrections only.
    @test !is_decoding_completed_for_positioning(decoder)
    @test !is_sat_healthy(decoder)
end

@testset "Galileo E6-B validity interval table" begin
    # ICD Table 23, with 15 reserved.
    @test GNSSDecoder.e6b_validity_interval(0) == 5
    @test GNSSDecoder.e6b_validity_interval(5) == 60
    @test GNSSDecoder.e6b_validity_interval(10) == 300
    @test GNSSDecoder.e6b_validity_interval(14) == 3600
    @test isnothing(GNSSDecoder.e6b_validity_interval(15))
end

@testset "Galileo E6-B page layout matches HAS SIS ICD Annex C" begin
    bits = e6b_annex_c_page_bits()
    @test length(bits) == 486
    # The ICD's own page carries the ICD's own checksum: CRC-24Q over all 486
    # bits (462 message + 24 CRC) must vanish. This is the strongest available
    # check on the CRC scope, since the checksum is not one we computed.
    @test crc24q(bits) == 0
    # ...and the 462-bit scope reproduces the broadcast checksum exactly.
    stored = UInt32(0)
    for b in bits[463:486]
        stored = (stored << 1) | UInt32(b)
    end
    @test crc24q(bits[1:462]) == stored
    # Reserved field is all ones in this page (ICD Annex C).
    @test all(bits[1:14])

    # Re-encode the ICD's page and drive it through the production decode path:
    # sync → 123×8 deinterleave → Viterbi → CRC → header → encoded-page octets.
    decoder = GalileoE6BDecoderState(1)
    decoder = e6b_decode_stream(decoder, e6b_symbol_stream([bits]))
    @test decoder.cache.page_counter[] == 1
    # MS = 15, so one page cannot complete the message; it must be held.
    groups = decoder.cache.page_groups
    @test haskey(groups, HAS_ANNEX_C_HEADER.message_id)
    group = groups[HAS_ANNEX_C_HEADER.message_id]
    @test group.message_type == HAS_ANNEX_C_HEADER.message_type
    @test group.message_size == HAS_ANNEX_C_HEADER.message_size
    @test length(group.page_ids) == 1
    @test group.page_ids[1] == HAS_ANNEX_C_HEADER.page_id
    @test group.octets[1, :] == hex2bytes(HAS_ANNEX_C_ENCODED_PAGE_HEX)
    # HAS status 0 = test mode: decoded, but not reported as usable.
    @test decoder.raw_data.HAS_status == GNSSDecoder.has_test_mode
    @test !is_sat_healthy(decoder)
    @test isnothing(decoder.data.message)
end

@testset "Galileo E6-B FEC round-trip" begin
    # Encoding a page and decoding it back through the production Viterbi path
    # must recover the exact 486 info bits — confirms the 123×8 deinterleave, the
    # G2 sign flip, and the AFF3CT Viterbi configuration are mutually consistent.
    decoder = GalileoE6BDecoderState(1)
    rng = MersenneTwister(20260821)
    for _ = 1:3
        octets = rand(rng, UInt8, 53)
        header = e6b_header(;
            has_status = 1,
            message_type = 1,
            message_id = 7,
            message_size = 4,
            page_id = 200,
        )
        bits = e6b_page_bits(header, octets)
        on_air = GNSSDecoder.interleave(galileo_conv_encode(bits), 123, 8)
        recovered = GNSSDecoder.galileo_viterbi(
            decoder.cache.viterbi,
            on_air,
            GNSSDecoder.E6B_INTERLEAVER_COLUMNS,
            GNSSDecoder.UInt512,
        )
        expected = GNSSDecoder.UInt512(0)
        for b in bits
            expected = (expected << 1) | GNSSDecoder.UInt512(b)
        end
        @test recovered == expected
    end
end

@testset "Galileo E6-B HPVRS reassembly matches HAS SIS ICD Annex D" begin
    # The step between "k valid pages collected" and "a HAS message": the
    # Reed-Solomon erasure decode. Checked against the ICD's own statement of the
    # reassembled message octets, independently of how they are later parsed —
    # so a reassembly error cannot hide behind a compensating parse error.
    for example in (HAS_EXAMPLE_1, HAS_EXAMPLE_2)
        cache = GNSSDecoder.GalileoE6BCache()
        group = nothing
        for (index, page_id) in enumerate(example.page_ids)
            octets = hex2bytes(example.encoded_pages_hex[index])
            group = GNSSDecoder.e6b_collect_page!(
                cache,
                example.message_id,
                example.message_type,
                example.message_size,
                page_id,
            ) do row
                row .= octets
            end
        end
        # Only the last page completes the group.
        @test !isnothing(group)
        @test length(group.page_ids) == example.message_size
        @test GNSSDecoder.e6b_reassemble_message(group) ==
              has_octets(example.message_octets_hex...)
    end
end

@testset "Galileo E6-B HAS message matches HAS SIS ICD Annex D example 1" begin
    decoder = GalileoE6BDecoderState(1)
    pages = e6b_example_pages(HAS_EXAMPLE_1)
    decoder = e6b_decode_stream(decoder, e6b_symbol_stream(pages))
    @test decoder.cache.page_counter[] == 15
    # A completed message clears its page group.
    @test isempty(decoder.cache.page_groups)
    @test decoder.raw_data.HAS_status == GNSSDecoder.has_operational_mode
    @test is_sat_healthy(decoder)
    # Corrections are published immediately (nothing left to corroborate).
    @test decoder.data == decoder.raw_data

    message = decoder.data.message
    @test !isnothing(message)
    @test message.message_id == 15
    @test message.message_type == 1
    @test message.message_size == 15
    # ICD Annex D, message header: TOH 0, flags [1 1 0 0 1 1], Mask ID 0,
    # IOD Set ID 11.
    @test message.TOH == 0
    @test message.mask_id == 0
    @test message.IOD_set_id == 11
    @test !isnothing(message.mask)
    @test !isnothing(message.orbit_corrections)
    @test isnothing(message.clock_corrections)
    @test isnothing(message.clock_subset_corrections)
    @test !isnothing(message.code_biases)
    @test !isnothing(message.phase_biases)

    # ICD Annex D, mask: N_sys = 2; GPS (GNSS ID 0) with 31 satellites, signals
    # {0, 7}, Cell Mask present, Nav Message 0; Galileo (GNSS ID 2) with 22
    # satellites, signals {1, 4, 7, 13}, no Cell Mask, Nav Message 0.
    mask = message.mask
    @test mask.mask_id == 0
    @test length(mask.satellite_masks) == 2
    gps, galileo = mask.satellite_masks
    @test gps.GNSS_ID == 0
    @test gps.SVIDs == vcat(1:10, 12:32)
    @test gps.signal_indices == [0, 7]
    @test !isnothing(gps.cell_mask)
    @test size(gps.cell_mask) == (31, 2)
    @test gps.nav_message_index == 0
    @test galileo.GNSS_ID == 2
    @test galileo.SVIDs ==
          [1, 2, 3, 4, 5, 7, 8, 9, 11, 12, 13, 15, 19, 21, 24, 25, 26, 27, 30, 31, 33, 36]
    @test galileo.signal_indices == [1, 4, 7, 13]
    @test isnothing(galileo.cell_mask)
    @test galileo.nav_message_index == 0
    @test sum(length(m.SVIDs) for m in mask.satellite_masks) == 53
    # The mask is remembered under its Mask ID for later messages.
    @test decoder.data.masks[0] == mask

    # ICD Annex D, orbit corrections: Validity Interval Index 10 = 300 s, one
    # entry per masked satellite, GPS first.
    orbit = message.orbit_corrections
    @test orbit.validity_interval == 300
    @test orbit.mask_id == 0 && orbit.IOD_set_id == 11 && orbit.TOH == 0
    @test length(orbit.corrections) == 53
    let c = orbit.corrections[1]     # G01: 96, 0.0500, 0.4160, 0.2960
        @test (c.GNSS_ID, c.SVID, c.IOD_ref) == (0, 1, 96)
        @test c.δ_radial ≈ 0.05
        @test c.δ_in_track ≈ 0.416
        @test c.δ_cross_track ≈ 0.296
    end
    let c = orbit.corrections[2]     # G02: the "data not available" sentinels
        @test (c.GNSS_ID, c.SVID, c.IOD_ref) == (0, 2, 0)
        @test isnothing(c.δ_radial)
        @test isnothing(c.δ_in_track)
        @test isnothing(c.δ_cross_track)
    end
    let c = orbit.corrections[32]    # E01: 18, -0.0825, 0.4480, -0.3760
        @test (c.GNSS_ID, c.SVID, c.IOD_ref) == (2, 1, 18)
        @test c.δ_radial ≈ -0.0825
        @test c.δ_in_track ≈ 0.448
        @test c.δ_cross_track ≈ -0.376
    end
    let c = orbit.corrections[53]    # E36: 18, -0.1500, -0.0240, -0.0720
        @test (c.GNSS_ID, c.SVID, c.IOD_ref) == (2, 36, 18)
        @test c.δ_radial ≈ -0.15
        @test c.δ_in_track ≈ -0.024
        @test c.δ_cross_track ≈ -0.072
    end
    # Every Galileo entry in this example references IODnav 15, 17 or 18.
    @test all(c.IOD_ref in (15, 17, 18) for c in orbit.corrections if c.GNSS_ID == 2)

    # ICD Annex D, code biases: Validity Interval Index 14 = 3600 s. GPS carries
    # one or two biases per satellite as the Cell Mask dictates (23 satellites
    # with two, 8 with one); Galileo carries four for each of its 22.
    code = message.code_biases
    @test code.validity_interval == 3600
    @test length(code.corrections) == 142
    @test count(c -> c.GNSS_ID == 0, code.corrections) == 54
    @test count(c -> c.GNSS_ID == 2, code.corrections) == 88
    let c = code.corrections[1]      # G01 signal 0: 3.74 m
        @test (c.GNSS_ID, c.SVID, c.signal_index) == (0, 1, 0)
        @test c.bias ≈ 3.74
    end
    let c = code.corrections[2]      # G01 signal 7: 5.72 m
        @test (c.GNSS_ID, c.SVID, c.signal_index) == (0, 1, 7)
        @test c.bias ≈ 5.72
    end
    let c = code.corrections[3]      # G02: only signal 0 (Cell Mask), -4.38 m
        @test (c.GNSS_ID, c.SVID, c.signal_index) == (0, 2, 0)
        @test c.bias ≈ -4.38
    end
    let e01 = filter(c -> c.GNSS_ID == 2 && c.SVID == 1, code.corrections)
        @test [c.signal_index for c in e01] == [1, 4, 7, 13]
        @test [c.bias for c in e01] ≈ [0.08, 0.14, 0.14, 1.04]
    end

    # ICD Annex D, phase biases: Validity Interval Index 5 = 60 s, and every
    # single value is the -10.24 cycle "data not available" sentinel with a zero
    # Phase Discontinuity Indicator. Exactly the case a decoder that folds the
    # sentinel to 0.0 would silently corrupt.
    phase = message.phase_biases
    @test phase.validity_interval == 60
    @test length(phase.corrections) == 142
    @test all(isnothing(c.bias) for c in phase.corrections)
    @test all(c.phase_discontinuity_indicator == 0 for c in phase.corrections)
    @test [c.signal_index for c in phase.corrections[1:2]] == [0, 7]
end

@testset "Galileo E6-B HAS message matches HAS SIS ICD Annex D example 2" begin
    # Example 2 is a clock-only message that references the mask of another
    # message by Mask ID, so it can only be parsed once example 1's mask is in.
    decoder = GalileoE6BDecoderState(1)
    decoder =
        e6b_decode_stream(decoder, e6b_symbol_stream(e6b_example_pages(HAS_EXAMPLE_1)))
    decoder =
        e6b_decode_stream(decoder, e6b_symbol_stream(e6b_example_pages(HAS_EXAMPLE_2)))
    @test decoder.cache.page_counter[] == 17

    message = decoder.data.message
    @test message.message_id == HAS_EXAMPLE_2.message_id
    @test message.message_size == 2
    # ICD Annex D: TOH 7, flags [0 0 1 0 0 0], Mask ID 0, IOD Set ID 11.
    @test message.TOH == 7
    @test message.mask_id == 0
    @test message.IOD_set_id == 11
    @test isnothing(message.mask)          # referenced, not carried
    @test isnothing(message.orbit_corrections)
    @test !isnothing(message.clock_corrections)
    @test isnothing(message.clock_subset_corrections)
    @test isnothing(message.code_biases)
    @test isnothing(message.phase_biases)
    # The earlier message's blocks survive; only the clock block is replaced.
    @test !isnothing(decoder.data.orbit_corrections)
    @test !isnothing(decoder.data.code_biases)
    @test decoder.data.clock_corrections === message.clock_corrections

    clock = message.clock_corrections
    @test clock.validity_interval == 60     # Validity Interval Index 5
    @test length(clock.corrections) == 53
    # ICD Annex D prints the raw Delta Clock Multiplier fields (2 for GPS, 0 for
    # Galileo) and the *un-multiplied* corrections. Table 29 maps those raw
    # fields to multipliers 3 and 1, so the usable corrections are the printed
    # values times 3 and 1 respectively.
    @test all(c.multiplier == 3 for c in clock.corrections if c.GNSS_ID == 0)
    @test all(c.multiplier == 1 for c in clock.corrections if c.GNSS_ID == 2)
    let c = clock.corrections[1]            # G01: -6.4100 × 3
        @test (c.GNSS_ID, c.SVID) == (0, 1)
        @test c.δ_clock ≈ -6.41 * 3
        @test !c.do_not_use
    end
    let c = clock.corrections[2]            # G02: -10.2400 = data not available
        @test (c.GNSS_ID, c.SVID) == (0, 2)
        @test isnothing(c.δ_clock)
        @test !c.do_not_use
    end
    let c = clock.corrections[3]            # G03: -6.0775 × 3
        @test c.δ_clock ≈ -6.0775 * 3
    end
    let c = clock.corrections[32]           # E01: 0.0800 × 1
        @test (c.GNSS_ID, c.SVID) == (2, 1)
        @test c.δ_clock ≈ 0.08
    end
    let c = clock.corrections[53]           # E36: -0.1025 × 1
        @test (c.GNSS_ID, c.SVID) == (2, 36)
        @test c.δ_clock ≈ -0.1025
    end
end

@testset "Galileo E6-B holds a message until its mask arrives" begin
    # Reverse the order: the clock-only message arrives first, with no mask to
    # size its block. The ICD gives no recovery for this, but masks recur every
    # few messages, so the message is held and re-parsed once one lands.
    decoder = GalileoE6BDecoderState(1)
    decoder =
        e6b_decode_stream(decoder, e6b_symbol_stream(e6b_example_pages(HAS_EXAMPLE_2)))
    @test isnothing(decoder.data.clock_corrections)
    @test isnothing(decoder.data.message)
    @test !isnothing(decoder.cache.pending_message[])

    decoder =
        e6b_decode_stream(decoder, e6b_symbol_stream(e6b_example_pages(HAS_EXAMPLE_1)))
    @test isnothing(decoder.cache.pending_message[])
    # Both messages are now applied: the mask message's blocks and the held
    # clock block.
    @test !isnothing(decoder.data.orbit_corrections)
    @test !isnothing(decoder.data.clock_corrections)
    @test length(decoder.data.clock_corrections.corrections) == 53
    @test decoder.data.clock_corrections.corrections[1].δ_clock ≈ -6.41 * 3

    # ...and the released orphan is merged *underneath* the message that unlocked
    # it, not on top. The orphan completed earlier, so it must not win the
    # single-valued header fields, nor replace a block of the same kind that the
    # newer message carried. Example 1 is the message that just completed here, so
    # it — not the held example 2 — is the one published as "most recent".
    @test decoder.data.message.message_id == HAS_EXAMPLE_1.message_id
    @test decoder.data.message.TOH == 0
    @test decoder.data.message.mask_id == 0
    @test decoder.data.message.IOD_set_id == 11
    @test !isnothing(decoder.data.message.orbit_corrections)
    # Every block each message carried is still published — the orphan's clock
    # block and the newer message's orbit/code/phase blocks do not collide.
    @test !isnothing(decoder.data.code_biases)
    @test !isnothing(decoder.data.phase_biases)
    # Arrival order aside, the outcome matches feeding them the other way round.
    forward = GalileoE6BDecoderState(1)
    forward =
        e6b_decode_stream(forward, e6b_symbol_stream(e6b_example_pages(HAS_EXAMPLE_1)))
    @test decoder.data.orbit_corrections == forward.data.orbit_corrections
    @test decoder.data.code_biases == forward.data.code_biases
    @test decoder.data.message == forward.data.message
end

"""
Pack a bit vector into whole octets, zero-padding the tail — the shape
`parse_has_message` reads a reassembled HAS message in.
"""
function e6b_message_octets(bits::Vector{Bool})
    padded = vcat(bits, falses(mod(-length(bits), 8)))
    map(1:8:length(padded)) do i
        reduce((acc, b) -> (acc << 1) | UInt8(b), padded[i:(i+7)]; init = UInt8(0))
    end
end

"""
The 32 MT1 header bits (HAS SIS ICD Table 12).
"""
function e6b_mt1_header_bits(;
    TOH = 0,
    mask = false,
    orbit = false,
    clock_full = false,
    clock_subset = false,
    code_bias = false,
    phase_bias = false,
    mask_id = 0,
    IOD_set_id = 0,
)
    bits = Bool[]
    galileo_push_field!(bits, TOH, 12)
    append!(bits, (mask, orbit, clock_full, clock_subset, code_bias, phase_bias))
    galileo_push_field!(bits, 0, 4)          # reserved
    galileo_push_field!(bits, mask_id, 5)
    galileo_push_field!(bits, IOD_set_id, 5)
    @assert length(bits) == 32
    bits
end

"""
One constellation's entry in a Mask block (ICD Table 16), with no Cell Mask.
"""
function e6b_system_mask_bits(; GNSS_ID, satellite_mask, signal_mask, nav_message = 0)
    bits = Bool[]
    galileo_push_field!(bits, GNSS_ID, 4)
    galileo_push_field!(bits, satellite_mask, 40)
    galileo_push_field!(bits, signal_mask, 16)
    push!(bits, false)                   # Cell Mask Availability Flag
    galileo_push_field!(bits, nav_message, 3)
    bits
end

@testset "Galileo E6-B rejects the reserved Nsys and Nsys_sub values" begin
    # One GPS satellite (index 0 = PRN 1) on one signal, for the positive controls.
    one_gps_system = e6b_system_mask_bits(;
        GNSS_ID = 0,
        satellite_mask = UInt64(1) << 39,
        signal_mask = UInt16(1) << 15,
    )

    @testset "Nsys = 0 is reserved" begin
        # ICD §5.2.1: Nsys assumes "values from 1 to 15 (value \"0\" is Reserved)".
        # An accepted zero would cache a mask correcting nothing under a real Mask
        # ID, which would then satisfy later Mask Flag = 0 messages instead of
        # holding them until the genuine mask arrives.
        bits = vcat(
            e6b_mt1_header_bits(; mask = true, mask_id = 7),
            galileo_push_field!(Bool[], 0, 4),      # Nsys = 0
            galileo_push_field!(Bool[], 0, 6),      # reserved
        )
        @test isnothing(
            GNSSDecoder.parse_has_message(e6b_message_octets(bits), 1, 1, 1, nothing),
        )
        # Positive control: the same message with Nsys = 1 parses, so the
        # rejection above is the zero and not the framing.
        bits = vcat(
            e6b_mt1_header_bits(; mask = true, mask_id = 7),
            galileo_push_field!(Bool[], 1, 4),
            one_gps_system,
            galileo_push_field!(Bool[], 0, 6),
        )
        message = GNSSDecoder.parse_has_message(e6b_message_octets(bits), 1, 1, 1, nothing)
        @test message isa GalileoHASMessage
        @test message.mask.mask_id == 7
        @test message.mask.satellite_masks[1].SVIDs == [1]
    end

    @testset "Nsys_sub = 0 ends the Clock Subset block" begin
        # The ICD does not reserve Nsys_sub = 0 (Table 32), but a subset block
        # correcting no constellation is vacuous and would otherwise publish an
        # empty `clock_subset_corrections` over a good one.
        bits = vcat(
            e6b_mt1_header_bits(; mask = true, clock_subset = true, mask_id = 3),
            galileo_push_field!(Bool[], 1, 4),
            one_gps_system,
            galileo_push_field!(Bool[], 0, 6),
            galileo_push_field!(Bool[], 5, 4),      # Validity Interval Index
            galileo_push_field!(Bool[], 0, 4),      # Nsys_sub = 0
        )
        message = GNSSDecoder.parse_has_message(e6b_message_octets(bits), 1, 1, 1, nothing)
        # The mask parsed before the bad block survives; the block itself does not.
        @test message isa GalileoHASMessage
        @test !isnothing(message.mask)
        @test isnothing(message.clock_subset_corrections)
        # Positive control: Nsys_sub = 1 for the masked GPS satellite decodes.
        bits = vcat(
            e6b_mt1_header_bits(; mask = true, clock_subset = true, mask_id = 3),
            galileo_push_field!(Bool[], 1, 4),
            one_gps_system,
            galileo_push_field!(Bool[], 0, 6),
            galileo_push_field!(Bool[], 5, 4),      # Validity Interval Index 5 = 60 s
            galileo_push_field!(Bool[], 1, 4),      # Nsys_sub = 1
            galileo_push_field!(Bool[], 0, 4),      # GNSS ID 0 = GPS
            galileo_push_field!(Bool[], 2, 2),      # Delta Clock Multiplier field 2 → ×3
            Bool[true],                         # satellite submask: the one satellite
            galileo_push_field!(Bool[], 400, 13),   # Delta Clock C0 = 400 × 0.0025 m
        )
        message = GNSSDecoder.parse_has_message(e6b_message_octets(bits), 1, 1, 1, nothing)
        subset = message.clock_subset_corrections
        @test !isnothing(subset)
        @test subset.validity_interval == 60
        @test length(subset.corrections) == 1
        @test subset.corrections[1].SVID == 1
        @test subset.corrections[1].multiplier == 3
        @test subset.corrections[1].δ_clock ≈ 400 * 0.0025 * 3
    end
end

@testset "Galileo E6-B arms the framework symbol counter" begin
    # `num_bits_after_valid_syncro_sequence` is documented as symbols since the
    # last valid sync, `nothing` meaning "not yet synchronised". A locked E6-B
    # decoder must not keep reporting `nothing` — that is how a consumer tells a
    # tracking decoder from a dead one.
    decoder = GalileoE6BDecoderState(1)
    @test isnothing(decoder.num_bits_after_valid_syncro_sequence)
    decoder =
        e6b_decode_stream(decoder, e6b_symbol_stream(e6b_example_pages(HAS_EXAMPLE_1)))
    # Re-armed to the trailing next-page sync that `drain_after_sync!` keeps.
    @test decoder.num_bits_after_valid_syncro_sequence == decoder.constants.preamble_length
    # A page whose contents are unusable still framed correctly, so it re-arms too.
    dummy = e6b_page_bits(0xAF3BC3, hex2bytes(HAS_EXAMPLE_1.encoded_pages_hex[1]))
    decoder = e6b_decode_stream(decoder, e6b_symbol_stream([dummy]))
    @test decoder.num_bits_after_valid_syncro_sequence == decoder.constants.preamble_length
end

@testset "Galileo E6-B expires a held message with its page groups" begin
    # The orphan slot has the same 150 s deadline as an incomplete page group:
    # past it, the corrections it carries have outlived their validity intervals,
    # so a late mask must not resurrect them.
    dummy = e6b_page_bits(0xAF3BC3, hex2bytes(HAS_EXAMPLE_1.encoded_pages_hex[1]))
    decoder = GalileoE6BDecoderState(1)
    decoder =
        e6b_decode_stream(decoder, e6b_symbol_stream(e6b_example_pages(HAS_EXAMPLE_2)))
    held = decoder.cache.pending_message[]
    @test !isnothing(held)
    # The slot carries its own Message Type now, rather than assuming MT1.
    @test held.message_type == 1
    @test held.message_id == HAS_EXAMPLE_2.message_id
    @test held.message_size == HAS_EXAMPLE_2.message_size

    decoder = e6b_decode_stream(decoder, e6b_symbol_stream(fill(dummy, 150)))
    @test isnothing(decoder.cache.pending_message[])
    # The mask now arrives, and finds nothing to unlock.
    decoder =
        e6b_decode_stream(decoder, e6b_symbol_stream(e6b_example_pages(HAS_EXAMPLE_1)))
    @test isnothing(decoder.data.clock_corrections)
    @test !isnothing(decoder.data.orbit_corrections)
end

@testset "Galileo E6-B cache equality sees in-flight state" begin
    # A cache holding most of a message must not compare equal to an empty one:
    # the two behave completely differently on the next page.
    empty_state = GalileoE6BDecoderState(1)
    partial = e6b_decode_stream(
        GalileoE6BDecoderState(1),
        e6b_symbol_stream(e6b_example_pages(HAS_EXAMPLE_1)[1:14]),
    )
    @test length(partial.cache.page_groups[HAS_EXAMPLE_1.message_id].page_ids) == 14
    @test partial.cache != empty_state.cache
    @test partial != empty_state

    # ...and, the other way round, two decoders fed the same pages must compare
    # equal. `GalileoHASPageGroup` is mutable, so the default `==` is reference
    # equality: without its own definition this holds only for the empty store,
    # and any decoder with a message in flight would differ from its own twin.
    twin = e6b_decode_stream(
        GalileoE6BDecoderState(1),
        e6b_symbol_stream(e6b_example_pages(HAS_EXAMPLE_1)[1:14]),
    )
    @test partial.cache.page_groups == twin.cache.page_groups
    @test partial.cache == twin.cache
    @test partial == twin
end

@testset "Galileo E6-B page acceptance rules" begin
    example = HAS_EXAMPLE_1
    octets = hex2bytes(example.encoded_pages_hex[1])

    @testset "dummy pages are discarded" begin
        decoder = GalileoE6BDecoderState(1)
        dummy = e6b_page_bits(0xAF3BC3, octets)
        decoder = e6b_decode_stream(decoder, e6b_symbol_stream([dummy, dummy]))
        # Received and counted (they are seconds of signal), but never stored.
        @test decoder.cache.page_counter[] == 2
        @test isempty(decoder.cache.page_groups)
        @test isnothing(decoder.raw_data.HAS_status)
    end

    @testset "reserved page ID and message type are discarded" begin
        for header in (
            e6b_header(;
                has_status = 1,
                message_type = 1,
                message_id = 3,
                message_size = 4,
                page_id = 0,
            ),
            e6b_header(;
                has_status = 1,
                message_type = 2,
                message_id = 3,
                message_size = 4,
                page_id = 40,
            ),
            # Page IDs message_size+1 … 32 are the message's zero padding and are
            # never transmitted (ICD §6.3).
            e6b_header(;
                has_status = 1,
                message_type = 1,
                message_id = 3,
                message_size = 4,
                page_id = 30,
            ),
        )
            decoder = GalileoE6BDecoderState(1)
            decoder = e6b_decode_stream(
                decoder,
                e6b_symbol_stream([e6b_page_bits(header, octets)]),
            )
            @test decoder.cache.page_counter[] == 1
            @test isempty(decoder.cache.page_groups)
            # The status is still read from a structurally valid page.
            @test decoder.raw_data.HAS_status == GNSSDecoder.has_operational_mode
        end
    end

    @testset "reserved HAS status is discarded" begin
        decoder = GalileoE6BDecoderState(1)
        header = e6b_header(;
            has_status = 2,
            message_type = 1,
            message_id = 3,
            message_size = 4,
            page_id = 40,
        )
        decoder =
            e6b_decode_stream(decoder, e6b_symbol_stream([e6b_page_bits(header, octets)]))
        @test decoder.raw_data.HAS_status == GNSSDecoder.has_status_reserved
        @test isempty(decoder.cache.page_groups)
        @test !is_sat_healthy(decoder)
    end

    @testset "duplicate page IDs are ignored" begin
        decoder = GalileoE6BDecoderState(1)
        header = e6b_header(;
            has_status = 1,
            message_type = 1,
            message_id = 3,
            message_size = 4,
            page_id = 40,
        )
        page = e6b_page_bits(header, octets)
        decoder = e6b_decode_stream(decoder, e6b_symbol_stream([page, page, page]))
        @test length(decoder.cache.page_groups[3].page_ids) == 1
    end

    @testset "a new message size restarts the group" begin
        decoder = GalileoE6BDecoderState(1)
        first_page = e6b_page_bits(
            e6b_header(;
                has_status = 1,
                message_type = 1,
                message_id = 3,
                message_size = 4,
                page_id = 40,
            ),
            octets,
        )
        second_page = e6b_page_bits(
            e6b_header(;
                has_status = 1,
                message_type = 1,
                message_id = 3,
                message_size = 5,
                page_id = 41,
            ),
            octets,
        )
        decoder = e6b_decode_stream(decoder, e6b_symbol_stream([first_page, second_page]))
        group = decoder.cache.page_groups[3]
        @test group.message_size == 5
        @test length(group.page_ids) == 1
        @test group.page_ids[1] == 41
    end

    @testset "incomplete messages time out after 150 s" begin
        decoder = GalileoE6BDecoderState(1)
        page = e6b_page_bits(
            e6b_header(;
                has_status = 1,
                message_type = 1,
                message_id = 3,
                message_size = 3,
                page_id = 40,
            ),
            octets,
        )
        dummy = e6b_page_bits(0xAF3BC3, octets)
        decoder = e6b_decode_stream(decoder, e6b_symbol_stream([page]))
        @test haskey(decoder.cache.page_groups, 3)
        # 149 further pages (seconds) keep it alive; the 150th expires it.
        decoder = e6b_decode_stream(decoder, e6b_symbol_stream(fill(dummy, 149)))
        @test haskey(decoder.cache.page_groups, 3)
        decoder = e6b_decode_stream(decoder, e6b_symbol_stream([dummy]))
        @test isempty(decoder.cache.page_groups)
    end
end

@testset "Galileo E6-B HAS status 'do not use' discards everything" begin
    decoder = GalileoE6BDecoderState(1)
    decoder =
        e6b_decode_stream(decoder, e6b_symbol_stream(e6b_example_pages(HAS_EXAMPLE_1)))
    @test !isnothing(decoder.data.orbit_corrections)
    # ICD Table 9: "Users shall stop using HAS from all satellites and discard
    # previously received messages."
    stop = e6b_page_bits(
        e6b_header(;
            has_status = 3,
            message_type = 1,
            message_id = 3,
            message_size = 4,
            page_id = 40,
        ),
        hex2bytes(HAS_EXAMPLE_1.encoded_pages_hex[1]),
    )
    decoder = e6b_decode_stream(decoder, e6b_symbol_stream([stop]))
    @test decoder.data.HAS_status == GNSSDecoder.has_do_not_use
    @test !is_sat_healthy(decoder)
    @test isnothing(decoder.data.orbit_corrections)
    @test isnothing(decoder.data.code_biases)
    @test isnothing(decoder.data.message)
    @test isnothing(decoder.data.masks)
    @test isempty(decoder.cache.page_groups)
end

@testset "Galileo E6-B resolves the 180-degree ambiguity" begin
    pages = e6b_example_pages(HAS_EXAMPLE_1)
    upright = GalileoE6BDecoderState(1)
    upright = e6b_decode_stream(upright, e6b_symbol_stream(pages))
    inverted = GalileoE6BDecoderState(1)
    inverted = e6b_decode_stream(inverted, -e6b_symbol_stream(pages))
    @test !upright.is_shifted_by_180_degrees
    @test inverted.is_shifted_by_180_degrees
    @test inverted.data.message == upright.data.message
    @test inverted.data.orbit_corrections == upright.data.orbit_corrections
end

@testset "Galileo E6-B reset_decoder_state" begin
    decoder = GalileoE6BDecoderState(1)
    decoder =
        e6b_decode_stream(decoder, e6b_symbol_stream(e6b_example_pages(HAS_EXAMPLE_1)))
    # Leave one page of a second message in flight.
    partial = e6b_page_bits(
        e6b_header(;
            has_status = 1,
            message_type = 1,
            message_id = 9,
            message_size = 4,
            page_id = 77,
        ),
        hex2bytes(HAS_EXAMPLE_1.encoded_pages_hex[1]),
    )
    decoder = e6b_decode_stream(decoder, e6b_symbol_stream([partial]))
    @test haskey(decoder.cache.page_groups, 9)

    decoder = reset_decoder_state(decoder)
    # In-flight pages go: their 150 s timeout is counted in received pages, so
    # they would otherwise look fresh forever across an outage.
    @test isempty(decoder.cache.page_groups)
    @test isempty(decoder.cache.soft_buffer)
    @test isnothing(decoder.num_bits_after_valid_syncro_sequence)
    @test isnothing(decoder.data.orbit_corrections)
    # Masks survive — they change only when the corrected set does.
    @test haskey(decoder.raw_data.masks, 0)
    @test decoder.raw_data.HAS_status == GNSSDecoder.has_operational_mode

    # And a clock-only message decodes straight away on the surviving mask.
    decoder =
        e6b_decode_stream(decoder, e6b_symbol_stream(e6b_example_pages(HAS_EXAMPLE_2)))
    @test !isnothing(decoder.data.clock_corrections)
end
