# Galileo E5a F/NAV decoder tests.
#
# Ground truth is a Spirent constellation-simulator capture (PRN 21), shipped as
# `test/data/galileo_e5a_fnav_pages.bin`: 152 decoded F/NAV pages, each the
# 244-bit page (238 info bits + 6 zero tail bits) packed MSB-first into 31 bytes
# (4 trailing pad bits). The companion Spirent text dump (with per-field decoded
# values and verified CRCs) was used to transcribe the expected values asserted
# below; every numeric field here was cross-checked against that dump.
#
# The capture provides decoded pages, not raw E5a-I symbols, so the test
# *re-encodes* each page through the transmit FEC chain (rate-1/2 K=7 NSC
# convolutional code with G2 inverted, 61×8 block interleave, 12-symbol sync
# prefix) to synthesise the on-air soft-symbol stream, then drives the full
# `decode` path (sync → deinterleave → Viterbi → CRC → parse). The re-encoder is
# verified to round-trip against the decoder, and — critically — the decoded
# field values are compared against the *independent* Spirent ground truth, so a
# shared encode/decode error cannot mask a field-layout mistake. CRC validity is
# also checked directly against Spirent's (independently computed) checksums.

const GALILEO_E5A_SYNC = (1, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0)
const GALILEO_E5A_FNAV_PAGES_PATH = joinpath(@__DIR__, "data", "galileo_e5a_fnav_pages.bin")

# Read the 152 packed pages back into 244-bit page values (UInt256).
function read_e5a_fnav_pages(path)
    raw = read(path)
    npages = length(raw) ÷ 31
    pages = Vector{UInt256}(undef, npages)
    for p = 1:npages
        v = UInt256(0)
        for b = 1:31
            v = (v << 8) | UInt256(raw[(p-1)*31+b])
        end
        pages[p] = v >> 4   # 248 packed bits → 244 page bits
    end
    pages
end

# The 238 information bits (page type + data + CRC) of a 244-bit page.
e5a_info_bits(page244::UInt256) = page244 >> 6

# Transmit-side rate-1/2 K=7 NSC convolutional encoder, G1 = 0o171, G2 = 0o133
# with the G2 output inverted (Galileo OS SIS ICD §4.1.4). Takes the 238 info
# bits, appends 6 zero tail bits, and returns 488 ±1 soft symbols (bit 1 → −1).
function e5a_conv_encode(info::Vector{Bool})
    polys = (0o171, 0o133)
    reg = zeros(Bool, 6)
    out = Float32[]
    for b in vcat(info, zeros(Bool, 6))
        taps = vcat([b], reg)               # [current, last-6]; tap j matches poly MSB
        for (k, p) in enumerate(polys)
            acc = false
            for j = 1:7
                ((p >> (7 - j)) & 1 == 1) && (acc ⊻= taps[j])
            end
            sym = (k == 2) ? !acc : acc     # invert G2
            push!(out, sym ? -1.0f0 : 1.0f0)
        end
        reg = vcat([b], reg[1:(end-1)])
    end
    out
end

# One F/NAV page → 500 on-air soft symbols: 12-symbol sync + 488 interleaved
# encoded symbols (8 rows × 61 columns block interleaver).
function e5a_page_to_symbols(page244::UInt256)
    info = e5a_info_bits(page244)
    info_vec = Bool[GNSSDecoder.get_bit(info, 238, k) for k = 1:238]
    sync = Float32[b == 1 ? -1.0f0 : 1.0f0 for b in GALILEO_E5A_SYNC]
    vcat(sync, GNSSDecoder.interleave(e5a_conv_encode(info_vec), 61, 8))
end

# Concatenate the given pages into a soft-symbol stream, plus a trailing sync so
# the final page's window is closed by the next page-sync pattern.
function e5a_symbol_stream(pages)
    sync = Float32[b == 1 ? -1.0f0 : 1.0f0 for b in GALILEO_E5A_SYNC]
    stream = Float32[]
    for page in pages
        append!(stream, e5a_page_to_symbols(page))
    end
    append!(stream, sync)
    stream
end

@testset "Galileo E5a constructor" begin
    decoder = GalileoE5aDecoderState(21)
    @test decoder.prn == 21
    @test decoder.data isa GNSSDecoder.GalileoE5aData
    @test decoder.constants.syncro_sequence_length == 500
    @test decoder.constants.preamble_length == 12
    @test decoder.constants.preamble == 0b101101110000
    @test GNSSDecoder.num_bits_buffered(decoder) == 0
    @test isnothing(decoder.num_bits_after_valid_syncro_sequence)

    # F/NAV rides on the E5a-I (data) component; dispatch from the GNSSSignals
    # system type must match the direct constructor.
    @test GalileoE5aDecoderState(21) == GNSSDecoderState(GNSSSignals.GalileoE5aI(), 21)
end

@testset "Galileo E5a FEC round-trip" begin
    # Encoding a page and decoding it back through the production Viterbi path
    # must recover the exact 238 info bits. Confirms the deinterleave (61×8), the
    # G2 sign flip, and the AFF3CT Viterbi configuration are mutually consistent.
    pages = read_e5a_fnav_pages(GALILEO_E5A_FNAV_PAGES_PATH)
    decoder = GalileoE5aDecoderState(21)
    for idx in (1, 5, 10, 76, 152)
        info = e5a_info_bits(pages[idx])
        info_vec = Bool[GNSSDecoder.get_bit(info, 238, k) for k = 1:238]
        on_air = GNSSDecoder.interleave(e5a_conv_encode(info_vec), 61, 8)
        recovered = GNSSDecoder.galileo_e5a_viterbi(decoder.cache.viterbi_decoder, on_air)
        @test recovered == GNSSDecoder.UInt256(info)
    end
end

@testset "Galileo E5a CRC matches Spirent" begin
    # Every Spirent page is flagged CRC-CORRECT; our CRC-24Q over the 238 info
    # bits must therefore vanish for all 152 pages. Independent of the FEC chain.
    pages = read_e5a_fnav_pages(GALILEO_E5A_FNAV_PAGES_PATH)
    @test length(pages) == 152
    @test all(pages) do page
        info = e5a_info_bits(page)
        crc_bits = Bool[GNSSDecoder.get_bit(info, 238, k) for k = 1:238]
        crc24q(crc_bits) == 0
    end
end

@testset "Galileo E5a ephemeris frame decoding" begin
    pages = read_e5a_fnav_pages(GALILEO_E5A_FNAV_PAGES_PATH)

    # Expected first-frame navigation data (Spirent PRN 21, IODnav 48, TOW 259230).
    # Every value cross-checked against the Spirent text dump (word types 1-4).
    expected = GNSSDecoder.GalileoE5aData(;
        WN = 1082,
        TOW = 259230,
        SVID = 21,
        t_0e = 259200.0,
        M_0 = 2.583259511699057,
        e = 9.999994654208422e-5,
        sqrt_A = 5440.588203430176,
        Ω_0 = 1.4358151884944024,
        i_0 = 0.9885709926241616,
        ω = 3.1416165752262324e-5,
        i_dot = 3.142988060925546e-11,
        Ω_dot = 3.141595145762181e-8,
        Δn = 3.142988060925545e-10,
        C_uc = 1.000240445137024e-6,
        C_us = 2.000480890274048e-6,
        C_rc = 5.0,
        C_rs = 6.0,
        C_ic = 3.0007213354110718e-6,
        C_is = 3.999099135398865e-6,
        SISA_e1_e5a = 25,
        t_0c = 259200.0,
        a_f0 = 9.999802568927407e-5,
        a_f1 = 1.0000036354540498e-9,
        a_f2 = 1.734723475976807e-18,
        IOD_nav1 = 0x0000000000000030,
        IOD_nav2 = 0x0000000000000030,
        IOD_nav3 = 0x0000000000000030,
        IOD_nav4 = 0x0000000000000030,
        num_pages_after_last_TOW = 2,
        num_bits_after_valid_syncro_sequence_after_last_TOW = nothing,
        signal_health_e5a = GNSSDecoder.signal_ok,
        data_validity_status_e5a = GNSSDecoder.navigation_data_valid,
        broadcast_group_delay_e1_e5a = -9.313225746154785e-10,
        a_i0 = 100.0,
        a_i1 = 1.0,
        a_i2 = 0.100006103515625,
        iono_storm_flag_region1 = true,
        iono_storm_flag_region2 = false,
        iono_storm_flag_region3 = true,
        iono_storm_flag_region4 = false,
        iono_storm_flag_region5 = false,
        A_0_utc = 1.000240445137024e-6,
        A_1_utc = 8.881784197001252e-16,
        Δt_LS = 18,
        t_0t = 259200,
        WN_0t = 58,
        WN_LSF = 56,
        DN = 2,
        Δt_LSF = 18,
        A_0G = -2.9103830456733704e-11,
        A_1G = -4.440892098500626e-16,
        t_0G = 918000,
        WN_0G = 63,
    )

    # Feed the first five pages (word types 1-5): word type 4 completes the
    # ephemeris/clock set and triggers validation.
    stream = e5a_symbol_stream(pages[1:5])
    decoder = GalileoE5aDecoderState(21)
    decoder = decode(decoder, stream, length(stream))

    @test decoder.data == expected
    @test decoder.is_shifted_by_180_degrees == false
    @test is_sat_healthy(decoder) == true
    # 1012 = symbols elapsed since the leading sync edge of word type 4 (the page
    # carrying the validated TOW): two further pages (WT4→WT5→trailing sync) of
    # 500 symbols each, less the 488 — i.e. preamble (12) + 2·500.
    @test decoder.num_bits_after_valid_syncro_sequence == 1012

    # Inverting the whole stream (180° Costas ambiguity) decodes identically but
    # flags the polarity flip.
    decoder_inv = GalileoE5aDecoderState(21)
    decoder_inv = decode(decoder_inv, -stream, length(stream))
    @test decoder_inv.data == expected
    @test decoder_inv.is_shifted_by_180_degrees == true
    @test is_sat_healthy(decoder_inv) == true
end

@testset "Galileo E5a almanac decoding" begin
    pages = read_e5a_fnav_pages(GALILEO_E5A_FNAV_PAGES_PATH)

    # Word type 5 (t=240 s) carries SVID-19 (full) + SVID-20 (first half); word
    # type 6 (t=290 s) completes SVID-20 and carries SVID-21 (full). The two are
    # five pages apart in the broadcast schedule; the partial persists in the
    # cache across the intervening WT1-4 pages. WN_a/t_0a are broadcast only in
    # WT5 and shared by all three almanacs, including SVID-21 from WT6.
    stream = e5a_symbol_stream(pages[1:31])
    decoder = GalileoE5aDecoderState(21)
    decoder = decode(decoder, stream, length(stream))
    almanacs = decoder.data.almanacs
    @test !isnothing(almanacs)

    # SVID 19 — full almanac from word type 5.
    @test almanacs[19] == GNSSDecoder.GalileoAlmanac(;
        SVID = 19,
        Δsqrt_A = 0.0,
        e = 0.0,
        ω = 0.0,
        δi = 0.0,
        Ω_0 = 2.815621736164101,
        Ω_dot = 0.0,
        M_0 = 0.9006384700873591,
        a_f0 = 0.0,
        a_f1 = 0.0,
        signal_health_e5a = GNSSDecoder.signal_ok,
        IOD_a = 0,
        WN_a = 2,
        t_0a = 259200,
    )

    # SVID 20 — split across word types 5 and 6 (note the stitched Ω_0).
    @test almanacs[20] == GNSSDecoder.GalileoAlmanac(;
        SVID = 20,
        Δsqrt_A = 0.0,
        e = 0.0,
        ω = 0.0,
        δi = 0.0,
        Ω_0 = 2.815621736164101,
        Ω_dot = 0.0,
        M_0 = 1.6860366334848091,
        a_f0 = 0.0,
        a_f1 = 0.0,
        signal_health_e5a = GNSSDecoder.signal_ok,
        IOD_a = 0,
        WN_a = 2,
        t_0a = 259200,
    )

    # SVID 21 — full almanac from word type 6; WN_a/t_0a inherited from WT5.
    @test almanacs[21] == GNSSDecoder.GalileoAlmanac(;
        SVID = 21,
        Δsqrt_A = 0.0,
        e = 0.0001068115234375,
        ω = 0.0,
        δi = 0.011121360712170923,
        Ω_0 = 1.4358060174609635,
        Ω_dot = 3.1452738704244004e-8,
        M_0 = 2.5832236467994254,
        a_f0 = 9.918212890625e-5,
        a_f1 = 1.000444171950221e-9,
        signal_health_e5a = GNSSDecoder.signal_ok,
        IOD_a = 0,
        WN_a = 2,
        t_0a = 259200,
    )
end

@testset "Galileo E5a almanac WT6 without preceding WT5" begin
    pages = read_e5a_fnav_pages(GALILEO_E5A_FNAV_PAGES_PATH)
    # Broadcast schedule: WT5 at pages 5/15/25, WT6 at pages 10/20/30. Feed pages
    # 1-24 to lock the decoder — page 24 is a WT4, so the SVID-2 chain partial is
    # empty at that point — then append the WT6 at page 30 while omitting its
    # paired WT5 at page 25. This is the mid-stream acquisition / IOD-cutover case:
    # the WT6 carries SVID-3 in full, but its reference epoch (WN_a/t_0a) lived
    # only in the missing WT5. Per the decode-everything policy the record is still
    # stored — nothing decodable is discarded — with the epoch left `nothing`. (The
    # completion of the SVID-2 chain is likewise skipped, so this WT6 stores exactly
    # one almanac.)
    partial_stream = e5a_symbol_stream(vcat(pages[1:24], pages[30]))
    decoder = GalileoE5aDecoderState(21)
    decoder = decode(decoder, partial_stream, length(partial_stream))
    almanacs = decoder.raw_data.almanacs
    @test !isnothing(almanacs)
    partial = only(almanacs)
    @test !isnothing(partial.Δsqrt_A)   # orbit decoded from the WT6
    @test isnothing(partial.WN_a)       # but the reference epoch is absent
    @test isnothing(partial.t_0a)

    # The full run (WT5 then WT6, pages 1-31) instead yields the epoch too.
    full_stream = e5a_symbol_stream(pages[1:31])
    full_decoder = GalileoE5aDecoderState(21)
    full_decoder = decode(full_decoder, full_stream, length(full_stream))
    full = full_decoder.data.almanacs[21]
    @test full.WN_a == 2
    @test full.t_0a == 259200
end

@testset "Galileo E5a almanac epoch back-patch from a later WT5" begin
    pages = read_e5a_fnav_pages(GALILEO_E5A_FNAV_PAGES_PATH)
    # Same lock prefix + lone WT6 (page 30) as above, so SVID-21 lands with its
    # orbit/clock/health decoded but WN_a/t_0a still `nothing`. Then append a WT5
    # (page 25) with a matching IOD_a — but *never* re-send the WT6. Because the
    # reference epoch is shared across every almanac of that IOD_a, the WT5
    # back-fills it into the stranded SVID-21 record, completing an almanac whose
    # only WT6 was seen once. A "store only when complete" policy would have
    # discarded that WT6 and could never reassemble SVID-21 from the WT5 alone.
    stream = e5a_symbol_stream(vcat(pages[1:24], pages[30], pages[25]))
    decoder = GalileoE5aDecoderState(21)
    decoder = decode(decoder, stream, length(stream))
    almanacs = decoder.raw_data.almanacs
    @test !isnothing(almanacs)
    @test haskey(almanacs, 21)
    completed = almanacs[21]
    @test !isnothing(completed.Δsqrt_A)   # orbit still from the one-shot WT6
    @test completed.WN_a == 2             # epoch back-filled by the later WT5
    @test completed.t_0a == 259200
end

@testset "Galileo E5a reset" begin
    pages = read_e5a_fnav_pages(GALILEO_E5A_FNAV_PAGES_PATH)
    stream = e5a_symbol_stream(pages[1:5])
    decoder = GalileoE5aDecoderState(21)
    decoder = decode(decoder, stream, length(stream))
    @test !isnothing(decoder.data.TOW)

    decoder = reset_decoder_state(decoder)
    @test GNSSDecoder.num_bits_buffered(decoder) == 0
    @test isnothing(decoder.raw_data.TOW)
    @test isnothing(decoder.data.TOW)
    @test isnothing(decoder.num_bits_after_valid_syncro_sequence)
    # Ephemeris is retained in raw_data across a reset (fast reacquisition).
    @test decoder.raw_data.sqrt_A == 5440.588203430176
end
