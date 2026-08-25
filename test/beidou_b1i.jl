# BeiDou B1I decoder tests: D1 NAV (MEO/IGSO) and D2 NAV (GEO), exercised
# through the public `decode` API on synthesized transmit chains (see
# `dnav_test_utils.jl`). Field values are checked against the ICD scalings
# applied to the raw integers fed into the encoder (BDS-SIS-ICD-B1I-3.0 §5.2.4).

using GNSSSignals: Hz

@testset "BeiDou B1I" begin
    PI = GNSSDecoder.GNSS_PI

    # Raw field integers for the D1 fundamental set (subframes 1-3).
    SOW0 = 345600
    d1_sf1 = (
        SOW = SOW0,
        SatH1 = 0,
        AODC = 1,
        URAI = 2,
        WN = 810,
        t_0c_raw = 5400,       # t_0c = 43200 s
        T_GD1_raw = -52,       # -5.2 ns
        T_GD2_raw = 37,        #  3.7 ns
        α_raw = (10, -20, 30, -40),
        β_raw = (5, -6, 7, -8),
        a_2_raw = -3,
        a_0_raw = 123456,
        a_1_raw = -2345,
        AODE = 4,
    )
    t_0e_raw = 70000           # t_0e = 560000 s; MSB pair = 2, LSB 15 = 4464
    d1_sf2 = (
        SOW = SOW0 + 6,
        Δn_raw = 1234,
        C_uc_raw = -9876,
        M_0_raw = 305419896,
        e_raw = 180150000,
        C_us_raw = 8765,
        C_rc_raw = -1000,
        C_rs_raw = 2000,
        sqrt_A_raw = 2769679000,
        t_0e_msb2 = t_0e_raw >> 15,
    )
    d1_sf3 = (
        SOW = SOW0 + 12,
        t_0e_lsb15 = t_0e_raw & 0x7FFF,
        i_0_raw = 644245094,
        C_ic_raw = -777,
        Ω_dot_raw = -8000,
        C_is_raw = 555,
        i_dot_raw = -99,
        Ω_0_raw = -1234567890,
        ω_raw = 987654321,
    )

    sf1 = dnav_test_encode_subframe(dnav_test_d1_subframe1_content(; d1_sf1...))
    sf2 = dnav_test_encode_subframe(dnav_test_d1_subframe2_content(; d1_sf2...))
    sf3 = dnav_test_encode_subframe(dnav_test_d1_subframe3_content(; d1_sf3...))

    # BCH(15,11,1) is a perfect code and so detects nothing — the decoder
    # therefore requires one broadcast repetition of a dataset before
    # promoting it (see `dnav_confirm_data`). Every promotion test feeds the
    # fundamental set twice.
    #
    # `Δ` shifts the cycle's SOWs. The stream is contiguous, so consecutive
    # subframes must stay exactly 6 s apart for the SOW-continuity screen:
    # a three-subframe cycle repeats at Δ = 18 s, not at the 30 s a real
    # five-subframe D1 frame would take.
    D1_CYCLE_SPAN = 18
    d1_cycle(Δ) = (
        dnav_test_encode_subframe(
            dnav_test_d1_subframe1_content(; merge(d1_sf1, (SOW = SOW0 + Δ,))...),
        ),
        dnav_test_encode_subframe(
            dnav_test_d1_subframe2_content(; merge(d1_sf2, (SOW = SOW0 + Δ + 6,))...),
        ),
        dnav_test_encode_subframe(
            dnav_test_d1_subframe3_content(; merge(d1_sf3, (SOW = SOW0 + Δ + 12,))...),
        ),
    )
    d1_two_cycles = (d1_cycle(0)..., d1_cycle(D1_CYCLE_SPAN)...)
    # SOW of the last subframe in `d1_two_cycles`.
    D1_LAST_SOW = SOW0 + D1_CYCLE_SPAN + 12

    function check_d1_fundamental(data; SOW = D1_LAST_SOW)
        @test data.SOW == SOW
        @test data.sat_h1 === false
        @test data.AODC == 1
        @test data.urai == 2
        @test data.ura == 2.0^(2 / 2 + 1)
        @test data.WN == 810
        @test data.t_0c == 43200
        @test data.T_GD1 ≈ -52 * 1.0e-10
        @test data.T_GD2 ≈ 37 * 1.0e-10
        @test data.α_0 == 10 / 2^30
        @test data.α_1 == -20 / 2^27
        @test data.α_2 == 30 / 2^24
        @test data.α_3 == -40 / 2^24
        @test data.β_0 == 5.0 * 2^11
        @test data.β_1 == -6.0 * 2^14
        @test data.β_2 == 7.0 * 2^16
        @test data.β_3 == -8.0 * 2^16
        @test data.a_2 == -3 / 2.0^66
        @test data.a_0 == 123456 / 2.0^33
        @test data.a_1 == -2345 / 2.0^50
        @test data.AODE == 4
        @test data.t_0e == 560000
        @test data.Δn == 1234 * PI / 2.0^43
        @test data.C_uc == -9876 / 2.0^31
        @test data.M_0 == 305419896 * PI / 2.0^31
        @test data.e == 180150000 / 2.0^33
        @test data.C_us == 8765 / 2.0^31
        @test data.C_rc == -1000 / 2.0^6
        @test data.C_rs == 2000 / 2.0^6
        @test data.sqrt_A == 2769679000 / 2.0^19
        @test data.i_0 == 644245094 * PI / 2.0^31
        @test data.C_ic == -777 / 2.0^31
        @test data.Ω_dot == -8000 * PI / 2.0^43
        @test data.C_is == 555 / 2.0^31
        @test data.i_dot == -99 * PI / 2.0^43
        @test data.Ω_0 == -1234567890 * PI / 2.0^31
        @test data.ω == 987654321 * PI / 2.0^31
    end

    @testset "D1 subframes 1-3 decode and promote (PRN 20)" begin
        state = BeiDouB1IDecoderState(20)
        symbols = dnav_test_soft_symbols(d1_two_cycles...)
        state = decode(state, symbols, length(symbols))
        @test is_decoding_completed_for_positioning(state)
        @test is_sat_healthy(state)
        check_d1_fundamental(state.data)
        @test state.is_shifted_by_180_degrees == false
    end

    @testset "D1 promotion waits for a broadcast repetition" begin
        # One cycle stages the dataset but must not publish it: BCH(15,11,1)
        # cannot detect a two-error word, so repetition voting is the only
        # thing between a mis-correction and `state.data`.
        state = BeiDouB1IDecoderState(20)
        symbols = dnav_test_soft_symbols(d1_cycle(0)...)
        state = decode(state, symbols, length(symbols))
        @test !is_decoding_completed_for_positioning(state)
        @test state.data == BeiDouDNAVData()
    end

    @testset "D1 time counter is anchored to this subframe's SOW epoch" begin
        # BDS-SIS-ICD-B1I-3.0 §5.2.4.3: the SOW stamps the leading edge of the
        # preamble of *this* subframe, so at promotion the counter must span
        # the whole subframe plus the trailing preamble (311 symbols = 6.22 s
        # at D1's 50 sps). GPS LNAV's HOW TOW stamps the *next* subframe and
        # anchors to the preamble alone; using that here would report the
        # time one whole subframe early.
        state = BeiDouB1IDecoderState(20)
        symbols = dnav_test_soft_symbols(d1_two_cycles...)
        state = decode(state, symbols, length(symbols))
        @test state.num_bits_after_valid_syncro_sequence == 311
        @test state.data.SOW + state.num_bits_after_valid_syncro_sequence / 50 ≈
              D1_LAST_SOW + 6.22
    end

    @testset "D1 decode with 180-degree phase shift" begin
        state = BeiDouB1IDecoderState(20)
        symbols = dnav_test_soft_symbols(d1_two_cycles...; polarity = -1.0f0)
        state = decode(state, symbols, length(symbols))
        @test is_decoding_completed_for_positioning(state)
        @test state.is_shifted_by_180_degrees == true
        check_d1_fundamental(state.data)
    end

    @testset "D1 single bit error is corrected by BCH(15,11,1)" begin
        # Flip one transmitted bit inside word 3 of subframe 1 (bit 65 of 300,
        # part of the WN field's BCH block).
        sf1_err = sf1 ⊻ UInt320(1) << (300 - 65)
        state = BeiDouB1IDecoderState(20)
        symbols = dnav_test_soft_symbols(sf1_err, sf2, sf3, d1_cycle(D1_CYCLE_SPAN)...)
        state = decode(state, symbols, length(symbols))
        @test is_decoding_completed_for_positioning(state)
        @test state.data.WN == 810
    end

    @testset "A subframe whose SOW fails the screen contributes nothing" begin
        # BCH(15,11,1) is perfect: `dnav_bch_decode` always returns *some*
        # codeword and can never reject one, so the SOW-continuity screen is
        # this message's only integrity gate. A subframe that fails it must
        # not write its payload either — otherwise a false lock injects
        # ephemeris fields behind a rejected timestamp.
        bad_sf2 = dnav_test_encode_subframe(
            dnav_test_d1_subframe2_content(;
                merge(d1_sf2, (SOW = SOW0 + 1000, M_0_raw = 1))...,
            ),
        )
        state = BeiDouB1IDecoderState(20)
        symbols = dnav_test_soft_symbols(sf1, bad_sf2, sf3)
        state = decode(state, symbols, length(symbols))
        # Neither the rejected subframe's payload nor a promotion built on it.
        @test isnothing(state.raw_data.M_0)
        @test state.data == BeiDouDNAVData()
        # And rejecting must not re-open the gate: the *next* subframe is still
        # screened against the last accepted SOW, two subframes back.
        @test state.raw_data.SOW == SOW0 + 12
    end

    @testset "t_0e halves are only paired across adjacent subframes 2 and 3" begin
        # D1 splits t_0e between subframes 2 and 3 with no IOD to pair them
        # by. If subframe 2 is missed, the stored MSBs belong to an earlier
        # frame and must not be joined to fresh LSBs — that would fabricate a
        # t_0e belonging to neither ephemeris, plausible enough to key the
        # voting dataset.
        state = BeiDouB1IDecoderState(20)
        # Cycle 1 in full, then cycle 2 with subframe 2 replaced by a
        # subframe 1 (so subframe 3 does not follow a subframe 2).
        c2 = d1_cycle(D1_CYCLE_SPAN)
        stand_in_sf1 = dnav_test_encode_subframe(
            dnav_test_d1_subframe1_content(;
                merge(d1_sf1, (SOW = SOW0 + D1_CYCLE_SPAN + 6,))...,
            ),
        )
        symbols = dnav_test_soft_symbols(d1_cycle(0)..., c2[1], stand_in_sf1, c2[3])
        state = decode(state, symbols, length(symbols))
        # Subframe 2's MSBs were dropped when it was skipped, so no t_0e is
        # assembled and nothing is promoted.
        @test isnothing(state.raw_data.t_0e)
        @test state.data == BeiDouDNAVData()
    end

    @testset "No sync on corrupted trailing preamble" begin
        state = BeiDouB1IDecoderState(20)
        symbols = dnav_test_soft_symbols(sf1)
        symbols[301] = -symbols[301]  # first bit of the trailing preamble
        symbols[303] = -symbols[303]
        state = decode(state, symbols, length(symbols))
        @test isnothing(state.raw_data.WN)
        @test state.data == BeiDouDNAVData()
    end

    @testset "D1 subframe 4/5 pages: almanac, health, time offsets" begin
        state = BeiDouB1IDecoderState(20)

        alm = (
            sqrt_A_raw = 10_460_000,
            a_1_raw = -12,
            a_0_raw = 345,
            Ω_0_raw = -4_000_000,
            e_raw = 6789,
            δi_raw = -321,
            t_oa_raw = 147,
            Ω_dot_raw = -654,
            ω_raw = 2_222_222,
            M_0_raw = -3_333_333,
        )
        # Subframe 4 page 1 → SV 1 almanac; announces AmEpID = 11.
        sf4_p1 = dnav_test_encode_subframe(
            dnav_test_d1_almanac_page_content(;
                FraID = 4,
                SOW = D1_LAST_SOW + 6,
                Pnum = 1,
                am_field = 0b11,
                alm...,
            ),
        )
        # Subframe 5 page 8: health codes for SV 20-30, WN_a, t_oa.
        health_codes = UInt16[k == 3 ? 0x100 : 0x000 for k = 1:11]  # SV 22 unhealthy
        sf5_p8 = dnav_test_encode_subframe(
            dnav_test_d1_subframe5_page8_content(;
                SOW = D1_LAST_SOW + 12,
                health_codes,
                WN_a = 55,
                t_oa_raw = 147,
            ),
        )
        sf5_p9 = dnav_test_encode_subframe(
            dnav_test_d1_subframe5_page9_content(;
                SOW = D1_LAST_SOW + 18,
                A_0GPS_raw = -15,
                A_1GPS_raw = 25,
                A_0Gal_raw = -35,
                A_1Gal_raw = 45,
                A_0GLO_raw = -55,
                A_1GLO_raw = 65,
            ),
        )
        sf5_p10 = dnav_test_encode_subframe(
            dnav_test_d1_subframe5_page10_content(;
                SOW = D1_LAST_SOW + 24,
                Δt_LS_raw = 4,
                Δt_LSF_raw = 5,
                WN_LSF = 123,
                A_0UTC_raw = -777777,
                A_1UTC_raw = 8888,
                DN = 6,
            ),
        )
        # Subframe 5 page 11 with AmID = 01 → expanded almanac of SV 31.
        sf5_p11 = dnav_test_encode_subframe(
            dnav_test_d1_almanac_page_content(;
                FraID = 5,
                SOW = D1_LAST_SOW + 30,
                Pnum = 11,
                am_field = 0b01,
                alm...,
            ),
        )
        # One contiguous stream. `dnav_test_soft_symbols` appends the *next*
        # subframe's 11-bit preamble, so feeding the fundamental set and the
        # pages in two calls would duplicate those 11 symbols and put every
        # later subframe 311 rather than 300 symbols from its predecessor —
        # which the SOW-continuity screen correctly rejects.
        symbols = dnav_test_soft_symbols(
            d1_two_cycles...,
            sf4_p1,
            sf5_p8,
            sf5_p9,
            sf5_p10,
            sf5_p11,
        )
        state = decode(state, symbols, length(symbols))

        for sv_id in (1, 31)
            @test !isnothing(state.data.almanac) && haskey(state.data.almanac, sv_id)
            entry = state.data.almanac[sv_id]
            @test entry.sqrt_A == 10_460_000 / 2.0^11
            @test entry.a_1 == -12 / 2.0^38
            @test entry.a_0 == 345 / 2.0^20
            @test entry.Ω_0 == -4_000_000 * PI / 2.0^23
            @test entry.e == 6789 / 2.0^21
            @test entry.δi == -321 * PI / 2.0^19
            @test entry.t_oa == 147 * 2^12
            @test entry.Ω_dot == -654 * PI / 2.0^38
            @test entry.ω == 2_222_222 * PI / 2.0^23
            @test entry.M_0 == -3_333_333 * PI / 2.0^23
        end
        @test state.data.AmEpID == 0b11
        @test state.data.health[22] == 0x100
        @test state.data.health[20] == 0x000
        @test state.data.health[30] == 0x000
        @test state.data.WN_a == 55
        @test state.data.t_oa == 147 * 2^12
        @test state.data.A_0GPS == -15 * 1.0e-10
        @test state.data.A_1GPS == 25 * 1.0e-10
        @test state.data.A_0Gal == -35 * 1.0e-10
        @test state.data.A_1Gal == 45 * 1.0e-10
        @test state.data.A_0GLO == -55 * 1.0e-10
        @test state.data.A_1GLO == 65 * 1.0e-10
        @test state.data.Δt_LS == 4
        @test state.data.Δt_LSF == 5
        @test state.data.WN_LSF == 123
        @test state.data.A_0UTC == -777777 / 2.0^30
        @test state.data.A_1UTC == 8888 / 2.0^50
        @test state.data.DN == 6
    end

    @testset "D2 subframe 1 pages 1-10 (GEO PRN 4)" begin
        # Raw values for the D2 fundamental set, split across pages per
        # Figures 5-14-1..-10.
        WN = 812
        t_0c_raw = 5400
        a_0_raw = UInt32(123456) & 0xFFFFFF          # 24 bits
        a_1_raw = mod(-2345, 1 << 22)                # 22 bits, two's complement
        a_2_raw = mod(-3, 1 << 11)                   # 11 bits
        C_uc_raw = mod(-9876, 1 << 18)
        M_0_raw = UInt32(305419896)
        C_us_raw = mod(8765, 1 << 18)
        e_raw = UInt32(180150000)
        sqrt_A_raw = UInt32(2769679000)
        C_ic_raw = mod(-777, 1 << 18)
        C_is_raw = mod(555, 1 << 18)
        t_0e_raw17 = 70000
        i_0_raw = UInt32(644245094)
        C_rc_raw = mod(-1000, 1 << 18)
        C_rs_raw = mod(2000, 1 << 18)
        Ω_dot_raw = mod(-8000, 1 << 24)
        Ω_0_raw = mod(-1234567890, Int64(1) << 32)
        ω_raw = UInt32(987654321)
        i_dot_raw = mod(-99, 1 << 14)
        Δn_raw = mod(1234, 1 << 16)
        α_raw = (mod(10, 256), mod(-20, 256), mod(30, 256), mod(-40, 256))
        β_raw = (mod(5, 256), mod(-6, 256), mod(7, 256), mod(-8, 256))

        SOW_D2 = 259200
        tails = [
            # Page 1: SatH1, AODC, URAI, WN, toc, TGD1, TGD2
            Tuple{Int,Int}[
                (0, 1),
                (1, 5),
                (2, 4),
                (WN, 13),
                (t_0c_raw, 17),
                (mod(-52, 1 << 10), 10),
                (mod(37, 1 << 10), 10),
            ],
            # Page 2: ionosphere
            Tuple{Int,Int}[
                (α_raw[1], 8),
                (α_raw[2], 8),
                (α_raw[3], 8),
                (α_raw[4], 8),
                (β_raw[1], 8),
                (β_raw[2], 8),
                (β_raw[3], 8),
                (β_raw[4], 8),
            ],
            # Page 3: Rev(6+22+10 = 38), a0(24), a1 MSB 4 (Figure 5-14-3)
            Tuple{Int,Int}[(0, 38), (Int(a_0_raw), 24), (Int(a_1_raw >> 18), 4)],
            # Page 4: a1 mid 6 + LSB 12, a2 10+1, AODE, Δn, Cuc MSB 14
            Tuple{Int,Int}[
                (Int(a_1_raw >> 12) & 0x3F, 6),
                (Int(a_1_raw) & 0xFFF, 12),
                (Int(a_2_raw >> 1), 10),
                (Int(a_2_raw) & 1, 1),
                (4, 5),
                (Δn_raw, 16),
                (Int(C_uc_raw >> 4), 14),
            ],
            # Page 5: Cuc LSB 4, M0, Cus, e MSB 10
            Tuple{Int,Int}[
                (Int(C_uc_raw) & 0xF, 4),
                (Int(M_0_raw), 32),
                (C_us_raw, 18),
                (Int(e_raw >> 22), 10),
            ],
            # Page 6: e mid 6 + LSB 16, √A, Cic MSB 10
            Tuple{Int,Int}[
                (Int(e_raw >> 16) & 0x3F, 6),
                (Int(e_raw) & 0xFFFF, 16),
                (Int(sqrt_A_raw), 32),
                (Int(C_ic_raw >> 8), 10),
            ],
            # Page 7: Cic mid 6 + LSB 2, Cis, toe 2+15, i0 MSB 7 + mid 14
            Tuple{Int,Int}[
                (Int(C_ic_raw >> 2) & 0x3F, 6),
                (Int(C_ic_raw) & 0x3, 2),
                (C_is_raw, 18),
                (t_0e_raw17 >> 15, 2),
                (t_0e_raw17 & 0x7FFF, 15),
                (Int(i_0_raw >> 25), 7),
                (Int(i_0_raw >> 11) & 0x3FFF, 14),
            ],
            # Page 8: i0 mid 6 + LSB 5, Crc 17+1, Crs, Ω̇ MSB 3 + mid 16
            Tuple{Int,Int}[
                (Int(i_0_raw >> 5) & 0x3F, 6),
                (Int(i_0_raw) & 0x1F, 5),
                (Int(C_rc_raw >> 1), 17),
                (Int(C_rc_raw) & 1, 1),
                (C_rs_raw, 18),
                (Int(Ω_dot_raw >> 21), 3),
                (Int(Ω_dot_raw >> 5) & 0xFFFF, 16),
            ],
            # Page 9: Ω̇ LSB 5, Ω0 1+22+9, ω MSB 13 + mid 14
            Tuple{Int,Int}[
                (Int(Ω_dot_raw) & 0x1F, 5),
                (Int(Ω_0_raw >> 31), 1),
                (Int(Ω_0_raw >> 9) & 0x3FFFFF, 22),
                (Int(Ω_0_raw) & 0x1FF, 9),
                (Int(ω_raw >> 19), 13),
                (Int(ω_raw >> 5) & 0x3FFF, 14),
            ],
            # Page 10: ω LSB 5, IDOT 1+13
            Tuple{Int,Int}[
                (Int(ω_raw) & 0x1F, 5),
                (Int(i_dot_raw >> 13), 1),
                (Int(i_dot_raw) & 0x1FFF, 13),
            ],
        ]

        state = BeiDouB1IDecoderState(4)   # PRN 4 is a GEO satellite ⇒ D2
        @test GNSSDecoder.is_beidou_geo(4)
        # One 3-second frame per page: subframe 1 carries the page, subframes
        # 2-5 (integrity/differential, out of scope) are fed as zero-filled
        # filler so the SOW spacing between subframe-1 decodes is genuine.
        # Two 30 s broadcast cycles: as for D1, the decoder publishes a
        # dataset only once the broadcast has repeated it (BCH(15,11,1)
        # detects nothing). A single cycle is exercised separately below.
        D2_CYCLE_SPAN = 3 * length(tails)
        function feed_d2_cycle(state, Δ)
            for (page, tail) in enumerate(tails)
                sow = SOW_D2 + Δ + 3 * (page - 1)
                frame = UInt320[dnav_test_encode_subframe(
                    dnav_test_d2_page_content(page, tail; SOW = sow),
                ),]
                for fra_id = 2:5
                    filler = dnav_test_content(
                        (DNAV_TEST_PREAMBLE, 11),
                        (0, 4),
                        (fra_id, 3),
                        (sow, 20),
                        (0, 186),
                    )
                    push!(frame, dnav_test_encode_subframe(filler))
                end
                symbols = dnav_test_soft_symbols(frame...)
                # Drop the helper's trailing preamble except after the last
                # frame: the next frame's own subframe 1 supplies it.
                if page < length(tails)
                    symbols = symbols[1:1500]
                end
                state = decode(state, symbols, length(symbols))
            end
            state
        end
        state = feed_d2_cycle(state, 0)
        # One cycle stages but does not publish.
        @test state.data == BeiDouDNAVData()
        state = feed_d2_cycle(state, D2_CYCLE_SPAN)

        @test is_decoding_completed_for_positioning(state)
        @test is_sat_healthy(state)
        data = state.data
        @test data.SOW == SOW_D2 + D2_CYCLE_SPAN + 27
        # D2 reads a SOW in subframe 1 only (§5.3.3.1(2)), so subframes 2-5
        # must leave the counter running rather than re-anchor it to a SOW
        # up to four subframes old. The stream ends four filler subframes
        # after the promoting subframe 1, so the counter reads
        # 311 + 4*300 and the implied time is SOW + 1511/500 s.
        @test state.num_bits_after_valid_syncro_sequence == 311 + 4 * 300
        @test data.SOW + state.num_bits_after_valid_syncro_sequence / 500 ≈
              SOW_D2 + D2_CYCLE_SPAN + 27 + 3.022
        @test data.WN == WN
        @test data.AODC == 1
        @test data.urai == 2
        @test data.t_0c == 43200
        @test data.T_GD1 ≈ -52 * 1.0e-10
        @test data.T_GD2 ≈ 37 * 1.0e-10
        @test data.α_0 == 10 / 2^30
        @test data.α_3 == -40 / 2^24
        @test data.β_1 == -6.0 * 2^14
        @test data.a_0 == 123456 / 2.0^33
        @test data.a_1 == -2345 / 2.0^50
        @test data.a_2 == -3 / 2.0^66
        @test data.AODE == 4
        @test data.Δn == 1234 * PI / 2.0^43
        @test data.C_uc == -9876 / 2.0^31
        @test data.M_0 == 305419896 * PI / 2.0^31
        @test data.e == 180150000 / 2.0^33
        @test data.C_us == 8765 / 2.0^31
        @test data.sqrt_A == 2769679000 / 2.0^19
        @test data.C_ic == -777 / 2.0^31
        @test data.C_is == 555 / 2.0^31
        @test data.t_0e == 560000
        @test data.i_0 == 644245094 * PI / 2.0^31
        @test data.C_rc == -1000 / 2.0^6
        @test data.C_rs == 2000 / 2.0^6
        @test data.Ω_dot == -8000 * PI / 2.0^43
        @test data.Ω_0 == -1234567890 * PI / 2.0^31
        @test data.ω == 987654321 * PI / 2.0^31
        @test data.i_dot == -99 * PI / 2.0^43
    end

    @testset "Signal metadata and constructors" begin
        state = BeiDouB1IDecoderState(20)
        @test get_signal_type(state) == BeiDouB1I
        @test GNSSDecoderState(BeiDouB1I(), 20) == state
        @test get_data_frequency(state) == get_data_frequency(BeiDouB1I)
        @test get_time_system_name(state) == "BeiDou Time"
        @test get_constellation_name(state) == get_constellation_name(BeiDouB1I)

        # The state-level accessor reports the rate the state actually
        # demodulates: GEO PRNs broadcast D2 at 500 sps, ten times the
        # type-level D1 answer (BDS-SIS-ICD-B1I-3.0 §5.1.1).
        geo_state = BeiDouB1IDecoderState(3)
        @test get_data_frequency(geo_state) == 10 * get_data_frequency(BeiDouB1I)
        @test get_data_frequency(BeiDouB3IDecoderState(59)) ==
              10 * get_data_frequency(BeiDouB3I)
        # Stated absolutely too, so the accessor cannot drift with whatever the
        # signal type happens to answer.
        @test get_data_frequency(state) == 50Hz
        @test get_data_frequency(geo_state) == 500Hz

        # Every BeiDou constructor rejects an out-of-range PRN rather than
        # building a decoder that can never sync.
        for ctor in (
            BeiDouB1IDecoderState,
            BeiDouB3IDecoderState,
            BeiDouB1CDecoderState,
            BeiDouB2aDecoderState,
            BeiDouB2bDecoderState,
        )
            @test_throws ArgumentError ctor(0)
            @test_throws ArgumentError ctor(64)
        end

        # One `is_sat_healthy` method serves both signals: B1I and B3I read the
        # same SatH1 out of the same message, unlike GPS L5I/L2CM which select
        # different health bits from a shared container.
        @test which(is_sat_healthy, (typeof(state),)) ===
              which(is_sat_healthy, (typeof(BeiDouB3IDecoderState(25)),))
    end

    @testset "reset_decoder_state clears SOW and D2 pages" begin
        state = BeiDouB1IDecoderState(20)
        symbols = dnav_test_soft_symbols(d1_two_cycles...)
        state = decode(state, symbols, length(symbols))
        state = reset_decoder_state(state)
        @test isnothing(state.raw_data.SOW)
        @test state.raw_data.WN == 810              # ephemeris survives the reset
        @test state.data == BeiDouDNAVData()
        @test isnothing(state.num_bits_after_valid_syncro_sequence)
        @test isempty(state.cache.d2_pages)
        # And decoding resumes cleanly after the reset.
        state = decode(state, symbols, length(symbols))
        @test is_decoding_completed_for_positioning(state)
        @test state.data.WN == 810
    end
end

@testset "BeiDou D1/D2 data equality is structural" begin
    # `BeiDouDNAVData` carries the `almanac` and `health` `Dictionary` fields, so
    # the default struct `==` would be reference equality and two decoders fed
    # the same subframes would never compare equal once either was populated.
    mk(M_0) = BeiDouDNAVData(;
        WN = 810,
        almanac = Dictionary([7], [BeiDouDNAVAlmanac(; M_0)]),
        health = Dictionary([7], [0x0000]),
    )
    @test mk(0.1) == mk(0.1)
    @test mk(0.1) != mk(0.2)
    @test BeiDouDNAVData() == BeiDouDNAVData()

    # The voting wrapper needs its own method rather than inheriting that one:
    # Julia's default struct `==` is `===`, which compares fields with `===`
    # instead of dispatching to their `==`. Without it the cache's `old_data`
    # vector — and so the whole state — stays reference-compared.
    voted(M_0) = [GNSSDecoder.VotedBeiDouDNAVData(3, mk(M_0))]
    @test voted(0.1) == voted(0.1)
    @test voted(0.1) != voted(0.2)
    @test GNSSDecoder.VotedBeiDouDNAVData(4, mk(0.1)) !=
          GNSSDecoder.VotedBeiDouDNAVData(3, mk(0.1))

    # ...so that it reaches the decoder state, which compares field by field.
    twin(M_0) =
        GNSSDecoderState(BeiDouB1IDecoderState(6); raw_data = mk(M_0), data = mk(M_0))
    @test twin(0.1) == twin(0.1)
    @test twin(0.1) != twin(0.2)

    # Subframe voting is deliberately *not* this comparison: `dnav_compare_data`
    # weighs the ephemeris fields the vote is about, and the almanac is not one
    # of them, so these two still count as the same candidate.
    @test GNSSDecoder.dnav_compare_data(mk(0.1), mk(0.2))
end
