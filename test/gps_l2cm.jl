using Test
using Random
using GNSSDecoder
using GNSSSignals
using GNSSDecoder: crc24q

# GPS L2C broadcasts the same CNAV message as GPS L5I (IS-GPS-200N §30 ≡
# IS-GPS-705J §20.3.3), so these tests reuse the shared CNAV transmit-chain
# helpers (CNAVTestEncoder, fec_encode_soft, build_mt10) from
# `cnav_test_utils.jl`. They focus on what is L2C-specific: the GPSL2CM
# dispatch, that the shared decoder runs through a GPSL2CM decoder state, and
# that `is_sat_healthy` reads the L2 health bit (53) rather than the L5 bit (54).

# ---------------------------------------------------------------------------
# Spirent fixture.
#
# `test/data/gps_l2c_prn25_nav_bits.bin` holds the first 26 consecutive
# 300-bit CNAV messages of PRN 25 from a Spirent L2C simulation, in the
# simulator's pre-FEC block layout (38 bytes per message, 300 bits MSB-first,
# 4 padding LSBs in the last byte) — the same GSS-CNAVDATA format as the L5I
# fixture. They were extracted from the GSS-CNAVDATA output (`l2c_nav_data.cnv`,
# L2C pre-FEC bits: a 16-byte header followed by 1178-byte epochs of 31
# satellite channels × 38 bytes in round-robin order, PRN 1..31) by
# demultiplexing satellite channel 25. The 26 messages cover one full
# 24-message broadcast cycle, i.e. every CNAV message type the simulator emits
# (10-15 and 30-37). The test FEC-encodes them continuously (the fixture is
# pre-FEC) and the golden field values below come from Spirent's own reference
# decode (`l2c_nav_data.txt`). Note L2C messages are 12 s apart, so the message
# TOW count increments by 2 between messages (vs 1 on L5I's 6 s messages),
# while the ×6 scaling of the count is unchanged.
# ---------------------------------------------------------------------------

"""
Load the L2C fixture's 26 messages as 300-bit `Vector{Bool}`s.
"""
function load_l2c_fixture_messages()
    raw = read(joinpath(@__DIR__, "data", "gps_l2c_prn25_nav_bits.bin"))
    @assert length(raw) == 26 * 38
    map(0:25) do i
        block = raw[(38i+1):(38i+38)]
        bits = [(byte >> (7 - j)) & 0x01 == 0x01 for byte in block for j = 0:7]
        bits[1:300]
    end
end

@testset "GPS L2C (CNAV)" begin
    PI = 3.1415926535898

    @testset "Constructor" begin
        state_a = GPSL2CMDecoderState(7)
        state_b = GNSSDecoderState(GPSL2CM(), 7)
        @test state_a.prn == 7
        @test state_a == state_b
        @test state_a.data == GPSCNAVData()          # shared CNAV container
        @test state_a.constants isa GNSSDecoder.GPSL2CMConstants
        @test isnothing(state_a.num_bits_after_valid_syncro_sequence)
        @test !state_a.is_shifted_by_180_degrees
    end

    @testset "Shares the L5I CNAV decode path" begin
        # One message followed by the next message's preamble, prefixed with an
        # odd number of idle symbols (arbitrary symbol-pair phase) — identical
        # exercise to the L5I synthetic test, here through a GPSL2CM state.
        enc = CNAVTestEncoder()
        prefix = fec_encode_soft(enc, falses(40))[1:79]  # odd-length prefix
        stream = vcat(
            prefix,
            fec_encode_soft(enc, [build_mt10(); build_mt10(; tow_count = 1235)]),
        )

        state = decode(GPSL2CMDecoderState(9), stream, length(stream))
        d = state.raw_data
        @test d.last_message_id == 10
        @test d.TOW == 1234 * 6          # TOW count × 6 on L2C too (IS-GPS-200N §30.3.3)
        @test d.WN == 2345
        @test d.M_0 ≈ 12345 * 2.0^-32 * PI
        @test d.e ≈ 1000 * 2.0^-34
        @test d.ω ≈ -54321 * 2.0^-32 * PI
        @test !state.is_shifted_by_180_degrees
    end

    @testset "is_sat_healthy reads the L2 health bit (not L5)" begin
        # MT10 with L2 OK (bit 53 = 0) but L5 bad (bit 54 = 1). The exact same
        # symbol stream must read healthy on L2C and unhealthy on L5I — this is
        # the one decode-level difference between the two signals.
        enc = CNAVTestEncoder()
        msgs = [
            build_mt10(; l2_health = false, l5_health = true)
            build_mt10(; tow_count = 1235, l2_health = false, l5_health = true)
        ]
        stream = fec_encode_soft(enc, msgs)

        state_l2c = decode(GPSL2CMDecoderState(9), stream, length(stream))
        @test state_l2c.raw_data.l2_health == false
        @test state_l2c.raw_data.l5_health == true
        # Health is reported from validated `data`; the synthetic MT10 alone is
        # not a complete positioning set, so promote raw_data to compare the
        # bit selection directly.
        @test is_sat_healthy(
            GNSSDecoder.GNSSDecoderState(state_l2c; data = state_l2c.raw_data),
        )

        state_l5i = decode(GPSL5IDecoderState(9), stream, length(stream))
        @test !is_sat_healthy(
            GNSSDecoder.GNSSDecoderState(state_l5i; data = state_l5i.raw_data),
        )

        # And the converse: L2 bad, L5 OK ⇒ L2C unhealthy, L5I healthy.
        enc2 = CNAVTestEncoder()
        msgs2 = [
            build_mt10(; l2_health = true, l5_health = false)
            build_mt10(; tow_count = 1235, l2_health = true, l5_health = false)
        ]
        stream2 = fec_encode_soft(enc2, msgs2)
        s2_l2c = decode(GPSL2CMDecoderState(9), stream2, length(stream2))
        @test !is_sat_healthy(GNSSDecoder.GNSSDecoderState(s2_l2c; data = s2_l2c.raw_data))
        s2_l5i = decode(GPSL5IDecoderState(9), stream2, length(stream2))
        @test is_sat_healthy(GNSSDecoder.GNSSDecoderState(s2_l5i; data = s2_l5i.raw_data))
    end

    @testset "Phase-inverted stream decodes with is_shifted_by_180_degrees" begin
        enc = CNAVTestEncoder()
        stream =
            fec_encode_soft(enc, [falses(20); build_mt10(); build_mt10(; tow_count = 1235)])
        state = decode(GPSL2CMDecoderState(9), -stream, length(stream))
        @test state.is_shifted_by_180_degrees
        @test state.raw_data.last_message_id == 10
        @test state.raw_data.WN == 2345
    end

    @testset "reset_decoder_state clears in-flight state, keeps decoded data" begin
        enc = CNAVTestEncoder()
        stream = fec_encode_soft(enc, [build_mt10(); build_mt10(; tow_count = 1235)])
        state = decode(GPSL2CMDecoderState(9), stream, length(stream))
        @test state.raw_data.WN == 2345

        state = reset_decoder_state(state)
        @test isempty(state.cache.soft_buffer)
        @test isnothing(state.raw_data.TOW)
        @test isnothing(state.num_bits_after_valid_syncro_sequence)
        @test !state.is_shifted_by_180_degrees
        @test state.raw_data.WN == 2345  # ephemeris preserved
    end

    # --- Spirent-derived recording fixture (committed; always runs) ---------

    """
    Assert the decoded `state` matches Spirent's L2C reference decode (PRN 25).
    """
    function assert_l2c_spirent_golden(state)
        d = state.data
        # Header of the last decoded message (#25 of the recording, message
        # type 10, message TOW count 43250). L2C messages are 12 s apart, so the
        # count steps by 2 per message (43202, 43204, …); the ×6 scaling holds.
        @test d.last_message_id == 10
        @test d.TOW == 43250 * 6
        @test d.alert_flag == false

        # Message type 10 — ephemeris 1 (semicircle fields stored in radians, × π).
        @test d.WN == 2106
        @test d.l1_health == false
        @test d.l2_health == false
        @test d.l5_health == false
        @test d.t_op == 259200
        @test d.ura_ed_index == 0
        @test d.t_0e == 264600
        @test d.ΔA ≈ 922.048828125
        @test d.A_dot ≈ 0.0
        @test d.Δn_0 ≈ 3.7252334550430533e-9 * π
        @test d.Δn_0_dot ≈ 1.000000082740371e-11 * π
        @test d.M_0 ≈ -0.2567902144510299 * π
        @test d.e ≈ 0.099999999976716936
        @test d.ω ≈ -0.3913338293787092 * π
        @test d.integrity_status_flag == false
        @test d.l2c_phasing == false

        # Message type 11 — ephemeris 2.
        @test d.Ω_0 ≈ 0.51429651258513331 * π
        @test d.i_0 ≈ 0.34116418519988656 * π
        @test d.ΔΩ_dot ≈ 3.7252334550430533e-9 * π
        @test d.i_dot ≈ 9.9987573776161298e-11 * π
        @test d.C_is ≈ 4.0000304579734802e-6
        @test d.C_ic ≈ 2.9997900128364563e-6
        @test d.C_rs ≈ 6.0
        @test d.C_rc ≈ 5.0
        @test d.C_us ≈ 1.9995495676994324e-6
        @test d.C_uc ≈ 1.0002404451370239e-6

        # Message types 30-37 — shared clock block.
        @test d.ura_ned0_index == 0
        @test d.ura_ned1_index == 0
        @test d.ura_ned2_index == 0
        @test d.t_0c == 264600
        @test d.a_f0 ≈ 6.4874766394495964e-6
        @test d.a_f1 ≈ 1.0324008314910316e-9
        @test d.a_f2 ≈ 4.4322184811207421e-16

        # Message type 30 — group delay + ionosphere.
        @test d.T_GD ≈ 9.8953023552894592e-10
        @test d.ISC_L1CA ≈ 0.0
        @test d.ISC_L2C ≈ 0.0
        @test d.ISC_L5I5 ≈ 0.0
        @test d.ISC_L5Q5 ≈ 0.0
        @test d.α_0 ≈ 4.6566128730773926e-9
        @test d.α_1 ≈ 1.4901161193847656e-8
        @test d.α_2 ≈ -5.9604644775390625e-8
        @test d.α_3 ≈ -5.9604644775390625e-8
        @test d.β_0 ≈ 79872.0
        @test d.β_1 ≈ 65536.0
        @test d.β_2 ≈ -65536.0
        @test d.β_3 ≈ -393216.0
        @test d.WN_op == 58

        # Message type 32 — EOP. Spirent's reference decoder predates the
        # 2⁻²³ ΔUT_GPS definition (IS-GPS-200N Table 30-VII) and prints the
        # 31-bit field as ΔUT1 with the old 2⁻²⁴ scale (-0.25315952301025391 s);
        # our decoder follows the current ICD, i.e. exactly twice that.
        @test d.t_EOP == 259200
        @test d.PM_X ≈ 0.10612583160400391
        @test d.PM_X_dot ≈ 0.0013251304626464844
        @test d.PM_Y ≈ 0.4459381103515625
        @test d.PM_Y_dot ≈ 0.00040197372436523438
        @test d.ΔUT_GPS ≈ 2 * -0.25315952301025391
        @test d.ΔUT_GPS_dot ≈ -0.00046199560165405273  # rate scale unchanged

        # Message type 33 — UTC.
        @test d.A0_UTC ≈ 9.5364521257579327e-7
        @test d.A1_UTC ≈ 1.8185453143360064e-12
        @test d.A2_UTC ≈ 0.0
        @test d.Δt_LS == 18
        @test d.t_ot == 507904
        @test d.WN_ot == 2106
        @test d.WN_LSF == 2104
        @test d.DN == 2
        @test d.Δt_LSF == 18

        # Message type 35 — GGTO.
        @test d.t_GGTO == 259200
        @test d.WN_GGTO == 2106
        @test d.GNSS_ID == 0
        @test d.A0_GGTO ≈ 0.0
        @test d.A1_GGTO ≈ 0.0
        @test d.A2_GGTO ≈ 0.0

        # Reduced almanacs: message type 31 broadcasts PRNs 1-4, message
        # type 12 PRNs 5-11 (semicircle fields × π).
        @test !isnothing(d.reduced_almanacs) && length(d.reduced_almanacs) == 11
        ra1 = d.reduced_almanacs[1]
        @test ra1.WN_a == 2106
        @test ra1.t_oa == 507904
        @test ra1.δA ≈ 1024.0
        @test ra1.Ω_0 ≈ 0.296875 * π
        @test ra1.Φ_0 ≈ -0.140625 * π
        @test !ra1.l1_health && !ra1.l2_health && !ra1.l5_health
        ra5 = d.reduced_almanacs[5]
        @test ra5.δA ≈ 1024.0
        @test ra5.Ω_0 ≈ -0.765625 * π
        @test ra5.Φ_0 ≈ -0.234375 * π
        ra11 = d.reduced_almanacs[11]
        @test ra11.Ω_0 ≈ -0.125 * π
        @test ra11.Φ_0 ≈ 0.171875 * π

        # Midi almanac (message type 37): one SV per message; the recording's
        # first cycle carries PRN 26.
        @test !isnothing(d.midi_almanacs) && length(d.midi_almanacs) == 1
        ma = d.midi_almanacs[26]
        @test ma.WN_a == 2106
        @test ma.t_oa == 507904
        @test ma.e ≈ 0.0
        @test ma.δi ≈ 0.00555419921875 * π
        @test ma.Ω_dot ≈ 0.0
        # Spirent's dump prints the truncated "5153"; the transmitted 17-bit
        # field is 82459 = 5153.6875.
        @test ma.sqrt_A ≈ 5153.6875
        @test ma.Ω_0 ≈ 0.6190185546875 * π
        @test ma.ω ≈ 0.308013916015625 * π
        @test ma.M_0 ≈ 0.364288330078125 * π
        @test ma.a_f0 ≈ 0.0
        @test ma.a_f1 ≈ 0.0
        @test !ma.l1_health && !ma.l2_health && !ma.l5_health

        # Differential corrections: message type 13 carries CDC packets for
        # PRNs 1-6, message type 14 EDC packets for PRNs 1-2, and message
        # type 34 one CDC (PRN 7) + EDC (PRN 3) pair.
        @test !isnothing(d.clock_corrections) && length(d.clock_corrections) == 7
        cdc1 = d.clock_corrections[1]
        @test cdc1.t_op_D == 259200
        @test cdc1.t_OD == 259200
        @test cdc1.δa_f0 ≈ 0.0
        @test cdc1.δa_f1 ≈ 0.0
        @test cdc1.UDRA_index == 0
        cdc7 = d.clock_corrections[7]  # from message type 34
        @test cdc7.t_op_D == 259200
        @test cdc7.t_OD == 259200
        @test !isnothing(d.ephemeris_corrections) && length(d.ephemeris_corrections) == 3
        edc3 = d.ephemeris_corrections[3]  # from message type 34
        @test edc3.t_op_D == 259200
        @test edc3.t_OD == 259200
        @test edc3.Δα ≈ 0.0
        @test edc3.Δβ ≈ 0.0
        @test edc3.Δγ ≈ 0.0
        @test edc3.Δi ≈ 0.0
        @test edc3.ΔΩ ≈ 0.0
        @test edc3.ΔA ≈ 0.0
        @test edc3.UDRA_dot_index == 0

        # Text messages (the simulator broadcasts "Spirent Communications L2C
        # test text message." across pages; the first cycle carries page 1).
        @test d.text_mt15 == "Spirent Communications L2C te"
        @test d.text_page_mt15 == 1
        @test d.text_mt36 == "Spirent Communicat"
        @test d.text_page_mt36 == 1

        # No message type 40 in the recording.
        @test isnothing(d.ism)

        @test GNSSDecoder.is_decoding_completed_for_positioning(d)
        @test is_sat_healthy(state)
    end

    @testset "Spirent recording fixture (PRN 25)" begin
        messages = load_l2c_fixture_messages()
        @test all(crc24q(m) == 0 for m in messages)
        stream = fec_encode_soft(CNAVTestEncoder(), reduce(vcat, messages))
        @test length(stream) == 26 * 600

        # Of the 26 encoded messages the last cannot complete (its sync
        # window needs the next message's first 16 symbols), so 25 decode.
        state = decode(GPSL2CMDecoderState(25), stream, length(stream))
        @test !state.is_shifted_by_180_degrees
        assert_l2c_spirent_golden(state)

        # Polarity-inverted stream must decode identically with the
        # 180°-flip flag set.
        state_inv = decode(GPSL2CMDecoderState(25), -stream, length(stream))
        @test state_inv.is_shifted_by_180_degrees
        @test state_inv.data == state.data

        # Mid-stream acquisition at an arbitrary (odd) symbol offset: the
        # decoder must slide to the next message boundary and decode the
        # remaining messages, losing only the partially received first one.
        offset = 351
        state_mid =
            decode(GPSL2CMDecoderState(25), stream[(offset+1):end], length(stream) - offset)
        assert_l2c_spirent_golden(state_mid)

        # Noisy soft symbols: moderate Gaussian noise on the ±1 LLRs must
        # not cost any message at this SNR.
        noisy = stream .+ 0.4f0 .* randn(MersenneTwister(7), Float32, length(stream))
        state_noisy = decode(GPSL2CMDecoderState(25), noisy, length(noisy))
        assert_l2c_spirent_golden(state_noisy)
    end

    @testset "decode_once stops at the first complete positioning set" begin
        messages = load_l2c_fixture_messages()
        stream = fec_encode_soft(CNAVTestEncoder(), reduce(vcat, messages))
        state = decode(GPSL2CMDecoderState(25), stream, length(stream); decode_once = true)
        # Message types 10, 11, 30 are the first three of the recording; the
        # positioning set is complete after message 3 (type 30) and decoding
        # stops there (message 4 is type 15 — its text must not have been
        # parsed into the validated `data`).
        @test GNSSDecoder.is_decoding_completed_for_positioning(state.data)
        @test state.data.last_message_id == 30
        @test isnothing(state.data.text_mt15)
        @test state.raw_data.text_mt15 == "Spirent Communications L2C te"
    end
end
