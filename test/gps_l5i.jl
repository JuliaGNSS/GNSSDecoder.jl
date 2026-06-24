using Test
using Random
using GNSSDecoder
using GNSSSignals
using GNSSDecoder: crc24q, gps_l5i_viterbi

# ---------------------------------------------------------------------------
# Reference CNAV transmit chain (test-only).
#
# The rate-1/2, K=7 convolutional encoder (G1 = 0o171, G2 = 0o133,
# IS-GPS-705J §3.3.3.1.1, Figure 3-7) is intentionally re-implemented here so
# the tests catch regressions in either the decoder's Viterbi or the
# message parsing. The encoder runs *continuously* across message boundaries
# — no tail bits, no reset — exactly like the satellite. Symbols are emitted
# as ±1 Float32 (bit 0 ⇒ +1, bit 1 ⇒ -1), the package-wide LLR convention.
# ---------------------------------------------------------------------------

"""
Continuous K=7 rate-1/2 FEC encoder state (6-bit register, [s1..s6], s1 most recent).
"""
mutable struct L5ITestEncoder
    register::UInt8
end
L5ITestEncoder() = L5ITestEncoder(0x00)

"""
Encode one bit; returns the (G1, G2) symbol pair as `Bool`s.
"""
function encode_bit!(enc::L5ITestEncoder, bit::Bool)
    u = UInt8(bit)
    s = enc.register
    s1 = (s >> 5) & 0x01
    s2 = (s >> 4) & 0x01
    s3 = (s >> 3) & 0x01
    s5 = (s >> 1) & 0x01
    s6 = s & 0x01
    y1 = u ⊻ s1 ⊻ s2 ⊻ s3 ⊻ s6  # G1 = 0o171
    y2 = u ⊻ s2 ⊻ s3 ⊻ s5 ⊻ s6  # G2 = 0o133
    enc.register = ((u << 5) | (s >> 1)) & 0x3f
    return (y1 == 0x01, y2 == 0x01)
end

"""
FEC-encode a bit stream into ±1 Float32 soft symbols (bit 0 ⇒ +1, bit 1 ⇒ -1).
"""
function fec_encode_soft(enc::L5ITestEncoder, bits::AbstractVector{Bool})
    soft = Vector{Float32}(undef, 2 * length(bits))
    for (i, b) in enumerate(bits)
        y1, y2 = encode_bit!(enc, b)
        soft[2i-1] = y1 ? -1.0f0 : 1.0f0
        soft[2i] = y2 ? -1.0f0 : 1.0f0
    end
    return soft
end

"""
Write `len` bits of `val` MSB-first into 1-based position `start` of `bits`.
"""
function setbits!(bits::BitVector, start::Int, len::Int, val::Integer)
    mask = (UInt64(1) << len) - UInt64(1)
    v = UInt64(unsigned(Int64(val)) & mask)
    @inbounds for i = 0:(len-1)
        bits[start+i] = ((v >> (len - 1 - i)) & UInt64(1)) == UInt64(1)
    end
    return bits
end

"""
Build a CNAV message type 10 with chosen field values and a valid CRC.
"""
function build_mt10(; prn = 9, tow_count = 1234)
    bits = falses(300)
    setbits!(bits, 1, 8, 0b10001011)   # preamble
    setbits!(bits, 9, 6, prn)
    setbits!(bits, 15, 6, 10)          # message type 10
    setbits!(bits, 21, 17, tow_count)
    setbits!(bits, 38, 1, 0)           # alert
    setbits!(bits, 39, 13, 2345)       # WN
    setbits!(bits, 52, 1, 0)           # L1 health
    setbits!(bits, 53, 1, 0)           # L2 health
    setbits!(bits, 54, 1, 0)           # L5 health
    setbits!(bits, 55, 11, 100)        # t_op (×300 = 30000)
    setbits!(bits, 66, 5, 0)           # URA_ED index
    setbits!(bits, 71, 11, 100)        # t_0e (×300 = 30000)
    setbits!(bits, 82, 26, 0)          # ΔA
    setbits!(bits, 108, 25, 0)         # A_dot
    setbits!(bits, 133, 17, 100)       # Δn_0
    setbits!(bits, 150, 23, 0)         # Δn_0_dot
    setbits!(bits, 173, 33, 12345)     # M_0
    setbits!(bits, 206, 33, 1000)      # e
    setbits!(bits, 239, 33, -54321)    # ω
    setbits!(bits, 272, 1, 1)          # integrity status flag
    setbits!(bits, 273, 1, 0)          # L2C phasing
    setbits!(bits, 277, 24, crc24q(collect(bits[1:276])))
    return collect(bits)
end

# ---------------------------------------------------------------------------
# Spirent fixture.
#
# `test/data/gps_l5i_prn25_nav_bits.bin` holds the first 26 consecutive
# 300-bit CNAV messages of PRN 25 from a Spirent GSS7000 L5 simulation, in
# the simulator's pre-FEC block layout (Spirent FAQ14061: 38 bytes per
# message, 300 bits MSB-first, 4 padding LSBs in the last byte). They were
# extracted from the GSS-CNAVDATA output (`l5_nav_data.cnv`, file type 64 =
# L5C pre-FEC bits, one block per satellite channel per 6-second epoch in
# round-robin order; the full recording carries 31 channels = PRN 1..31) by
# demultiplexing satellite channel 24. The 26 messages cover one full
# 24-message broadcast cycle, i.e. every CNAV message type the simulator
# emits (10-15 and 30-37). The test FEC-encodes them continuously (the
# fixture is pre-FEC, so the encoder above recreates the transmitted symbol
# stream) and golden field values below come from Spirent's own reference
# decode (`l5_nav_data.txt`).
# ---------------------------------------------------------------------------

"""
Load the fixture's 26 messages as 300-bit `Vector{Bool}`s.
"""
function load_l5i_fixture_messages()
    raw = read(joinpath(@__DIR__, "data", "gps_l5i_prn25_nav_bits.bin"))
    @assert length(raw) == 26 * 38
    map(0:25) do i
        block = raw[(38i+1):(38i+38)]
        bits = [(byte >> (7 - j)) & 0x01 == 0x01 for byte in block for j = 0:7]
        bits[1:300]
    end
end

@testset "GPS L5I (CNAV)" begin
    PI = 3.1415926535898

    @testset "Constructor" begin
        state_a = GPSL5IDecoderState(13)
        state_b = GNSSDecoderState(GPSL5I(), 13)
        @test state_a.prn == 13
        @test state_a == state_b
        @test state_a.data == GPSL5IData()
        @test isnothing(state_a.num_bits_after_valid_syncro_sequence)
        @test !state_a.is_shifted_by_180_degrees
    end

    @testset "Window Viterbi round-trips the continuous FEC" begin
        Random.seed!(1234)
        bits = rand(Bool, 308)
        soft = fec_encode_soft(L5ITestEncoder(), bits)
        @test gps_l5i_viterbi(soft) == bits

        # Soft decisions: moderate Gaussian noise must still decode exactly
        # (rate-1/2 K=7 at this SNR has plenty of margin).
        noisy = soft .+ 0.4f0 .* randn(MersenneTwister(42), Float32, length(soft))
        @test gps_l5i_viterbi(noisy) == bits

        # The encoder runs across message boundaries: a window cut from the
        # *middle* of a longer stream (unknown initial register) must decode
        # exactly as well.
        enc = L5ITestEncoder()
        fec_encode_soft(enc, rand(Bool, 100))  # advance the register
        soft_mid = fec_encode_soft(enc, bits)
        @test gps_l5i_viterbi(soft_mid) == bits
    end

    @testset "Synthetic message type 10 round trip" begin
        mt10 = build_mt10()
        @test crc24q(mt10) == 0

        # One message followed by the next message's preamble, with a prefix
        # of an odd number of idle symbols: sync must cope with an arbitrary
        # symbol-pair phase.
        enc = L5ITestEncoder()
        prefix = fec_encode_soft(enc, falses(40))[1:79]  # odd-length prefix
        stream = vcat(prefix, fec_encode_soft(enc, [mt10; build_mt10(; tow_count = 1235)]))

        state = decode(GPSL5IDecoderState(9), stream, length(stream))
        d = state.raw_data
        @test d.last_message_id == 10
        # Only the first message completes (the second one's sync window
        # would need a third message's preamble symbols).
        @test d.TOW == 1234 * 6
        @test d.alert_flag == false
        @test d.WN == 2345
        @test d.l1_health == false
        @test d.l2_health == false
        @test d.l5_health == false
        @test d.t_op == 30_000
        @test d.ura_ed_index == 0
        @test d.t_0e == 30_000
        @test d.ΔA == 0.0
        @test d.A_dot == 0.0
        @test d.Δn_0 ≈ 100 * 2.0^-44 * PI
        @test d.M_0 ≈ 12345 * 2.0^-32 * PI
        @test d.e ≈ 1000 * 2.0^-34
        @test d.ω ≈ -54321 * 2.0^-32 * PI
        @test d.integrity_status_flag == true
        @test d.l2c_phasing == false
        @test !state.is_shifted_by_180_degrees
        # Positioning needs message types 11 and 30 as well.
        @test !GNSSDecoder.is_decoding_completed_for_positioning(d)
        @test !is_sat_healthy(state)
    end

    @testset "Phase-inverted stream decodes with is_shifted_by_180_degrees" begin
        enc = L5ITestEncoder()
        stream =
            fec_encode_soft(enc, [falses(20); build_mt10(); build_mt10(; tow_count = 1235)])
        state = decode(GPSL5IDecoderState(9), -stream, length(stream))
        @test state.is_shifted_by_180_degrees
        @test state.raw_data.last_message_id == 10
        @test state.raw_data.WN == 2345
    end

    @testset "reset_decoder_state clears in-flight state, keeps decoded data" begin
        enc = L5ITestEncoder()
        stream = fec_encode_soft(enc, [build_mt10(); build_mt10(; tow_count = 1235)])
        state = decode(GPSL5IDecoderState(9), stream, length(stream))
        @test state.raw_data.WN == 2345

        state = reset_decoder_state(state)
        @test isempty(state.cache.soft_buffer)
        @test isnothing(state.raw_data.TOW)
        @test isnothing(state.data.TOW)
        @test isnothing(state.num_bits_after_valid_syncro_sequence)
        @test !state.is_shifted_by_180_degrees
        # raw_data keeps WN, ephemeris fields, etc.
        @test state.raw_data.WN == 2345
    end

    # --- Spirent-derived recording fixture (committed; always runs) ---------

    """
    Assert the decoded `state` matches Spirent's reference decode (PRN 25).
    """
    function assert_spirent_golden(state)
        d = state.data
        # Header of the last decoded message (#25 of the recording, message
        # type 10, message TOW count 43225).
        @test d.last_message_id == 10
        @test d.TOW == 43225 * 6
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
        @test d.e ≈ 0.030000000027939677
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

        # Message type 32 — EOP. Spirent's reference decoder predates
        # IS-GPS-705J and prints the 31-bit field as ΔUT1 with the old 2⁻²⁴
        # scale (-0.25315952301025391 s); 705J defines it as ΔUT_GPS with
        # scale 2⁻²³ (Table 20-VII), i.e. exactly twice that.
        @test d.t_EOP == 259200
        @test d.PM_X ≈ 0.10612583160400391
        @test d.PM_X_dot ≈ 0.0013251304626464844
        @test d.PM_Y ≈ 0.4459381103515625
        @test d.PM_Y_dot ≈ 0.00040197372436523438
        @test d.ΔUT_GPS ≈ 2 * -0.25315952301025391
        @test d.ΔUT_GPS_dot ≈ -0.00046199560165405273  # rate scale unchanged in 705J

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
        # field is 82459 = 5153.6875 (matches the same constellation's L1C-D
        # midi almanacs).
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

        # Text messages (the simulator broadcasts "Spirent Communications
        # L2C test text message." across pages; the first cycle carries
        # page 1 of each).
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
        messages = load_l5i_fixture_messages()
        @test all(crc24q(m) == 0 for m in messages)
        # The recording starts at the first message: TOW count 43201 in
        # message 1 stamps the start of message 2.
        stream = fec_encode_soft(L5ITestEncoder(), reduce(vcat, messages))
        @test length(stream) == 26 * 600

        # Of the 26 encoded messages the last cannot complete (its sync
        # window needs the next message's first 16 symbols), so 25 decode.
        state = decode(GPSL5IDecoderState(25), stream, length(stream))
        @test !state.is_shifted_by_180_degrees
        assert_spirent_golden(state)

        # Polarity-inverted stream must decode identically with the
        # 180°-flip flag set.
        state_inv = decode(GPSL5IDecoderState(25), -stream, length(stream))
        @test state_inv.is_shifted_by_180_degrees
        @test state_inv.data == state.data

        # Mid-stream acquisition at an arbitrary (odd) symbol offset: the
        # decoder must slide to the next message boundary and decode the
        # remaining messages, losing only the partially received first one.
        offset = 351
        state_mid =
            decode(GPSL5IDecoderState(25), stream[(offset+1):end], length(stream) - offset)
        assert_spirent_golden(state_mid)

        # Noisy soft symbols: moderate Gaussian noise on the ±1 LLRs must
        # not cost any message at this SNR.
        noisy = stream .+ 0.4f0 .* randn(MersenneTwister(7), Float32, length(stream))
        state_noisy = decode(GPSL5IDecoderState(25), noisy, length(noisy))
        assert_spirent_golden(state_noisy)
    end

    @testset "decode_once stops at the first complete positioning set" begin
        messages = load_l5i_fixture_messages()
        stream = fec_encode_soft(L5ITestEncoder(), reduce(vcat, messages))
        state = decode(GPSL5IDecoderState(25), stream, length(stream); decode_once = true)
        # Message types 10, 11, 30 are the first three of the recording; the
        # positioning set is complete after message 3 and decoding stops
        # there (message 4 is type 15 — its text must not have been parsed).
        # `decode_once` freezes the validated `data` snapshot at the first
        # complete set (after message 3, type 30) while `raw_data` keeps
        # following the stream.
        @test GNSSDecoder.is_decoding_completed_for_positioning(state.data)
        @test state.data.last_message_id == 30
        @test isnothing(state.data.text_mt15)
        @test state.raw_data.text_mt15 == "Spirent Communications L2C te"
    end
end
