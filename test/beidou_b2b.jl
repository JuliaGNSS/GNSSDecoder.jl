using Test
using GNSSDecoder
using GNSSSignals
using Dictionaries
using Random
using GNSSSignals: Hz

# The GF(2^6) reference arithmetic + systematic encoder for the transmit
# chain (`beidou_ldpc_encode`, `gf64_symbols_to_bits`, `gf64_bits_to_symbols`,
# `BEIDOU_LDPC_CODES`). Guarded includes inside make re-inclusion cheap.
isdefined(Main, :beidou_ldpc_encode) ||
    include(joinpath(@__DIR__, "..", "scripts", "generate_beidou_alist.jl"))

# ---- B-CNAV3 transmit chain (BDS-SIS-ICD-B2b-1.0 §6.2) ----------------------

"""
Write `value`'s low `len` bits MSB-first into `bits[start:start+len-1]` (1-based).
"""
function b2b_set_field!(bits::AbstractVector{Bool}, start::Int, len::Int, value::Integer)
    u = UInt64(value & ((Int64(1) << len) - 1))
    for i = 1:len
        bits[start+i-1] = (u >> (len - i)) & 1 == 1
    end
    bits
end

"""
Assemble one 486-bit B-CNAV3 message: MesType(6) SOW(20) + caller-filled data
via `fill!(bits)`, with the CRC-24Q appended over the leading 462 bits.
"""
function b2b_build_message(fill!::Function, mestype::Int, sow_field::Int)
    bits = zeros(Bool, 486)
    b2b_set_field!(bits, 1, 6, mestype)
    b2b_set_field!(bits, 7, 20, sow_field)
    fill!(bits)
    crc = GNSSDecoder.crc24q(view(bits, 1:462))
    b2b_set_field!(bits, 463, 24, Int(crc))
    @assert GNSSDecoder.crc24q(bits) == 0
    bits
end

"""
Encode a 486-bit message into one 1000-symbol B-CNAV3 frame of ±1 soft
symbols (positive ⇒ bit 0): Pre(16, 0xEB90) + PRN(6) + Rev(6) + LDPC(972).
"""
function b2b_frame_symbols(message_bits::Vector{Bool}; prn::Int, invert::Bool = false)
    @assert length(message_bits) == 486
    codeword =
        beidou_ldpc_encode(BEIDOU_LDPC_CODES.bcnv3, gf64_bits_to_symbols(message_bits))
    frame_bits = zeros(Bool, 1000)
    b2b_set_field!(frame_bits, 1, 16, 0xEB90)
    b2b_set_field!(frame_bits, 17, 6, prn)
    # bits 23-28 reserved (zeros)
    frame_bits[29:1000] = gf64_symbols_to_bits(codeword)
    sign = invert ? -1.0f0 : 1.0f0
    Float32[b ? -sign : sign for b in frame_bits]
end

# Reference field values for the three message types. Chosen to exercise
# signed (negative) fields, unsigned fields, and every scale factor.
const B2B_SOW_FIELD_BASE = 100_000       # broadcast value == seconds (LSB 1 s)

function b2b_mt10_message(; sow_field = B2B_SOW_FIELD_BASE)
    b2b_build_message(10, sow_field) do bits
        b2b_set_field!(bits, 31, 11, 1200)      # t_0e = 1200*300
        b2b_set_field!(bits, 42, 2, 0b11)       # MEO
        b2b_set_field!(bits, 44, 26, -5000)     # ΔA
        b2b_set_field!(bits, 70, 25, -321)      # A_dot
        b2b_set_field!(bits, 95, 17, 12345)     # Δn_0
        b2b_set_field!(bits, 112, 23, -54321)   # Δn_0_dot
        b2b_set_field!(bits, 135, 33, -1234567) # M_0
        b2b_set_field!(bits, 168, 33, 87654321) # e
        b2b_set_field!(bits, 201, 33, 2345678)  # ω
        b2b_set_field!(bits, 234, 33, -7654321) # Ω_0
        b2b_set_field!(bits, 267, 33, 3456789)  # i_0
        b2b_set_field!(bits, 300, 19, -98765)   # Ω_dot
        b2b_set_field!(bits, 319, 15, 4321)     # i_dot
        b2b_set_field!(bits, 334, 16, -1111)    # C_is
        b2b_set_field!(bits, 350, 16, 2222)     # C_ic
        b2b_set_field!(bits, 366, 24, -33333)   # C_rs
        b2b_set_field!(bits, 390, 24, 44444)    # C_rc
        b2b_set_field!(bits, 414, 21, -5555)    # C_us
        b2b_set_field!(bits, 435, 21, 6666)     # C_uc
        b2b_set_field!(bits, 456, 1, 0)         # DIF
        b2b_set_field!(bits, 457, 1, 1)         # SIF
        b2b_set_field!(bits, 458, 1, 0)         # AIF
        b2b_set_field!(bits, 459, 4, 9)         # SISMAI
    end
end

function b2b_mt30_message(; sow_field = B2B_SOW_FIELD_BASE + 1, HS = 0)
    b2b_build_message(30, sow_field) do bits
        b2b_set_field!(bits, 27, 13, 913)       # WN
        b2b_set_field!(bits, 44, 11, 1500)      # t_0c = 1500*300
        b2b_set_field!(bits, 55, 25, -400000)   # a_f0
        b2b_set_field!(bits, 80, 22, 123456)    # a_f1
        b2b_set_field!(bits, 102, 11, -321)     # a_f2
        b2b_set_field!(bits, 113, 12, -1000)    # T_GD_B2bI
        b2b_set_field!(bits, 125, 10, 200)      # α_bdgim_1 (unsigned)
        b2b_set_field!(bits, 135, 8, -50)       # α_bdgim_2 (signed)
        b2b_set_field!(bits, 143, 8, 60)        # α_bdgim_3 (unsigned)
        b2b_set_field!(bits, 151, 8, 70)        # α_bdgim_4 (unsigned)
        b2b_set_field!(bits, 159, 8, 80)        # α_bdgim_5 (unsigned, scale -2⁻³)
        b2b_set_field!(bits, 167, 8, -90)       # α_bdgim_6
        b2b_set_field!(bits, 175, 8, 100)       # α_bdgim_7
        b2b_set_field!(bits, 183, 8, -110)      # α_bdgim_8
        b2b_set_field!(bits, 191, 8, 120)       # α_bdgim_9
        b2b_set_field!(bits, 199, 16, -20000)   # A_0UTC
        b2b_set_field!(bits, 215, 13, 3000)     # A_1UTC
        b2b_set_field!(bits, 228, 7, -40)       # A_2UTC
        b2b_set_field!(bits, 235, 8, 4)         # Δt_LS
        b2b_set_field!(bits, 243, 16, 37800)    # t_0t = 37800*16
        b2b_set_field!(bits, 259, 13, 913)      # WN_0t
        b2b_set_field!(bits, 272, 13, 829)      # WN_LSF
        b2b_set_field!(bits, 285, 3, 6)         # DN
        b2b_set_field!(bits, 288, 8, 5)         # Δt_LSF
        b2b_set_field!(bits, 296, 16, 21000)    # t_EOP = 21000*16
        b2b_set_field!(bits, 312, 21, -300000)  # PM_X
        b2b_set_field!(bits, 333, 15, 8000)     # PM_X_dot
        b2b_set_field!(bits, 348, 21, 400000)   # PM_Y
        b2b_set_field!(bits, 369, 15, -9000)    # PM_Y_dot
        b2b_set_field!(bits, 384, 31, -12345678)# ΔUT1
        b2b_set_field!(bits, 415, 19, 200000)   # ΔUT1_dot
        b2b_set_field!(bits, 434, 11, 1750)     # t_op (raw)
        b2b_set_field!(bits, 445, 5, 12)        # SISAI_ocb
        b2b_set_field!(bits, 450, 3, 3)         # SISAI_oc1
        b2b_set_field!(bits, 453, 3, 5)         # SISAI_oc2
        b2b_set_field!(bits, 456, 5, 7)         # SISAI_oe
        b2b_set_field!(bits, 461, 2, HS)        # HS
    end
end

function b2b_mt40_message(; sow_field = B2B_SOW_FIELD_BASE + 2)
    b2b_build_message(40, sow_field) do bits
        b2b_set_field!(bits, 27, 3, 1)          # GNSS_ID = GPS
        b2b_set_field!(bits, 30, 13, 913)       # WN_0BGTO
        b2b_set_field!(bits, 43, 16, 11250)     # t_0BGTO = 11250*16
        b2b_set_field!(bits, 59, 16, -1234)     # A_0BGTO
        b2b_set_field!(bits, 75, 13, 987)       # A_1BGTO
        b2b_set_field!(bits, 88, 7, -12)        # A_2BGTO
        # Midi almanac for PRN 30 (IGSO):
        b2b_set_field!(bits, 95, 6, 30)         # PRN_a
        b2b_set_field!(bits, 101, 2, 0b10)      # IGSO
        b2b_set_field!(bits, 103, 13, 913)      # WN_a
        b2b_set_field!(bits, 116, 8, 100)       # t_0a = 100*4096
        b2b_set_field!(bits, 124, 11, 500)      # e
        b2b_set_field!(bits, 135, 11, -200)     # δi
        b2b_set_field!(bits, 146, 17, 103900)   # sqrt_A = 103900*2⁻⁴
        b2b_set_field!(bits, 163, 16, -10000)   # Ω_0
        b2b_set_field!(bits, 179, 11, -700)     # Ω_dot
        b2b_set_field!(bits, 190, 16, 20000)    # ω
        b2b_set_field!(bits, 206, 16, -30000)   # M_0
        b2b_set_field!(bits, 222, 11, -600)     # a_f0
        b2b_set_field!(bits, 233, 10, 300)      # a_f1
        b2b_set_field!(bits, 243, 8, 0)         # health
        # Reduced-almanac reference time (Table 7-15):
        b2b_set_field!(bits, 251, 13, 913)      # WN_a (reduced)
        b2b_set_field!(bits, 264, 8, 101)       # t_0a = 101*4096
        # Reduced almanac 1: PRN 7, MEO
        b2b_set_field!(bits, 272, 6, 7)
        b2b_set_field!(bits, 278, 2, 0b11)
        b2b_set_field!(bits, 280, 8, -20)       # δA = -20*512
        b2b_set_field!(bits, 288, 7, 30)        # Ω_0
        b2b_set_field!(bits, 295, 7, -25)       # Φ_0
        b2b_set_field!(bits, 302, 8, 0)         # health
        # Reduced almanac 2: PRN 8, GEO
        b2b_set_field!(bits, 310, 6, 8)
        b2b_set_field!(bits, 316, 2, 0b01)
        b2b_set_field!(bits, 318, 8, 15)
        b2b_set_field!(bits, 326, 7, -30)
        b2b_set_field!(bits, 333, 7, 20)
        b2b_set_field!(bits, 340, 8, 2)
        # Blocks 3-5 left empty (PRN_a = 0 marks an empty block).
    end
end

# Feed a sequence of frames (plus the trailing next-frame preamble that the
# both-ends sync rule needs) through the public streaming API.
function b2b_decode_frames(state, frames::Vector{Vector{Float32}}; invert::Bool = false)
    tail_message = b2b_mt10_message(; sow_field = B2B_SOW_FIELD_BASE + 3)
    tail = b2b_frame_symbols(tail_message; prn = state.prn, invert)[1:16]
    stream = vcat(frames..., tail)
    decode(state, stream, length(stream))
end

@testset "BeiDou B2b (B-CNAV3)" begin
    prn = 26
    PI = GNSSDecoder.GNSS_PI

    @testset "Signal metadata" begin
        state = BeiDouB2bDecoderState(prn)
        @test get_signal_type(state) == BeiDouB2bI
        @test get_data_frequency(state) == 1000Hz
        @test get_time_system_name(state) == "BeiDou Time"
        @test get_constellation_name(state) == "BeiDou"
        # Pin the carrier, not the label: GNSSSignals v4 renamed the 1207.14 MHz
        # band type B2b -> E5b (Galileo's name for the same carrier), so
        # `get_band_name` answers "B2b" on 3.9 and "E5b" on 4. The frequency is
        # what the decoder's callers actually key off, and it is unchanged.
        @test get_center_frequency(get_band(state)) == 1_207_140_000Hz
    end

    @testset "MT10 + MT30 + MT40 decode and promotion" begin
        state = BeiDouB2bDecoderState(prn)
        frames = [
            b2b_frame_symbols(b2b_mt10_message(); prn),
            b2b_frame_symbols(b2b_mt30_message(); prn),
            b2b_frame_symbols(b2b_mt40_message(); prn),
        ]
        state = b2b_decode_frames(state, frames)
        @test !state.is_shifted_by_180_degrees

        # MT10 — ephemeris + integrity:
        d = state.raw_data
        @test d.t_0e == 1200 * 300
        @test d.sat_type == 3
        @test d.ΔA == -5000 * 2.0^-9
        @test d.A_dot == -321 * 2.0^-21
        @test d.Δn_0 == 12345 * 2.0^-44 * PI
        @test d.Δn_0_dot == -54321 * 2.0^-57 * PI
        @test d.M_0 == -1234567 * 2.0^-32 * PI
        @test d.e == 87654321 * 2.0^-34
        @test d.ω == 2345678 * 2.0^-32 * PI
        @test d.Ω_0 == -7654321 * 2.0^-32 * PI
        @test d.i_0 == 3456789 * 2.0^-32 * PI
        @test d.Ω_dot == -98765 * 2.0^-44 * PI
        @test d.i_dot == 4321 * 2.0^-44 * PI
        @test d.C_is == -1111 * 2.0^-30
        @test d.C_ic == 2222 * 2.0^-30
        @test d.C_rs == -33333 * 2.0^-8
        @test d.C_rc == 44444 * 2.0^-8
        @test d.C_us == -5555 * 2.0^-30
        @test d.C_uc == 6666 * 2.0^-30
        @test d.DIF === false
        @test d.SIF === true
        @test d.AIF === false
        @test d.SISMAI == 9

        # MT30 — clock/TGD/iono/UTC/EOP/SISAI/health:
        @test d.WN == 913
        @test d.t_0c == 1500 * 300
        @test d.a_f0 == -400000 * 2.0^-34
        @test d.a_f1 == 123456 * 2.0^-50
        @test d.a_f2 == -321 * 2.0^-66
        @test d.T_GD_B2bI == -1000 * 2.0^-34
        @test d.α_bdgim_1 == 200 * 2.0^-3
        @test d.α_bdgim_2 == -50 * 2.0^-3
        @test d.α_bdgim_3 == 60 * 2.0^-3
        @test d.α_bdgim_4 == 70 * 2.0^-3
        @test d.α_bdgim_5 == 80 * -(2.0^-3)
        @test d.α_bdgim_6 == -90 * 2.0^-3
        @test d.α_bdgim_7 == 100 * 2.0^-3
        @test d.α_bdgim_8 == -110 * 2.0^-3
        @test d.α_bdgim_9 == 120 * 2.0^-3
        @test d.A_0UTC == -20000 * 2.0^-35
        @test d.A_1UTC == 3000 * 2.0^-51
        @test d.A_2UTC == -40 * 2.0^-68
        @test d.Δt_LS == 4
        @test d.t_0t == 37800 * 16
        @test d.WN_0t == 913
        @test d.WN_LSF == 829
        @test d.DN == 6
        @test d.Δt_LSF == 5
        @test d.t_EOP == 21000 * 16
        @test d.PM_X == -300000 * 2.0^-20
        @test d.PM_X_dot == 8000 * 2.0^-21
        @test d.PM_Y == 400000 * 2.0^-20
        @test d.PM_Y_dot == -9000 * 2.0^-21
        @test d.ΔUT1 == -12345678 * 2.0^-24
        @test d.ΔUT1_dot == 200000 * 2.0^-25
        @test d.t_op == 1750
        @test d.SISAI_ocb == 12
        @test d.SISAI_oc1 == 3
        @test d.SISAI_oc2 == 5
        @test d.SISAI_oe == 7
        @test d.HS == 0

        # MT40 — BGTO + almanacs:
        @test d.GNSS_ID == 1
        @test d.WN_0BGTO == 913
        @test d.t_0BGTO == 11250 * 16
        @test d.A_0BGTO == -1234 * 2.0^-35
        @test d.A_1BGTO == 987 * 2.0^-51
        @test d.A_2BGTO == -12 * 2.0^-68
        @test d.WN_a == 913
        @test d.t_0a == 101 * 4096
        midi = d.midi_almanacs[30]
        @test midi.PRN_a == 30
        @test midi.sat_type == 2
        @test midi.WN_a == 913
        @test midi.t_0a == 100 * 4096
        @test midi.e == 500 * 2.0^-16
        @test midi.δi == -200 * 2.0^-14 * PI
        @test midi.sqrt_A == 103900 * 2.0^-4
        @test midi.Ω_0 == -10000 * 2.0^-15 * PI
        @test midi.Ω_dot == -700 * 2.0^-33 * PI
        @test midi.ω == 20000 * 2.0^-15 * PI
        @test midi.M_0 == -30000 * 2.0^-15 * PI
        @test midi.a_f0 == -600 * 2.0^-20
        @test midi.a_f1 == 300 * 2.0^-37
        @test midi.health == 0
        red7 = d.reduced_almanacs[7]
        @test red7.sat_type == 3
        @test red7.δA == -20 * 512.0
        @test red7.Ω_0 == 30 * 2.0^-6 * PI
        @test red7.Φ_0 == -25 * 2.0^-6 * PI
        @test red7.health == 0
        red8 = d.reduced_almanacs[8]
        @test red8.sat_type == 1
        @test red8.δA == 15 * 512.0
        @test red8.health == 2
        @test length(d.reduced_almanacs) == 2

        # Promotion: SOW + WN + ephemeris (MT10) + clock (MT30) present.
        @test is_decoding_completed_for_positioning(state)
        # A reserved orbit type (0) leaves the semi-major axis unknowable: it
        # is what selects the `A_ref` the broadcast `ΔA` corrects, so the
        # satellite must not be used (ICD Table 7-6).
        @test !GNSSDecoder.is_ephemeris_decoded(BeiDouB2bData(d; sat_type = 0))
        @test all(
            GNSSDecoder.is_ephemeris_decoded(BeiDouB2bData(d; sat_type = t)) for t = 1:3
        )
        # The SOW stamps this frame's preamble, so the armed counter spans the
        # frame plus the next preamble: 1000 + 16 = 1016 symbols = 1.016 s.
        @test state.num_bits_after_valid_syncro_sequence == 1016
        @test state.data.SOW == B2B_SOW_FIELD_BASE + 2
        @test state.data.last_message_type == 40
        @test is_sat_healthy(state)
    end

    # The 1016 above is armed by `validate_data`, and by nothing else: a frame
    # that clears CRC but does not complete the positioning set must leave the
    # counter alone, because `get_time_of_week` publishes the *validated* SOW and
    # the two only describe a time as a pair. On B2b the completeness gate holds
    # once MT10 and MT30 are in, so this is defensive — it is the same split that
    # fixes a live once-per-message-cycle 3-second error on B2a, whose MT10/MT11
    # adjacency gate really does skip promotions (see src/beidou/b2b.jl).
    @testset "Only promotion arms the symbol counter" begin
        state = BeiDouB2bDecoderState(prn)
        # MT10 alone: ephemeris, but neither the clock nor the week number, so
        # `is_decoding_completed_for_positioning` is still false.
        state = b2b_decode_frames(
            state,
            [b2b_frame_symbols(b2b_mt10_message(; sow_field = 5000); prn)],
        )
        @test state.raw_data.SOW == 5000                              # frame decoded ...
        @test !is_decoding_completed_for_positioning(state)
        @test isnothing(state.num_bits_after_valid_syncro_sequence)   # ... armed nothing
        # A second decoded frame does not arm it either.
        state = b2b_decode_frames(
            state,
            [b2b_frame_symbols(b2b_mt10_message(; sow_field = 5001); prn)],
        )
        @test state.raw_data.SOW == 5001
        @test isnothing(state.num_bits_after_valid_syncro_sequence)
        # The MT30 that completes the set promotes, and arms the counter to the
        # frame plus next-frame preamble whose SOW it just published.
        state = b2b_decode_frames(
            state,
            [b2b_frame_symbols(b2b_mt30_message(; sow_field = 5002); prn)],
        )
        @test state.data.SOW == 5002
        @test state.num_bits_after_valid_syncro_sequence == 1016
    end

    @testset "SOW stamps each frame; unhealthy satellite is flagged" begin
        state = BeiDouB2bDecoderState(prn)
        frames = [
            b2b_frame_symbols(b2b_mt10_message(; sow_field = 2000); prn),
            b2b_frame_symbols(b2b_mt30_message(; sow_field = 2001, HS = 1); prn),
        ]
        state = b2b_decode_frames(state, frames)
        @test state.raw_data.SOW == 2001
        @test state.data.HS == 1
        @test is_decoding_completed_for_positioning(state)
        @test !is_sat_healthy(state)
    end

    # The 20-bit B-CNAV3 SOW is a count of *seconds* (LSB 1 s), not the 3-second
    # count that B-CNAV2 broadcasts on B2a in 18 bits. The Chinese-language
    # original (BDS-SIS-ICD-B2b-1.0 中文版, 表 7-2) prints scale factor 1; only
    # the English edition misprints it as 3 (see the decoder comment in
    # src/beidou/b2b.jl). Pinning the scale needs an assertion the encoder
    # cannot cancel out: a B2b frame is 1000 symbols = 1 s, so decoding N
    # consecutive frames must advance SOW by exactly N-1 seconds, and the top
    # of the 20-bit range must still be a legal second-of-week.
    @testset "SOW resolution is 1 s per frame (B-CNAV3 is not B2a's 3 s count)" begin
        state = BeiDouB2bDecoderState(prn)
        sows = Int[]
        for k = 0:3
            state = b2b_decode_frames(
                state,
                [b2b_frame_symbols(b2b_mt10_message(; sow_field = 2000 + k); prn)],
            )
            push!(sows, state.raw_data.SOW)
        end
        # One frame elapses one second, so the stamps are consecutive seconds.
        @test diff(sows) == [1, 1, 1]
        @test sows[1] == 2000

        # A field value that is only a valid second-of-week under LSB 1 s: with
        # the 3 s scale this would decode to 1_814_397 s, far beyond a week.
        state = BeiDouB2bDecoderState(prn)
        state = b2b_decode_frames(
            state,
            [b2b_frame_symbols(b2b_mt10_message(; sow_field = 604_799); prn)],
        )
        @test state.raw_data.SOW == 604_799
        @test 0 <= state.raw_data.SOW < 604_800
    end

    @testset "180-degree phase-shifted stream decodes identically" begin
        state = BeiDouB2bDecoderState(prn)
        frames = [
            b2b_frame_symbols(b2b_mt10_message(); prn, invert = true),
            b2b_frame_symbols(b2b_mt30_message(); prn, invert = true),
        ]
        state = b2b_decode_frames(state, frames; invert = true)
        @test state.is_shifted_by_180_degrees
        @test state.raw_data.t_0e == 1200 * 300
        @test state.raw_data.a_f0 == -400000 * 2.0^-34
        @test is_decoding_completed_for_positioning(state)
    end

    @testset "Foreign PRN in the unencoded PRN field is rejected" begin
        state = BeiDouB2bDecoderState(prn)
        frames = [b2b_frame_symbols(b2b_mt10_message(); prn = prn + 1)]
        state = b2b_decode_frames(state, frames)
        @test state.raw_data.last_message_type == 0
        @test isnothing(state.raw_data.SOW)
    end

    @testset "CRC failure drops the frame" begin
        # A message whose CRC field is deliberately wrong is still a valid
        # LDPC codeword after encoding, so the BP decode converges cleanly and
        # only the CRC-24Q gate can (and must) reject it.
        bad = b2b_mt10_message()
        bad[470] = !bad[470]  # corrupt the CRC field itself
        state = BeiDouB2bDecoderState(prn)
        state = b2b_decode_frames(state, [b2b_frame_symbols(bad; prn)])
        @test state.raw_data.last_message_type == 0
        @test isnothing(state.raw_data.SOW)
    end

    @testset "Reserved message type contributes header only" begin
        msg = b2b_build_message(63, 4000) do bits
        end
        state = BeiDouB2bDecoderState(prn)
        state = b2b_decode_frames(state, [b2b_frame_symbols(msg; prn)])
        @test state.raw_data.last_message_type == 63
        @test state.raw_data.SOW == 4000
        @test isnothing(state.raw_data.t_0e)
        @test !is_decoding_completed_for_positioning(state)
    end

    @testset "Noisy symbols decode through the LDPC" begin
        rng = Random.MersenneTwister(0xB2B)
        state = BeiDouB2bDecoderState(prn)
        frames = [
            b2b_frame_symbols(b2b_mt10_message(); prn),
            b2b_frame_symbols(b2b_mt30_message(); prn),
        ]
        noisy = [f .* 2.0f0 .+ 0.8f0 .* randn(rng, Float32, length(f)) for f in frames]
        # Keep the tail preamble noiseless so the hard-decision preamble gate
        # is deterministic; the LDPC/CRC path is what is under test.
        state = b2b_decode_frames(state, noisy)
        @test state.raw_data.t_0e == 1200 * 300
        @test state.raw_data.a_f0 == -400000 * 2.0^-34
        @test is_decoding_completed_for_positioning(state)
    end

    @testset "reset_decoder_state clears SOW and validated data" begin
        state = BeiDouB2bDecoderState(prn)
        frames = [
            b2b_frame_symbols(b2b_mt10_message(); prn),
            b2b_frame_symbols(b2b_mt30_message(); prn),
        ]
        state = b2b_decode_frames(state, frames)
        @test is_decoding_completed_for_positioning(state)
        state = reset_decoder_state(state)
        @test isnothing(state.raw_data.SOW)
        @test state.raw_data.t_0e == 1200 * 300  # ephemeris survives for warm restart
        @test isnothing(state.data.SOW)
        @test !is_decoding_completed_for_positioning(state)
        @test isnothing(state.num_bits_after_valid_syncro_sequence)
        @test GNSSDecoder.num_bits_buffered(state) == 0
    end

    @testset "GNSSDecoderState(::BeiDouB2bI, prn) dispatch" begin
        state = GNSSDecoderState(BeiDouB2bI(), prn)
        @test state isa GNSSDecoderState{BeiDouB2bData}
        @test state.prn == prn
    end
end
