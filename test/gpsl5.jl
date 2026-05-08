# GPS L5 CNAV decoder tests.
#
# Test data is generated locally by encoding chosen CNAV message bits with the
# rate-1/2 K=7 (G1=171, G2=133 octal) convolutional encoder defined in
# IS-GPS-705J §3.3.3.1.1, then feeding the resulting symbol stream through
# the decoder. The encoder is intentionally re-implemented here so that the
# test catches regressions in either the decoder's Viterbi or the message-type
# parsing.

# CRC-24Q matching IS-GPS-705J §20.3.5.1 (g_i = 1 for i ∈
# {0, 1, 3, 4, 5, 6, 7, 10, 11, 14, 17, 18, 23, 24}). Returns the 24-bit CRC
# of `bits` (vector of UInt8 0/1).
function _gpsl5_test_crc24q(bits::Vector{UInt8})
    poly = UInt32(0x1864cfb)
    reg = UInt32(0)
    for b in bits
        reg = (reg << 1) | UInt32(b & 1)
        if (reg & 0x01000000) != 0
            reg ⊻= poly
        end
    end
    for _ = 1:24
        reg <<= 1
        if (reg & 0x01000000) != 0
            reg ⊻= poly
        end
    end
    return reg & 0x00ffffff
end

mutable struct _GPSL5TestEncoder
    state::UInt8
end
_GPSL5TestEncoder() = _GPSL5TestEncoder(0x00)

function _gpsl5_test_encode_bit!(enc::_GPSL5TestEncoder, u::UInt8)
    s = enc.state
    s1 = (s >> 5) & 0x01
    s3 = (s >> 3) & 0x01
    s4 = (s >> 2) & 0x01
    s5 = (s >> 1) & 0x01
    s6 =  s       & 0x01
    y1 = u ⊻ s3 ⊻ s4 ⊻ s5 ⊻ s6   # G1 = 171
    y2 = u ⊻ s1 ⊻ s3 ⊻ s4 ⊻ s6   # G2 = 133
    enc.state = (((u & 0x01) << 5) | (s >> 1)) & 0x3f
    return (y1, y2)
end

function _gpsl5_test_fec_encode(bits::Vector{UInt8})
    enc = _GPSL5TestEncoder()
    syms = UInt8[]
    for b in bits
        y1, y2 = _gpsl5_test_encode_bit!(enc, UInt8(b))
        push!(syms, y1, y2)
    end
    return syms
end

function _gpsl5_test_set_bits!(bits::Vector{UInt8}, start::Int, length::Int, value::Integer)
    v = UInt64(value & ((UInt64(1) << length) - UInt64(1)))
    for i = 0:length-1
        bits[start + i] = UInt8((v >> (length - 1 - i)) & 0x01)
    end
end

# Wide setter for ≥33-bit fields.
function _gpsl5_test_set_bits_wide!(bits::Vector{UInt8}, start::Int, length::Int, value::Integer)
    v = (value < 0 ? typemax(UInt128) - UInt128(-value) + UInt128(1) : UInt128(value)) &
        ((UInt128(1) << length) - UInt128(1))
    for i = 0:length-1
        bits[start + i] = UInt8((v >> (length - 1 - i)) & 0x01)
    end
end

function _gpsl5_test_build_mt10()
    bits = zeros(UInt8, 300)
    _gpsl5_test_set_bits!(bits, 1, 8, 0x8b)        # preamble
    _gpsl5_test_set_bits!(bits, 9, 6, 9)           # PRN
    _gpsl5_test_set_bits!(bits, 15, 6, 10)         # MT 10
    _gpsl5_test_set_bits!(bits, 21, 17, 1234)      # TOW count → 7404 s
    _gpsl5_test_set_bits!(bits, 38, 1, 0)          # alert
    _gpsl5_test_set_bits!(bits, 39, 13, 2345)      # WN
    _gpsl5_test_set_bits!(bits, 52, 1, 0)          # L1 health
    _gpsl5_test_set_bits!(bits, 53, 1, 0)          # L2 health
    _gpsl5_test_set_bits!(bits, 54, 1, 0)          # L5 health
    _gpsl5_test_set_bits!(bits, 55, 11, 100)       # t_op (×300 = 30000)
    _gpsl5_test_set_bits!(bits, 66, 5, 0)          # URA_ED
    _gpsl5_test_set_bits!(bits, 71, 11, 100)       # t_oe (×300 = 30000)
    _gpsl5_test_set_bits!(bits, 82, 26, 0)         # ΔA
    _gpsl5_test_set_bits!(bits, 108, 25, 0)        # A_dot
    _gpsl5_test_set_bits!(bits, 133, 17, 100)      # Δn0
    _gpsl5_test_set_bits!(bits, 150, 23, 0)        # Δn0_dot
    _gpsl5_test_set_bits_wide!(bits, 173, 33, 12345)    # M_0
    _gpsl5_test_set_bits_wide!(bits, 206, 33, 1000)     # e
    _gpsl5_test_set_bits_wide!(bits, 239, 33, -54321)   # ω
    _gpsl5_test_set_bits!(bits, 272, 1, 1)         # integrity status flag
    _gpsl5_test_set_bits!(bits, 273, 1, 0)         # L2C phasing
    crc = _gpsl5_test_crc24q(bits[1:276])
    _gpsl5_test_set_bits!(bits, 277, 24, crc)
    return bits
end

@testset "GPS L5 constructor" begin
    gpsl5 = GPSL5()
    state_a = GPSL5DecoderState(13)
    state_b = GNSSDecoderState(gpsl5, 13)
    @test state_a.prn == state_b.prn
    @test state_a.raw_buffer == state_b.raw_buffer
    @test state_a.buffer == state_b.buffer
    @test state_a.raw_data == state_b.raw_data
    @test state_a.data == state_b.data
    @test state_a.constants == state_b.constants
    @test state_a.cache == state_b.cache
    @test state_a.num_bits_buffered == state_b.num_bits_buffered
    @test state_a.is_shifted_by_180_degrees == state_b.is_shifted_by_180_degrees
end

@testset "GPS L5 CRC-24Q" begin
    # All-zero message with computed CRC must verify clean.
    bits = zeros(UInt8, 300)
    crc = _gpsl5_test_crc24q(bits[1:276])
    _gpsl5_test_set_bits!(bits, 277, 24, crc)
    @test _gpsl5_test_crc24q(bits) == 0

    # Random message
    msg = _gpsl5_test_build_mt10()
    @test _gpsl5_test_crc24q(msg) == 0
end

@testset "GPS L5 round-trip MT10" begin
    mt10 = _gpsl5_test_build_mt10()
    @test _gpsl5_test_crc24q(mt10) == 0

    # Build symbol stream: warmup (encoder ramp-up) + message + flush bits.
    warmup = zeros(UInt8, 100)
    trailing = zeros(UInt8, 50)
    input = vcat(warmup, mt10, trailing)
    symbols = _gpsl5_test_fec_encode(input)

    state = GPSL5DecoderState(9)
    for s in symbols
        state = decode(state, UInt8(s & 1), 1)
    end

    @test state.cache.locked_part == 1
    @test state.raw_data.last_message_id == 10
    @test state.raw_data.WN == 2345
    @test state.raw_data.TOW == 1234 * 6
    @test state.raw_data.alert_flag == false
    @test state.raw_data.signal_health_l1 == false
    @test state.raw_data.signal_health_l2 == false
    @test state.raw_data.signal_health_l5 == false
    @test state.raw_data.t_op == 30_000
    @test state.raw_data.URA_ED == 0
    @test state.raw_data.t_oe == 30_000
    @test state.raw_data.ΔA == 0.0
    @test state.raw_data.A_dot == 0.0
    @test state.raw_data.integrity_status_flag == true
    @test state.raw_data.l2c_phasing == false
    @test state.raw_data.M_0 ≈ 12345 * 3.1415926535898 / Float64(1 << 32)
    @test state.raw_data.e ≈ 1000 / Float64(1 << 34)
    @test state.raw_data.ω ≈ -54321 * 3.1415926535898 / Float64(1 << 32)
    @test is_sat_healthy(state) == false  # data not yet promoted to .data (validate needs MT11 + MT30)
end

@testset "GPS L5 reset_decoder_state" begin
    mt10 = _gpsl5_test_build_mt10()
    input = vcat(zeros(UInt8, 100), mt10, zeros(UInt8, 50))
    symbols = _gpsl5_test_fec_encode(input)
    state = GPSL5DecoderState(9)
    for s in symbols
        state = decode(state, UInt8(s & 1), 1)
    end
    state = reset_decoder_state(state)
    @test state.raw_buffer == 0
    @test state.buffer == 0
    @test isnothing(state.raw_data.TOW)
    @test isnothing(state.data.TOW)
    @test state.cache == GNSSDecoder.GPSL5Cache()
    @test state.num_bits_buffered == 0
    @test isnothing(state.num_bits_after_valid_syncro_sequence)
    # raw_data still has WN, ephemeris fields, etc. preserved
    @test state.raw_data.WN == 2345
end

@testset "GPS L5 phase-inverted preamble" begin
    # Send the 180°-inverted symbol stream — decoder should still acquire lock
    # by matching the inverted preamble (0x74).
    mt10 = _gpsl5_test_build_mt10()
    input = vcat(zeros(UInt8, 100), mt10, zeros(UInt8, 50))
    symbols = _gpsl5_test_fec_encode(input)
    inverted = UInt8.(symbols .⊻ 0x01)
    state = GPSL5DecoderState(9)
    for s in inverted
        state = decode(state, UInt8(s & 1), 1)
    end
    @test state.cache.locked_part != 0
    @test state.raw_data.last_message_id == 10
    @test state.raw_data.WN == 2345
    @test state.raw_data.TOW == 1234 * 6
end
