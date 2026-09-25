using GNSSDecoder: crc24q

# ---------------------------------------------------------------------------
# Reference CNAV transmit chain (test-only), shared by the GPS L5I and GPS L2C
# decoder tests.
#
# The CNAV message is identical on GPS L5I (IS-GPS-705J §20.3) and GPS L2C
# (IS-GPS-200N §30), so a single reference encoder feeds both. The rate-1/2,
# K=7 convolutional encoder (G1 = 0o171, G2 = 0o133, IS-GPS-705J §3.3.3.1.1 ≡
# IS-GPS-200N §3.3.3.1.1, Figure 3-14) is intentionally re-implemented here so
# the tests catch regressions in either the decoder's Viterbi or the message
# parsing. The encoder runs *continuously* across message boundaries — no tail
# bits, no reset — exactly like the satellite. Symbols are emitted as ±1
# Float32 (bit 0 ⇒ +1, bit 1 ⇒ -1), the package-wide LLR convention.
# ---------------------------------------------------------------------------

"""
Continuous K=7 rate-1/2 FEC encoder state (6-bit register, [s1..s6], s1 most recent).
"""
mutable struct CNAVTestEncoder
    register::UInt8
end
CNAVTestEncoder() = CNAVTestEncoder(0x00)

"""
Encode one bit; returns the (G1, G2) symbol pair as `Bool`s.
"""
function encode_bit!(enc::CNAVTestEncoder, bit::Bool)
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
function fec_encode_soft(enc::CNAVTestEncoder, bits::AbstractVector{Bool})
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

Health bits default to OK (0). `l1_health`/`l2_health`/`l5_health` set the
three signal-health bits (52/53/54) — the L2/L5 split is what distinguishes
the GPS L2C and L5I health checks.
"""
function build_mt10(;
    prn = 9,
    tow_count = 1234,
    l1_health = false,
    l2_health = false,
    l5_health = false,
)
    bits = falses(300)
    setbits!(bits, 1, 8, 0b10001011)   # preamble
    setbits!(bits, 9, 6, prn)
    setbits!(bits, 15, 6, 10)          # message type 10
    setbits!(bits, 21, 17, tow_count)
    setbits!(bits, 38, 1, 0)           # alert
    setbits!(bits, 39, 13, 2345)       # WN
    setbits!(bits, 52, 1, l1_health)   # L1 health
    setbits!(bits, 53, 1, l2_health)   # L2 health
    setbits!(bits, 54, 1, l5_health)   # L5 health
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

"""
Build a CNAV message type 40 (Integrity Support Message, IS-GPS-705J
Fig 20-14a) with a valid CRC. The recordings carry none, so the allocation
tests append this one to reach the ISM parser.
"""
function build_mt40(; prn = 9, tow_count = 1234, mask = UInt64(0x5555_5555_5555_5555 >> 1))
    bits = falses(300)
    setbits!(bits, 1, 8, 0b10001011)   # preamble
    setbits!(bits, 9, 6, prn)
    setbits!(bits, 15, 6, 40)          # message type 40
    setbits!(bits, 21, 17, tow_count)
    setbits!(bits, 38, 1, 0)           # alert
    setbits!(bits, 39, 4, 1)           # GNSS ID
    setbits!(bits, 43, 13, 2345)       # WN_ISM
    setbits!(bits, 56, 6, 7)           # TOW_ISM
    setbits!(bits, 62, 4, 3)           # t_correl
    setbits!(bits, 66, 4, 4)           # b_nom
    setbits!(bits, 70, 4, 5)           # γ_nom
    setbits!(bits, 74, 4, 6)           # R_sat
    setbits!(bits, 78, 4, 7)           # P_const
    setbits!(bits, 82, 4, 8)           # MFD
    setbits!(bits, 86, 3, 2)           # service level
    setbits!(bits, 89, 63, mask)       # SV mask
    setbits!(bits, 277, 24, crc24q(collect(bits[1:276])))
    return collect(bits)
end

"""
    cnav_allocation_stream(messages) -> Vector{Float32}

Soft-symbol stream for the `decode!` allocation tests: the recording's
`messages` (every message type the simulator emits, 10-15 and 30-37), then a
synthetic message type 40 and one more message type 10 whose preamble lets the
message type 40 sync, all FEC-encoded continuously.
"""
function cnav_allocation_stream(messages)
    extra = [build_mt40(; prn = 25); build_mt10(; prn = 25)]
    fec_encode_soft(CNAVTestEncoder(), [reduce(vcat, messages); extra])
end

"""
Check that `decode_allocations` over `cnav_allocation_stream` allocated nothing
and that the stream really exercised every store.
"""
function test_cnav_decode_allocation_free(make_state, messages)
    stream = cnav_allocation_stream(messages)
    allocations = decode_allocations(make_state, stream)
    @test allocations.fresh == 0 skip = !CHECK_ALLOCATIONS
    @test allocations.warm == 0 skip = !CHECK_ALLOCATIONS
    @test allocations.reset == 0 skip = !CHECK_ALLOCATIONS
    state = allocations.state
    @test is_decoding_completed_for_positioning(state)
    d = state.data
    @test d.last_message_type == 40
    @test length(d.reduced_almanacs) == 11
    @test length(d.midi_almanacs) == 1
    @test length(d.clock_corrections) == 7
    @test length(d.ephemeris_corrections) == 3
    @test d.text_mt15 == "Spirent Communications L2C te"
    @test d.text_mt36 == "Spirent Communicat"
    @test d.ism.service_level == 2

    # The polarity-inverted stream (message complemented after the Viterbi)
    # is allocation-free as well.
    inverted = decode_allocations(make_state, -stream)
    @test inverted.fresh == 0 skip = !CHECK_ALLOCATIONS
    @test inverted.warm == 0 skip = !CHECK_ALLOCATIONS
    @test inverted.reset == 0 skip = !CHECK_ALLOCATIONS
    @test inverted.state.is_shifted_by_180_degrees
    @test inverted.state.data == d

    # `decode` keeps value semantics on top of it: the input state is untouched.
    fresh = make_state()
    decoded = decode(fresh, stream, length(stream))
    @test decoded.data == d
    @test fresh == make_state()
    @test fresh.raw_data == GPSCNAVData()
end
