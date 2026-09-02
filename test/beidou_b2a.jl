using Test
using GNSSDecoder
using GNSSSignals
using Dictionaries
using GNSSDecoder: crc24q, GNSS_PI

# ---------------------------------------------------------------------------
# Reference B-CNAV2 transmit chain (test-only).
#
# Frames are built exactly as BDS-SIS-ICD-B2a-1.0 §6.2 prescribes: a 288-bit
# message (PRN, MesType, SOW, data, CRC-24Q over the leading 264 bits) is
# encoded with the 64-ary LDPC(96,48) reference encoder from
# `scripts/generate_beidou_alist.jl` (an independent GF(2^6) implementation of
# the ICD Annex algorithm — the decoder under test only ever sees the binary
# image through Aff3ct), mapped to bits MSB first, and prefixed with the
# 24-symbol preamble 0xE24DE8. Symbols are ±1 Float32 (bit 0 ⇒ +1, bit 1 ⇒
# −1), the package-wide LLR convention.
# ---------------------------------------------------------------------------

# GF(2^6) reference arithmetic + systematic LDPC encoder. Guarded: runtests.jl
# may already have included it via test/beidou_ldpc.jl.
isdefined(@__MODULE__, :beidou_ldpc_encode) ||
    isdefined(Main, :beidou_ldpc_encode) ||
    include(joinpath(@__DIR__, "..", "scripts", "generate_beidou_alist.jl"))

"""
Write `len` bits of `val` MSB-first at 1-based position `start` of `bits`.
"""
function b2a_setbits!(bits::BitVector, start::Int, len::Int, val::Integer)
    v = UInt64(unsigned(Int64(val)) & ((UInt64(1) << len) - UInt64(1)))
    for i = 1:len
        bits[start+i-1] = (v >> (len - i)) & 0x1 == 0x1
    end
    return bits
end

"""
Assemble one 288-bit B-CNAV2 message: header + `set_data!`-filled body + CRC.

`set_data!(bits)` fills the message-specific fields (bits 31-264). The CRC-24Q
over bits 1-264 is appended at bits 265-288 (ICD §6.2.1); `corrupt_crc` flips
the CRC's last bit to exercise the decoder's CRC gate.
"""
function build_b2a_message(
    set_data!;
    prn::Int,
    mestype::Int,
    sow_seconds::Int,
    corrupt_crc::Bool = false,
)
    @assert sow_seconds % 3 == 0
    bits = falses(288)
    b2a_setbits!(bits, 1, 6, prn)
    b2a_setbits!(bits, 7, 6, mestype)
    b2a_setbits!(bits, 13, 18, sow_seconds ÷ 3)
    set_data!(bits)
    crc = crc24q(Vector{Bool}(bits[1:264]))
    b2a_setbits!(bits, 265, 24, Int(crc))
    corrupt_crc && (bits[288] = !bits[288])
    @assert corrupt_crc || crc24q(Vector{Bool}(bits)) == 0
    return bits
end

"""
LDPC-encode a 288-bit message and prepend the preamble → 600 ±1 Float32 symbols.
"""
function b2a_frame_symbols(message::BitVector)
    syms = gf64_bits_to_symbols(Vector{Bool}(message))
    codeword = beidou_ldpc_encode(BEIDOU_LDPC_CODES.bcnv2, syms)
    encoded_bits = gf64_symbols_to_bits(codeword)
    preamble_bits = [(0xE24DE8 >> (24 - i)) & 0x1 == 0x1 for i = 1:24]
    Float32[b ? -1.0f0 : 1.0f0 for b in vcat(preamble_bits, encoded_bits)]
end

"""
The 24 preamble symbols alone (to trigger sync on the final buffered frame).
"""
b2a_trailing_preamble() =
    Float32[(0xE24DE8>>(24-i))&0x1 == 0x1 ? -1.0f0 : 1.0f0 for i = 1:24]

# ---------------------------------------------------------------------------
# Message builders — raw field integers chosen distinct so any bit-position or
# scale mix-up between two fields fails the assertions below.
# ---------------------------------------------------------------------------

const B2A_PRN = 19
const B2A_WN = 890
const B2A_IODE = 21
const B2A_IODC = (2 << 8) | B2A_IODE  # 533; IODE == IODC & 0xff (ICD §7.4.3)

# Integrity-flag block used in every message so the last decoded frame leaves
# a known state: DIF_B2a=0 SIF_B2a=0 AIF_B2a=1, SISMAI=5, DIF/SIF/AIF_B1C=1/0/0.
function set_flags!(bits, start)
    b2a_setbits!(bits, start, 1, 0)      # DIF_B2a
    b2a_setbits!(bits, start + 1, 1, 0)  # SIF_B2a
    b2a_setbits!(bits, start + 2, 1, 1)  # AIF_B2a
    b2a_setbits!(bits, start + 3, 4, 5)  # SISMAI
    b2a_setbits!(bits, start + 7, 1, 1)  # DIF_B1C
    b2a_setbits!(bits, start + 8, 1, 0)  # SIF_B1C
    b2a_setbits!(bits, start + 9, 1, 0)  # AIF_B1C
end

function set_clock!(bits, start)
    b2a_setbits!(bits, start, 11, 1005)        # t_0c raw (×300 = 301500 s)
    b2a_setbits!(bits, start + 11, 25, -5000000) # a_f0
    b2a_setbits!(bits, start + 36, 22, -150000)  # a_f1
    b2a_setbits!(bits, start + 58, 11, 300)      # a_f2
end

build_b2a_mt10(; sow, iode = B2A_IODE, prn = B2A_PRN, corrupt_crc = false) =
    build_b2a_message(; prn, mestype = 10, sow_seconds = sow, corrupt_crc) do bits
        b2a_setbits!(bits, 31, 13, B2A_WN)
        set_flags!(bits, 44)
        b2a_setbits!(bits, 54, 8, iode)
        b2a_setbits!(bits, 62, 11, 1000)        # t_0e raw (×300 = 300000 s)
        b2a_setbits!(bits, 73, 2, 3)            # SatType = MEO
        b2a_setbits!(bits, 75, 26, -12345)      # ΔA
        b2a_setbits!(bits, 101, 25, 345)        # A_dot
        b2a_setbits!(bits, 126, 17, -2000)      # Δn_0
        b2a_setbits!(bits, 143, 23, 12000)      # Δn_0_dot
        b2a_setbits!(bits, 166, 33, -123456789) # M_0
        b2a_setbits!(bits, 199, 33, 123456789)  # e (unsigned)
        b2a_setbits!(bits, 232, 33, 987654321)  # ω
    end

build_b2a_mt11(; sow, prn = B2A_PRN) =
    build_b2a_message(; prn, mestype = 11, sow_seconds = sow) do bits
        b2a_setbits!(bits, 31, 2, 0)            # HS healthy
        set_flags!(bits, 33)
        b2a_setbits!(bits, 43, 33, -1987654321) # Ω_0
        b2a_setbits!(bits, 76, 33, 555666777)   # i_0
        b2a_setbits!(bits, 109, 19, -90000)     # Ω_dot
        b2a_setbits!(bits, 128, 15, 9000)       # i_dot
        b2a_setbits!(bits, 143, 16, -2100)      # C_is
        b2a_setbits!(bits, 159, 16, 1500)       # C_ic
        b2a_setbits!(bits, 175, 24, -400000)    # C_rs
        b2a_setbits!(bits, 199, 24, 300000)     # C_rc
        b2a_setbits!(bits, 223, 21, -800000)    # C_us
        b2a_setbits!(bits, 244, 21, 700000)     # C_uc
    end

build_b2a_mt30(; sow, iodc = B2A_IODC) =
    build_b2a_message(; prn = B2A_PRN, mestype = 30, sow_seconds = sow) do bits
        b2a_setbits!(bits, 31, 2, 0)
        set_flags!(bits, 33)
        set_clock!(bits, 43)
        b2a_setbits!(bits, 112, 10, iodc)
        b2a_setbits!(bits, 122, 12, -800)   # T_GD_B2ap
        b2a_setbits!(bits, 134, 12, 600)    # ISC_B2ad
        b2a_setbits!(bits, 146, 10, 800)    # α₁ (unsigned)
        b2a_setbits!(bits, 156, 8, -100)    # α₂
        b2a_setbits!(bits, 164, 8, 200)     # α₃ (unsigned)
        b2a_setbits!(bits, 172, 8, 90)      # α₄ (unsigned)
        b2a_setbits!(bits, 180, 8, 40)      # α₅ (unsigned, scale −2⁻³)
        b2a_setbits!(bits, 188, 8, -3)      # α₆
        b2a_setbits!(bits, 196, 8, 17)      # α₇
        b2a_setbits!(bits, 204, 8, -128)    # α₈
        b2a_setbits!(bits, 212, 8, 127)     # α₉
        b2a_setbits!(bits, 220, 12, -1024)  # T_GD_B1Cp
    end

# Three reduced-almanac slots: PRN 7, PRN 8, and an empty (PRN 0) slot.
build_b2a_mt31(; sow) =
    build_b2a_message(; prn = B2A_PRN, mestype = 31, sow_seconds = sow) do bits
        b2a_setbits!(bits, 31, 2, 0)
        set_flags!(bits, 33)
        set_clock!(bits, 43)
        b2a_setbits!(bits, 112, 10, B2A_IODC)
        b2a_setbits!(bits, 122, 13, B2A_WN) # WN_a
        b2a_setbits!(bits, 135, 8, 100)     # t_0a raw (×2¹² = 409600 s)
        b2a_setbits!(bits, 143, 6, 7)       # PRN_a
        b2a_setbits!(bits, 149, 2, 3)       # MEO
        b2a_setbits!(bits, 151, 8, -50)     # δA
        b2a_setbits!(bits, 159, 7, 31)      # Ω_0
        b2a_setbits!(bits, 166, 7, -64)     # Φ_0
        b2a_setbits!(bits, 173, 8, 0)       # health
        b2a_setbits!(bits, 181, 6, 8)       # PRN_a
        b2a_setbits!(bits, 187, 2, 2)       # IGSO
        b2a_setbits!(bits, 189, 8, 100)     # δA
        b2a_setbits!(bits, 197, 7, -31)     # Ω_0
        b2a_setbits!(bits, 204, 7, 45)      # Φ_0
        b2a_setbits!(bits, 211, 8, 128)     # health (clock unhealthy)
        # third slot left all-zero: PRN_a = 0 ⇒ empty, must be skipped
    end

build_b2a_mt32(; sow) =
    build_b2a_message(; prn = B2A_PRN, mestype = 32, sow_seconds = sow) do bits
        b2a_setbits!(bits, 31, 2, 0)
        set_flags!(bits, 33)
        set_clock!(bits, 43)
        b2a_setbits!(bits, 112, 10, B2A_IODC)
        b2a_setbits!(bits, 122, 16, 30000)      # t_EOP raw (×2⁴ = 480000 s)
        b2a_setbits!(bits, 138, 21, -700000)    # PM_X
        b2a_setbits!(bits, 159, 15, 8000)       # PM_X_dot
        b2a_setbits!(bits, 174, 21, 650000)     # PM_Y
        b2a_setbits!(bits, 195, 15, -9000)      # PM_Y_dot
        b2a_setbits!(bits, 210, 31, -500000000) # ΔUT1
        b2a_setbits!(bits, 241, 19, 140000)     # ΔUT1_dot
    end

build_b2a_mt33(; sow) =
    build_b2a_message(; prn = B2A_PRN, mestype = 33, sow_seconds = sow) do bits
        b2a_setbits!(bits, 31, 2, 0)
        set_flags!(bits, 33)
        set_clock!(bits, 43)
        b2a_setbits!(bits, 112, 3, 1)       # GNSS_ID = GPS
        b2a_setbits!(bits, 115, 13, 880)    # WN_0BGTO
        b2a_setbits!(bits, 128, 16, 20000)  # t_0BGTO raw (×2⁴ = 320000 s)
        b2a_setbits!(bits, 144, 16, -20000) # A_0BGTO
        b2a_setbits!(bits, 160, 13, 3000)   # A_1BGTO
        b2a_setbits!(bits, 173, 7, -60)     # A_2BGTO
        b2a_setbits!(bits, 180, 6, 9)       # reduced almanac PRN_a
        b2a_setbits!(bits, 186, 2, 1)       # GEO
        b2a_setbits!(bits, 188, 8, 100)     # δA
        b2a_setbits!(bits, 196, 7, -31)     # Ω_0
        b2a_setbits!(bits, 203, 7, 45)      # Φ_0
        b2a_setbits!(bits, 210, 8, 2)       # health
        b2a_setbits!(bits, 218, 10, B2A_IODC)
        b2a_setbits!(bits, 228, 13, B2A_WN) # WN_a
        b2a_setbits!(bits, 241, 8, 100)     # t_0a raw
    end

build_b2a_mt34(; sow) =
    build_b2a_message(; prn = B2A_PRN, mestype = 34, sow_seconds = sow) do bits
        b2a_setbits!(bits, 31, 2, 0)
        set_flags!(bits, 33)
        b2a_setbits!(bits, 43, 11, 1200)  # t_op (broadcast count; ×300 s on decode)
        b2a_setbits!(bits, 54, 5, 17)     # SISAI_ocb
        b2a_setbits!(bits, 59, 3, 5)      # SISAI_oc1
        b2a_setbits!(bits, 62, 3, 2)      # SISAI_oc2
        set_clock!(bits, 65)
        b2a_setbits!(bits, 134, 10, B2A_IODC)
        b2a_setbits!(bits, 144, 16, -30000) # A_0UTC
        b2a_setbits!(bits, 160, 13, 2500)   # A_1UTC
        b2a_setbits!(bits, 173, 7, -50)     # A_2UTC
        b2a_setbits!(bits, 180, 8, 4)       # Δt_LS
        b2a_setbits!(bits, 188, 16, 35000)  # t_0t raw (×2⁴ = 560000 s)
        b2a_setbits!(bits, 204, 13, B2A_WN) # WN_0t
        b2a_setbits!(bits, 217, 13, 902)    # WN_LSF
        b2a_setbits!(bits, 230, 3, 6)       # DN
        b2a_setbits!(bits, 233, 8, 5)       # Δt_LSF
    end

build_b2a_mt40(; sow) =
    build_b2a_message(; prn = B2A_PRN, mestype = 40, sow_seconds = sow) do bits
        b2a_setbits!(bits, 31, 2, 0)
        set_flags!(bits, 33)
        b2a_setbits!(bits, 43, 5, 11)     # SISAI_oe
        b2a_setbits!(bits, 48, 11, 1200)  # t_op
        b2a_setbits!(bits, 59, 5, 17)     # SISAI_ocb
        b2a_setbits!(bits, 64, 3, 5)      # SISAI_oc1
        b2a_setbits!(bits, 67, 3, 2)      # SISAI_oc2
        b2a_setbits!(bits, 70, 6, 23)     # midi PRN_a
        b2a_setbits!(bits, 76, 2, 2)      # IGSO
        b2a_setbits!(bits, 78, 13, B2A_WN)
        b2a_setbits!(bits, 91, 8, 100)    # t_0a raw
        b2a_setbits!(bits, 99, 11, 1500)  # e (unsigned)
        b2a_setbits!(bits, 110, 11, -800) # δi
        b2a_setbits!(bits, 121, 17, 105000) # √A (unsigned)
        b2a_setbits!(bits, 138, 16, -20000) # Ω_0
        b2a_setbits!(bits, 154, 11, -700)   # Ω_dot
        b2a_setbits!(bits, 165, 16, 15000)  # ω
        b2a_setbits!(bits, 181, 16, -15000) # M_0
        b2a_setbits!(bits, 197, 11, -900)   # a_f0
        b2a_setbits!(bits, 208, 10, 333)    # a_f1
        b2a_setbits!(bits, 218, 8, 0)       # health
    end

"""
Feed a sequence of 288-bit messages (plus trailing preamble) into `state`.
"""
function decode_b2a_frames(state, messages; invert = false)
    stream = reduce(vcat, [b2a_frame_symbols(m) for m in messages])
    stream = vcat(stream, b2a_trailing_preamble())
    invert && (stream = -stream)
    decode(state, stream, length(stream))
end

@testset "BeiDou B2a B-CNAV2" begin
    PI = GNSS_PI
    t0 = 302400  # arbitrary mid-week SOW, multiple of 3

    @testset "Full broadcast cycle decodes every message type" begin
        state = BeiDouB2aDecoderState(B2A_PRN)
        messages = [
            build_b2a_mt10(; sow = t0),
            build_b2a_mt11(; sow = t0 + 3),
            build_b2a_mt30(; sow = t0 + 6),
            build_b2a_mt31(; sow = t0 + 9),
            build_b2a_mt32(; sow = t0 + 12),
            build_b2a_mt33(; sow = t0 + 15),
            build_b2a_mt34(; sow = t0 + 18),
            build_b2a_mt40(; sow = t0 + 21),
        ]
        state = decode_b2a_frames(state, messages)

        # Header / time
        @test state.raw_data.last_message_type == 40
        @test state.raw_data.SOW == t0 + 21
        @test state.raw_data.WN == B2A_WN
        @test state.is_shifted_by_180_degrees == false
        # SOW epoch is the start of the newest decoded frame: one frame plus
        # one preamble ago.
        @test state.num_bits_after_valid_syncro_sequence == 624

        # Flags (identical in every message)
        @test state.raw_data.HS == 0
        @test state.raw_data.DIF_B2a == false
        @test state.raw_data.SIF_B2a == false
        @test state.raw_data.AIF_B2a == true
        @test state.raw_data.SISMAI == 5
        @test state.raw_data.DIF_B1C == true
        @test state.raw_data.SIF_B1C == false
        @test state.raw_data.AIF_B1C == false

        # MT10: ephemeris I
        @test state.raw_data.IODE == B2A_IODE
        @test state.raw_data.t_0e == 300000
        @test state.raw_data.sat_type == 3
        @test state.raw_data.ΔA == -12345 * 2.0^-9
        @test state.raw_data.A_dot == 345 * 2.0^-21
        @test state.raw_data.Δn_0 == -2000 * 2.0^-44 * PI
        @test state.raw_data.Δn_0_dot == 12000 * 2.0^-57 * PI
        @test state.raw_data.M_0 == -123456789 * 2.0^-32 * PI
        @test state.raw_data.e == 123456789 * 2.0^-34
        @test state.raw_data.ω == 987654321 * 2.0^-32 * PI
        @test state.raw_data.SOW_mt10 == t0

        # MT11: ephemeris II
        @test state.raw_data.Ω_0 == -1987654321 * 2.0^-32 * PI
        @test state.raw_data.i_0 == 555666777 * 2.0^-32 * PI
        @test state.raw_data.Ω_dot == -90000 * 2.0^-44 * PI
        @test state.raw_data.i_dot == 9000 * 2.0^-44 * PI
        @test state.raw_data.C_is == -2100 * 2.0^-30
        @test state.raw_data.C_ic == 1500 * 2.0^-30
        @test state.raw_data.C_rs == -400000 * 2.0^-8
        @test state.raw_data.C_rc == 300000 * 2.0^-8
        @test state.raw_data.C_us == -800000 * 2.0^-30
        @test state.raw_data.C_uc == 700000 * 2.0^-30
        @test state.raw_data.SOW_mt11 == t0 + 3

        # Clock block (identical in MT30-34)
        @test state.raw_data.IODC == B2A_IODC
        @test state.raw_data.t_0c == 301500
        @test state.raw_data.a_f0 == -5000000 * 2.0^-34
        @test state.raw_data.a_f1 == -150000 * 2.0^-50
        @test state.raw_data.a_f2 == 300 * 2.0^-66

        # MT30: group delay + BDGIM
        @test state.raw_data.T_GD_B2ap == -800 * 2.0^-34
        @test state.raw_data.ISC_B2ad == 600 * 2.0^-34
        @test state.raw_data.T_GD_B1Cp == -1024 * 2.0^-34
        @test state.raw_data.α_bdgim_1 == 800 * 2.0^-3
        @test state.raw_data.α_bdgim_2 == -100 * 2.0^-3
        @test state.raw_data.α_bdgim_3 == 200 * 2.0^-3
        @test state.raw_data.α_bdgim_4 == 90 * 2.0^-3
        @test state.raw_data.α_bdgim_5 == 40 * -2.0^-3
        @test state.raw_data.α_bdgim_6 == -3 * 2.0^-3
        @test state.raw_data.α_bdgim_7 == 17 * 2.0^-3
        @test state.raw_data.α_bdgim_8 == -128 * 2.0^-3
        @test state.raw_data.α_bdgim_9 == 127 * 2.0^-3

        # MT31 + MT33: reduced almanacs (PRN 0 slot skipped)
        @test !isnothing(state.raw_data.reduced_almanacs)
        @test sort(collect(keys(state.raw_data.reduced_almanacs))) == [7, 8, 9]
        alm7 = state.raw_data.reduced_almanacs[7]
        @test alm7.sat_type == 3
        @test alm7.δA == -50 * 2.0^9
        @test alm7.Ω_0 == 31 * 2.0^-6 * PI
        @test alm7.Φ_0 == -64 * 2.0^-6 * PI
        @test alm7.health == 0
        @test alm7.WN_a == B2A_WN
        @test alm7.t_0a == 100 * 2^12
        @test state.raw_data.reduced_almanacs[8].health == 128
        alm9 = state.raw_data.reduced_almanacs[9]
        @test alm9.sat_type == 1
        @test alm9.δA == 100 * 2.0^9
        @test alm9.health == 2

        # MT32: EOP
        @test state.raw_data.t_EOP == 30000 * 2^4
        @test state.raw_data.PM_X == -700000 * 2.0^-20
        @test state.raw_data.PM_X_dot == 8000 * 2.0^-21
        @test state.raw_data.PM_Y == 650000 * 2.0^-20
        @test state.raw_data.PM_Y_dot == -9000 * 2.0^-21
        @test state.raw_data.ΔUT1 == -500000000 * 2.0^-24
        @test state.raw_data.ΔUT1_dot == 140000 * 2.0^-25

        # MT33: BGTO
        @test state.raw_data.GNSS_ID == 1
        @test state.raw_data.WN_0BGTO == 880
        @test state.raw_data.t_0BGTO == 20000 * 2^4
        @test state.raw_data.A_0BGTO == -20000 * 2.0^-35
        @test state.raw_data.A_1BGTO == 3000 * 2.0^-51
        @test state.raw_data.A_2BGTO == -60 * 2.0^-68

        # MT34: SISAIoc + BDT-UTC
        @test state.raw_data.t_op == 1200 * 300
        @test state.raw_data.SISAI_ocb == 17
        @test state.raw_data.SISAI_oc1 == 5
        @test state.raw_data.SISAI_oc2 == 2
        @test state.raw_data.A_0UTC == -30000 * 2.0^-35
        @test state.raw_data.A_1UTC == 2500 * 2.0^-51
        @test state.raw_data.A_2UTC == -50 * 2.0^-68
        @test state.raw_data.Δt_LS == 4
        @test state.raw_data.t_0t == 35000 * 2^4
        @test state.raw_data.WN_0t == B2A_WN
        @test state.raw_data.WN_LSF == 902
        @test state.raw_data.DN == 6
        @test state.raw_data.Δt_LSF == 5

        # MT40: SISAIoe + midi almanac
        @test state.raw_data.SISAI_oe == 11
        @test !isnothing(state.raw_data.midi_almanacs)
        @test collect(keys(state.raw_data.midi_almanacs)) == [23]
        midi = state.raw_data.midi_almanacs[23]
        @test midi.sat_type == 2
        @test midi.WN_a == B2A_WN
        @test midi.t_0a == 100 * 2^12
        @test midi.e == 1500 * 2.0^-16
        @test midi.δi == -800 * 2.0^-14 * PI
        @test midi.sqrt_A == 105000 * 2.0^-4
        @test midi.Ω_0 == -20000 * 2.0^-15 * PI
        @test midi.Ω_dot == -700 * 2.0^-33 * PI
        @test midi.ω == 15000 * 2.0^-15 * PI
        @test midi.M_0 == -15000 * 2.0^-15 * PI
        @test midi.a_f0 == -900 * 2.0^-20
        @test midi.a_f1 == 333 * 2.0^-37
        @test midi.health == 0

        # Validation: complete after MT10 + MT11 + a clock message
        @test is_decoding_completed_for_positioning(state)
        @test state.data == state.raw_data
        @test is_sat_healthy(state)
    end

    @testset "Validation gating" begin
        # MT10 alone: not complete.
        state = BeiDouB2aDecoderState(B2A_PRN)
        state = decode_b2a_frames(state, [build_b2a_mt10(; sow = t0)])
        @test !is_decoding_completed_for_positioning(state)
        @test state.data == BeiDouB2aData()

        # + adjacent MT11: ephemeris paired, but still no clock.
        state = decode_b2a_frames(state, [build_b2a_mt11(; sow = t0 + 3)])
        @test !is_decoding_completed_for_positioning(state)

        # + MT30 with a *mismatched* IODC (IODE ≠ IODC & 0xff): not validated.
        state =
            decode_b2a_frames(state, [build_b2a_mt30(; sow = t0 + 6, iodc = B2A_IODC + 1)])
        @test !is_decoding_completed_for_positioning(state)
        @test state.data == BeiDouB2aData()

        # + MT30 with the matching IODC: validated.
        state = decode_b2a_frames(state, [build_b2a_mt30(; sow = t0 + 9)])
        @test is_decoding_completed_for_positioning(state)
        @test state.data.IODE == state.data.IODC & 0xff

        # A reserved orbit type (0) leaves the semi-major axis unknowable: it
        # is what selects the `A_ref` the broadcast `ΔA` corrects, so the
        # satellite must not be used (ICD Table 7-6).
        @test !GNSSDecoder.is_ephemeris_decoded(BeiDouB2aData(state.data; sat_type = 0))
        @test all(
            GNSSDecoder.is_ephemeris_decoded(BeiDouB2aData(state.data; sat_type = t)) for
            t = 1:3
        )
    end

    @testset "MT10/MT11 pairing requires adjacent frames" begin
        # MT10 and MT11 separated by an MT30 frame: no valid pair.
        state = BeiDouB2aDecoderState(B2A_PRN)
        state = decode_b2a_frames(
            state,
            [
                build_b2a_mt10(; sow = t0),
                build_b2a_mt30(; sow = t0 + 3),
                build_b2a_mt11(; sow = t0 + 6),
            ],
        )
        @test !is_decoding_completed_for_positioning(state)

        # A fresh MT10 right after that MT11 (order reversed) pairs up.
        state = decode_b2a_frames(state, [build_b2a_mt10(; sow = t0 + 9)])
        @test is_decoding_completed_for_positioning(state)
    end

    @testset "MT10/MT11 pairing wraps at the BDT week rollover" begin
        # Adjacency is modular: the last frame of the week (SOW 604797)
        # followed by the first of the next (SOW 0) is a valid pair, while a
        # straddling non-adjacent gap is not.
        @test GNSSDecoder._b2a_sow_adjacent(604797, 0)
        @test GNSSDecoder._b2a_sow_adjacent(0, 604797)
        @test GNSSDecoder._b2a_sow_adjacent(604794, 604797)
        @test !GNSSDecoder._b2a_sow_adjacent(604794, 0)
        @test !GNSSDecoder._b2a_sow_adjacent(604797, 3)

        state = BeiDouB2aDecoderState(B2A_PRN)
        state = decode_b2a_frames(
            state,
            [build_b2a_mt10(; sow = 604797), build_b2a_mt11(; sow = 0)],
        )
        @test GNSSDecoder.is_ephemeris_decoded(state.raw_data)
    end

    @testset "Skipped promotions keep tow + num_bits/rate on the true time" begin
        # `get_time_of_week` reads the *validated* `data.SOW`, so the symbol
        # counter must anchor to the last frame `validate_data` promoted, not to
        # every frame that clears CRC. A fresh ephemeris block breaks first the
        # MT10/MT11 adjacency gate and then the IODE == IODC & 0xff match, so
        # promotion is skipped for two frames; across that gap the published SOW
        # stands still and only the counter may carry the elapsed time. Re-arming
        # it per decoded frame made `now = tow + num_bits/get_data_frequency`
        # (src/gnss.jl) read exactly one 3-second frame low for the whole skip.
        new_iode = B2A_IODE + 1
        new_iodc = (3 << 8) | new_iode  # 790; IODE == IODC & 0xff again
        frames = b2a_frame_symbols.([
            build_b2a_mt10(; sow = t0),
            build_b2a_mt11(; sow = t0 + 3),
            build_b2a_mt30(; sow = t0 + 6),                   # promotes
            build_b2a_mt10(; sow = t0 + 9, iode = new_iode),  # no adjacent MT11 yet
            build_b2a_mt11(; sow = t0 + 12),                  # pairs, but IODC is stale
            build_b2a_mt30(; sow = t0 + 15, iodc = new_iodc), # promotes again
        ])
        # A frame is decoded once the *next* frame's 24 preamble symbols are
        # buffered, so the stream is cut at those sync instants: after the first
        # step each `decode` call below feeds exactly one frame, 600 symbols = 3 s.
        b2a_now(state) =
            get_time_of_week(state.data) + state.num_bits_after_valid_syncro_sequence / 200

        state = BeiDouB2aDecoderState(B2A_PRN)
        head = vcat(frames[1], frames[2], frames[3], frames[4][1:24])
        state = decode(state, head, length(head))
        @test state.data.SOW == t0 + 6
        @test state.num_bits_after_valid_syncro_sequence == 624
        t_promoted = b2a_now(state)

        # Frame 4 decodes but is not promoted: the new MT10 has no partner yet.
        step = vcat(frames[4][25:600], frames[5][1:24])
        state = decode(state, step, length(step))
        @test state.raw_data.SOW == t0 + 9          # the frame did decode ...
        @test state.raw_data.IODE == new_iode
        @test state.data.SOW == t0 + 6              # ... and nothing was promoted
        @test state.num_bits_after_valid_syncro_sequence == 624 + 600
        @test b2a_now(state) ≈ t_promoted + 3 atol = 1e-6

        # Frame 5 pairs the ephemeris again, but the clock set still carries the
        # old IODC, so promotion is skipped a second time.
        step = vcat(frames[5][25:600], frames[6][1:24])
        state = decode(state, step, length(step))
        @test state.raw_data.SOW == t0 + 12
        @test state.data.SOW == t0 + 6
        @test state.num_bits_after_valid_syncro_sequence == 624 + 1200
        @test b2a_now(state) ≈ t_promoted + 6 atol = 1e-6

        # The matching MT30 promotes: SOW and counter re-anchor together, so the
        # reported time carries on ticking without a step.
        step = vcat(frames[6][25:600], b2a_trailing_preamble())
        state = decode(state, step, length(step))
        @test state.data.SOW == t0 + 15
        @test state.data.IODE == new_iode
        @test state.data.IODC == new_iodc
        @test state.num_bits_after_valid_syncro_sequence == 624
        @test b2a_now(state) ≈ t_promoted + 9 atol = 1e-6
    end

    @testset "180-degree phase inversion" begin
        state = BeiDouB2aDecoderState(B2A_PRN)
        state = decode_b2a_frames(
            state,
            [build_b2a_mt10(; sow = t0), build_b2a_mt11(; sow = t0 + 3)];
            invert = true,
        )
        @test state.is_shifted_by_180_degrees == true
        @test state.raw_data.IODE == B2A_IODE
        @test state.raw_data.SOW == t0 + 3
        @test state.raw_data.M_0 == -123456789 * 2.0^-32 * PI
    end

    @testset "Wrong-PRN frames are rejected" begin
        state = BeiDouB2aDecoderState(B2A_PRN)
        state = decode_b2a_frames(state, [build_b2a_mt10(; sow = t0, prn = B2A_PRN + 1)])
        @test isnothing(state.raw_data.IODE)
        @test isnothing(state.raw_data.SOW)
    end

    @testset "CRC-corrupted frames are dropped" begin
        state = BeiDouB2aDecoderState(B2A_PRN)
        state = decode_b2a_frames(state, [build_b2a_mt10(; sow = t0, corrupt_crc = true)])
        @test isnothing(state.raw_data.SOW)
        @test isnothing(state.num_bits_after_valid_syncro_sequence)
    end

    @testset "reset_decoder_state" begin
        state = BeiDouB2aDecoderState(B2A_PRN)
        state = decode_b2a_frames(
            state,
            [
                build_b2a_mt10(; sow = t0),
                build_b2a_mt11(; sow = t0 + 3),
                build_b2a_mt30(; sow = t0 + 6),
            ],
        )
        @test is_decoding_completed_for_positioning(state)
        state = reset_decoder_state(state)
        @test isnothing(state.raw_data.SOW)
        @test state.raw_data.IODE == B2A_IODE   # ephemeris survives reset
        @test state.data == BeiDouB2aData()
        @test isnothing(state.num_bits_after_valid_syncro_sequence)
        @test length(GNSSDecoder.soft_buffer(state)) == 0
        # And the decoder re-locks from scratch.
        state = decode_b2a_frames(
            state,
            [build_b2a_mt11(; sow = t0 + 12), build_b2a_mt10(; sow = t0 + 15)],
        )
        @test is_decoding_completed_for_positioning(state)
    end

    @testset "decode_once stops after the positioning set" begin
        state = BeiDouB2aDecoderState(B2A_PRN)
        messages = [
            build_b2a_mt10(; sow = t0),
            build_b2a_mt11(; sow = t0 + 3),
            build_b2a_mt30(; sow = t0 + 6),
            build_b2a_mt32(; sow = t0 + 9),
        ]
        stream = vcat(reduce(vcat, b2a_frame_symbols.(messages)), b2a_trailing_preamble())
        state = decode(state, stream, length(stream); decode_once = true)
        @test is_decoding_completed_for_positioning(state)
    end

    @testset "Signal metadata accessors" begin
        state = BeiDouB2aDecoderState(B2A_PRN)
        @test get_signal_type(state) == BeiDouB2aI
        @test get_data_frequency(state) == get_data_frequency(BeiDouB2aI)
        @test get_time_system_name(state) == "BeiDou Time"
        @test get_constellation_name(state) == "BeiDou"
        # The pilot component does not map to a decoder; only the data
        # component does.
        @test GNSSDecoderState(BeiDouB2aI(), B2A_PRN) isa GNSSDecoderState{BeiDouB2aData}
    end

    @testset "Unhealthy satellite (HS = 1)" begin
        state = BeiDouB2aDecoderState(B2A_PRN)
        # Same full set but MT11 (and the clock messages) with HS = 1.
        mt11_unhealthy =
            build_b2a_message(; prn = B2A_PRN, mestype = 11, sow_seconds = t0 + 3) do bits
                b2a_setbits!(bits, 31, 2, 1)  # HS = unhealthy
                set_flags!(bits, 33)
                b2a_setbits!(bits, 43, 33, -1987654321)
                b2a_setbits!(bits, 76, 33, 555666777)
                b2a_setbits!(bits, 109, 19, -90000)
                b2a_setbits!(bits, 128, 15, 9000)
                b2a_setbits!(bits, 143, 16, -2100)
                b2a_setbits!(bits, 159, 16, 1500)
                b2a_setbits!(bits, 175, 24, -400000)
                b2a_setbits!(bits, 199, 24, 300000)
                b2a_setbits!(bits, 223, 21, -800000)
                b2a_setbits!(bits, 244, 21, 700000)
            end
        mt30_unhealthy =
            build_b2a_message(; prn = B2A_PRN, mestype = 30, sow_seconds = t0 + 6) do bits
                b2a_setbits!(bits, 31, 2, 1)
                set_flags!(bits, 33)
                set_clock!(bits, 43)
                b2a_setbits!(bits, 112, 10, B2A_IODC)
            end
        state = decode_b2a_frames(
            state,
            [build_b2a_mt10(; sow = t0), mt11_unhealthy, mt30_unhealthy],
        )
        @test is_decoding_completed_for_positioning(state)
        @test state.raw_data.HS == 1
        @test !is_sat_healthy(state)
    end
end
