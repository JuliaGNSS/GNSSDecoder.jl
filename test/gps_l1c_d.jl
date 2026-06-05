using Test
using Random
import Aff3ct
using Aff3ct: LDPCMatrix, LDPCEncoder, encode
using GNSSDecoder
using GNSSDecoder: BCH_TOI_CODEWORDS, crc24q, interleave!, GPSL1C_DData

# ---------------------------------------------------------------------------
# Synthetic CNAV-2 frame generator (test-only, mirrors the transmit chain).
#
# A frame is 1800 channel symbols: 52 BCH-encoded TOI symbols (subframe 1)
# followed by the 1748-symbol block-interleaved concatenation of the
# 1200-symbol LDPC-encoded subframe 2 and the 548-symbol LDPC-encoded
# subframe 3. Symbols are emitted as ±1 Float32 (bit 0 ⇒ +1, bit 1 ⇒ -1).
# ---------------------------------------------------------------------------

const _REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
const _SF2_ALIST = joinpath(_REPO_ROOT, "data", "cnv2_sf2.alist")
const _SF3_ALIST = joinpath(_REPO_ROOT, "data", "cnv2_sf3.alist")

"Write `len` bits of `val` MSB-first into 1-based position `start` of `bits`."
function _setbits!(bits::BitVector, start::Int, len::Int, val::Integer)
    mask = (UInt64(1) << len) - UInt64(1)
    v = UInt64(unsigned(val) & mask)
    @inbounds for i in 0:(len - 1)
        bits[start + i] = ((v >> (len - 1 - i)) & UInt64(1)) == UInt64(1)
    end
    return bits
end

"Append a 24-bit CRC-24Q (computed over `bits[1:msg_len]`) MSB-first at `bits[msg_len+1 .. msg_len+24]`."
function _append_crc!(bits::BitVector, msg_len::Int)
    crc = crc24q(collect(bits[1:msg_len]))
    @inbounds for i in 0:23
        bits[msg_len + 1 + i] = ((crc >> (23 - i)) & UInt32(1)) == UInt32(1)
    end
    return bits
end

"52 hard symbols (0/1) of the BCH(51,8) TOI codeword for `toi` (first symbol first)."
_toi_symbols(toi::Int) = Int[(BCH_TOI_CODEWORDS[toi + 1] >> i) & UInt64(1) for i in 0:51]

"Encode one frame to ±1 Float32 soft symbols given the TOI and the 1748-symbol interleaved payload."
function _frame_symbols(toi::Int, payload::Vector{Int})
    bits = vcat(_toi_symbols(toi), payload)
    Float32[b == 0 ? 1.0f0 : -1.0f0 for b in bits]
end

@testset "GPS L1C-D (CNAV-2)" begin
    sf2_H = LDPCMatrix(_SF2_ALIST)
    sf3_H = LDPCMatrix(_SF3_ALIST)
    sf2_enc = LDPCEncoder(sf2_H)
    sf3_enc = LDPCEncoder(sf3_H)
    rng = MersenneTwister(0x1C0D)

    # --- Hand-packed golden subframe-2 field values (IS-GPS-800G Fig 3.5-1) ---
    golden = (
        WN = 1234,        # 13 bits
        ITOW = 5,         # 8 bits, two-hour epochs
        top_raw = 200,    # 11 bits, scale 300
        health = 0,       # bit 33: 0 = OK
        ura_ed = -3,      # 5 bits, signed
        toe_raw = 100,    # 11 bits, scale 300
        e_raw = 0x0001_2345, # 33 bits, scale 2^-34
        M0_raw = -12_345, # 33 bits, signed, semicircles 2^-32
    )

    "Build a 600-bit subframe-2 info block (576 message bits + 24-bit CRC)."
    function build_sf2_bits(; corrupt::Bool = false)
        bits = falses(600)
        _setbits!(bits, 1, 13, golden.WN)
        _setbits!(bits, 14, 8, golden.ITOW)
        _setbits!(bits, 22, 11, golden.top_raw)
        bits[33] = golden.health == 1
        _setbits!(bits, 34, 5, golden.ura_ed)
        _setbits!(bits, 39, 11, golden.toe_raw)
        _setbits!(bits, 141, 33, golden.M0_raw)
        _setbits!(bits, 174, 33, golden.e_raw)
        # Fill the remaining payload bits pseudo-randomly so the test exercises
        # the full LDPC block, but leave the golden fields fixed.
        for i in 1:576
            if i in (1:21) || i in (22:32) || i == 33 || i in (34:38) ||
               i in (39:49) || i in (141:173) || i in (174:206)
                continue
            end
            bits[i] = rand(rng, Bool)
        end
        _append_crc!(bits, 576)
        @assert crc24q(collect(bits[1:600])) == 0
        corrupt && (bits[1] = !bits[1])  # break the CRC after it is computed
        return Int32.(collect(bits))
    end

    "Build a 274-bit subframe-3 info block (250 message bits + 24-bit CRC)."
    function build_sf3_bits()
        bits = falses(274)
        for i in 1:250
            bits[i] = rand(rng, Bool)
        end
        _append_crc!(bits, 250)
        @assert crc24q(collect(bits[1:274])) == 0
        return Int32.(collect(bits))
    end

    "LDPC-encode SF2+SF3 info blocks and block-interleave (38×46) into 1748 symbols."
    function build_payload(sf2_info::Vector{Int32}, sf3_info::Vector{Int32})
        x_sf2 = encode(sf2_enc, sf2_info)
        x_sf3 = encode(sf3_enc, sf3_info)
        @assert length(x_sf2) == 1200
        @assert length(x_sf3) == 548
        src = vcat(Int.(x_sf2), Int.(x_sf3))
        dst = Vector{Int}(undef, 1748)
        interleave!(dst, src, 38, 46)
        return dst
    end

    "Concatenate `n_frames` consecutive frames starting at `toi0`, sharing `payload`."
    function build_stream(toi0::Int, n_frames::Int, payload::Vector{Int})
        stream = Float32[]
        for k in 0:(n_frames - 1)
            append!(stream, _frame_symbols((toi0 + k) % 400, payload))
        end
        return stream
    end

    sf2_info = build_sf2_bits()
    sf3_info = build_sf3_bits()
    payload = build_payload(sf2_info, sf3_info)

    @testset "Sync, TOI tracking, and subframe-2 golden fields" begin
        toi0 = 137
        stream = build_stream(toi0, 4, payload)
        state = GPSL1C_DDecoderState(7)
        state = decode(state, stream, length(stream))

        # Sync found and TOI tracked monotonically (last validated frame).
        @test state.data.toi == (toi0 + 2) % 400
        @test state.is_shifted_by_180_degrees == false

        d = state.data
        @test d.WN == golden.WN
        @test d.ITOW == golden.ITOW
        @test d.t_op == golden.top_raw * 300
        @test d.l1c_health == false
        @test d.ura_ed_index == golden.ura_ed
        @test d.t_0e == golden.toe_raw * 300
        @test d.t_0c == golden.toe_raw * 300            # CNAV-2 shares one reference time
        @test d.e ≈ golden.e_raw * 2.0^-34
        @test d.M_0 ≈ golden.M0_raw * 2.0^-32 * state.constants.PI
        # Subframe 3 recorded as received pages (no field parsing — issue #39).
        @test d.num_sf3_pages_received >= 1
    end

    @testset "is_sat_healthy reflects the L1C health bit" begin
        # Healthy (bit 33 = 0).
        state = decode(GPSL1C_DDecoderState(7), build_stream(50, 4, payload), 4 * 1800)
        @test is_sat_healthy(state)

        # Unhealthy (bit 33 = 1): rebuild with health flag set.
        bad = falses(600)
        _setbits!(bad, 1, 13, golden.WN)
        _setbits!(bad, 14, 8, golden.ITOW)
        _setbits!(bad, 22, 11, golden.top_raw)
        bad[33] = true
        _setbits!(bad, 34, 5, golden.ura_ed)
        _setbits!(bad, 39, 11, golden.toe_raw)
        _setbits!(bad, 174, 33, golden.e_raw)
        _append_crc!(bad, 576)
        bad_payload = build_payload(Int32.(collect(bad)), sf3_info)
        bad_state =
            decode(GPSL1C_DDecoderState(7), build_stream(60, 4, bad_payload), 4 * 1800)
        @test bad_state.data.l1c_health == true
        @test !is_sat_healthy(bad_state)
    end

    @testset "Polarity-inverted stream decodes with is_shifted_by_180_degrees" begin
        toi0 = 200
        stream = build_stream(toi0, 4, payload)
        state = GPSL1C_DDecoderState(7)
        state = decode(state, -stream, length(stream))  # 180° phase flip ⇒ negate all symbols
        @test state.is_shifted_by_180_degrees == true
        @test state.data.toi == (toi0 + 2) % 400
        @test state.data.WN == golden.WN
        @test state.data.M_0 ≈ golden.M0_raw * 2.0^-32 * state.constants.PI
    end

    @testset "Corrupted SF2 CRC is silently dropped" begin
        sf2_bad = build_sf2_bits(corrupt = true)
        bad_payload = build_payload(sf2_bad, sf3_info)
        state = GPSL1C_DDecoderState(7)
        state = decode(state, build_stream(90, 4, bad_payload), 4 * 1800)
        # SF2 fields never populated; data never validated; no exception thrown.
        @test state.raw_data.WN === nothing
        @test state.data == GPSL1C_DData()
        # SF3 page still counted (CRC of SF3 is independent of the SF2 corruption).
        @test state.raw_data.num_sf3_pages_received >= 1
    end

    @testset "reset_decoder_state clears in-flight state, keeps decoded data" begin
        state = decode(GPSL1C_DDecoderState(7), build_stream(137, 4, payload), 4 * 1800)
        @test state.data.WN == golden.WN

        reset = reset_decoder_state(state)
        @test reset.data == GPSL1C_DData()                 # validated data cleared
        @test reset.raw_data.WN == golden.WN               # long-lived CED preserved
        @test reset.raw_data.toi === nothing               # in-flight TOI cleared
        @test isempty(reset.cache.soft_buffer)             # soft-symbol deque drained
        @test reset.num_bits_after_valid_syncro_sequence === nothing
        @test reset.is_shifted_by_180_degrees == false
    end

    # --- Optional Spirent GSS fixture (gated; not available in CI) -----------
    @testset "Spirent GSS-CNAVDATA fixture (gated)" begin
        fixture_dir = get(ENV, "GPS_L1C_D_FIXTURE_DIR", nothing)
        if isnothing(fixture_dir)
            @info "GPS_L1C_D_FIXTURE_DIR not set — skipping Spirent L1C-D fixture test"
            @test true
        else
            # Spirent GSS streamer writes post-FEC L1C symbols to a file named
            # `nav_data_fec.L1_cnv` (one signed char per channel symbol). Load
            # them as ±1 soft symbols and run the same end-to-end assertions.
            symfile = joinpath(fixture_dir, "nav_data_fec.L1_cnv")
            @test isfile(symfile)
            raw = read(symfile)
            soft = Float32[Float32(reinterpret(Int8, b)) for b in raw]
            state = GPSL1C_DDecoderState(get(ENV, "GPS_L1C_D_FIXTURE_PRN", "1") |> x -> parse(Int, x))
            state = decode(state, soft, length(soft))
            @test !isnothing(state.data.toi)
            @test !isnothing(state.data.WN)
            @test !isnothing(state.data.t_0e)
        end
    end
end
