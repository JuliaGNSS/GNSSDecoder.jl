using Test
using Random
import Aff3ct
using GNSSDecoder
using GNSSDecoder: BCH_TOI_CODEWORDS, crc24q, interleave!, GPSL1C_DData

# ---------------------------------------------------------------------------
# Synthetic CNAV-2 frame generator (test-only, mirrors the transmit chain).
#
# A frame is 1800 channel symbols: 52 BCH-encoded TOI symbols (subframe 1)
# followed by the 1748-symbol block-interleaved concatenation of the
# 1200-symbol LDPC-encoded subframe 2 and the 548-symbol LDPC-encoded
# subframe 3. Symbols are emitted as ±1 Float32 (bit 0 ⇒ +1, bit 1 ⇒ -1).
#
# The CNAV-2 LDPC codes are systematic with codeword layout [info | parity]
# (IS-GPS-800G §3.2.3.3) — verified bit-for-bit against a Spirent GSS L1C
# recording. Aff3ct's `LDPCEncoder` derives a *different* (LU-based) info-bit
# mapping from the parity-check matrix, so the transmit chain here encodes
# systematically itself: parse H = [A | B] from the alist and compute
# parity = B⁻¹A·info over GF(2).
# ---------------------------------------------------------------------------

const _REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
const _SF2_ALIST = joinpath(_REPO_ROOT, "data", "cnv2_sf2.alist")
const _SF3_ALIST = joinpath(_REPO_ROOT, "data", "cnv2_sf3.alist")

"""
Parse an alist file into a dense `Bool` parity-check matrix (M×N).
"""
function _alist_H(path::String)
    lines = readlines(path)
    N, M = parse.(Int, split(lines[1]))
    H = falses(M, N)
    # Per-column row indices start at line 5 (after dims, max degrees, and
    # the two degree lists); zero entries are padding.
    for col = 1:N
        for r in parse.(Int, split(lines[4+col]))
            r == 0 && continue
            H[r, col] = true
        end
    end
    return H
end

"""
Precompute the systematic parity generator `P = B⁻¹A` (M×K over GF(2)) from
`H = [A | B]`, so that `codeword = [info; P·info]` satisfies `H·codeword = 0`.
"""
function _systematic_parity_matrix(H::AbstractMatrix{Bool})
    M, N = size(H)
    K = N - M
    aug = Matrix{Bool}(hcat(H[:, (K+1):N], H[:, 1:K]))  # [B | A]
    for col = 1:M
        piv = findfirst(r -> aug[r, col], col:M)
        piv === nothing && error("parity submatrix B is singular at column $col")
        piv += col - 1
        piv != col && (aug[[col, piv], :] = aug[[piv, col], :])
        for r = 1:M
            r != col && aug[r, col] && (aug[r, :] .⊻= aug[col, :])
        end
    end
    return aug[:, (M+1):end]  # B⁻¹A
end

"""
Systematically LDPC-encode `info` bits: returns `[info; parity]` as `Vector{Int}`.
"""
function _ldpc_encode_systematic(parity_matrix::AbstractMatrix{Bool}, info::AbstractVector)
    u = Bool.(info .!= 0)
    parity =
        [reduce(⊻, u[parity_matrix[r, :]]; init = false) for r = 1:size(parity_matrix, 1)]
    return Int.(vcat(u, parity))
end

"""
Write `len` bits of `val` MSB-first into 1-based position `start` of `bits`.
"""
function _setbits!(bits::BitVector, start::Int, len::Int, val::Integer)
    mask = (UInt64(1) << len) - UInt64(1)
    v = UInt64(unsigned(val) & mask)
    @inbounds for i = 0:(len-1)
        bits[start+i] = ((v >> (len - 1 - i)) & UInt64(1)) == UInt64(1)
    end
    return bits
end

"""
Append a 24-bit CRC-24Q (computed over `bits[1:msg_len]`) MSB-first at `bits[msg_len+1 .. msg_len+24]`.
"""
function _append_crc!(bits::BitVector, msg_len::Int)
    crc = crc24q(collect(bits[1:msg_len]))
    @inbounds for i = 0:23
        bits[msg_len+1+i] = ((crc >> (23 - i)) & UInt32(1)) == UInt32(1)
    end
    return bits
end

"""
52 hard symbols (0/1) of the BCH(51,8) TOI codeword for `toi` (first symbol first).
"""
_toi_symbols(toi::Int) = Int[(BCH_TOI_CODEWORDS[toi+1] >> i) & UInt64(1) for i = 0:51]

"""
Encode one frame to ±1 Float32 soft symbols given the TOI and the 1748-symbol interleaved payload.
"""
function _frame_symbols(toi::Int, payload::Vector{Int})
    bits = vcat(_toi_symbols(toi), payload)
    Float32[b == 0 ? 1.0f0 : -1.0f0 for b in bits]
end

@testset "GPS L1C-D (CNAV-2)" begin
    sf2_parity_matrix = _systematic_parity_matrix(_alist_H(_SF2_ALIST))
    sf3_parity_matrix = _systematic_parity_matrix(_alist_H(_SF3_ALIST))
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

    """
    Build a 600-bit subframe-2 info block (576 message bits + 24-bit CRC).
    """
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
        for i = 1:576
            if i in (1:21) ||
               i in (22:32) ||
               i == 33 ||
               i in (34:38) ||
               i in (39:49) ||
               i in (141:173) ||
               i in (174:206)
                continue
            end
            bits[i] = rand(rng, Bool)
        end
        _append_crc!(bits, 576)
        @assert crc24q(collect(bits[1:600])) == 0
        corrupt && (bits[1] = !bits[1])  # break the CRC after it is computed
        return Int32.(collect(bits))
    end

    """
    Build a 274-bit subframe-3 info block (250 message bits + 24-bit CRC).
    """
    function build_sf3_bits()
        bits = falses(274)
        for i = 1:250
            bits[i] = rand(rng, Bool)
        end
        _append_crc!(bits, 250)
        @assert crc24q(collect(bits[1:274])) == 0
        return Int32.(collect(bits))
    end

    """
    LDPC-encode SF2+SF3 info blocks and block-interleave (38×46) into 1748 symbols.
    """
    function build_payload(sf2_info::Vector{Int32}, sf3_info::Vector{Int32})
        x_sf2 = _ldpc_encode_systematic(sf2_parity_matrix, sf2_info)
        x_sf3 = _ldpc_encode_systematic(sf3_parity_matrix, sf3_info)
        @assert length(x_sf2) == 1200
        @assert length(x_sf3) == 548
        src = vcat(x_sf2, x_sf3)
        dst = Vector{Int}(undef, 1748)
        interleave!(dst, src, 38, 46)
        return dst
    end

    """
    Concatenate `n_frames` consecutive frames starting at `toi0`, sharing `payload`.
    """
    function build_stream(toi0::Int, n_frames::Int, payload::Vector{Int})
        stream = Float32[]
        for k = 0:(n_frames-1)
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

    @testset "TOI discontinuity mid-stream resets without crashing" begin
        # A locked frame establishes a decoded TOI (137); a later frame whose
        # TOI does not follow by +1 (200) makes `decode_syncro_sequence` reset
        # the decoder mid-stream, which empties the soft-symbol deque. The
        # generic `decode` loop then runs `drain_after_sync!` unconditionally —
        # before the drain was clamped to `length(deque)` this `popfirst!`'d an
        # empty `CircularDeque` and threw `BoundsError`. This is the exact path
        # hit when re-acquiring a satellite after a signal outage.
        discontinuous = vcat(build_stream(137, 2, payload), build_stream(200, 2, payload))
        state = GPSL1C_DDecoderState(7)
        crashed = false
        try
            state = decode(state, discontinuous, length(discontinuous))
        catch
            crashed = true
        end
        @test !crashed                                     # regression: no BoundsError
        @test state.data == GPSL1C_DData()                 # validated data cleared by reset
        @test state.raw_data.toi === nothing               # in-flight TOI cleared
        @test state.raw_data.WN == golden.WN               # long-lived CED preserved

        # And the decoder recovers: a subsequent run of consecutive frames
        # re-acquires sync and decodes subframe 2 again. (The exact final TOI
        # depends on residual buffer alignment after the reset, so only assert
        # that sync was re-acquired and the golden fields decode.)
        state = decode(state, build_stream(300, 4, payload), 4 * 1800)
        @test state.data.toi !== nothing
        @test state.data.WN == golden.WN
        @test state.data.M_0 ≈ golden.M0_raw * 2.0^-32 * state.constants.PI
    end

    # --- Subframe 3 page-format parsing (IS-GPS-800J §3.5.4) -----------------
    #
    # Hand-pack a 274-bit SF3 info block (250 message bits + 24-bit CRC) with a
    # known page number + golden fields, LDPC-encode it alongside a valid SF2,
    # interleave, and round-trip through `decode`, then assert the decoded
    # `GPSL1C_DData` fields equal the known values with the ICD scale factors.

    """
    Finalise a 274-bit SF3 block: append CRC over bits 1-250 and return Int32 vector.
    """
    function _finish_sf3(bits::BitVector)
        @assert length(bits) == 274
        _append_crc!(bits, 250)
        @assert crc24q(collect(bits[1:274])) == 0
        return Int32.(collect(bits))
    end

    """
    Build a SF3 page with PRN=`prn` and 6-bit page number `page`, then run `fill!` on the bits.
    """
    function build_sf3_page(fill!::Function, prn::Int, page::Int)
        bits = falses(274)
        _setbits!(bits, 1, 8, prn)
        _setbits!(bits, 9, 6, page)
        fill!(bits)
        return _finish_sf3(bits)
    end

    """
    Decode a stream carrying `sf3_page` (with the standard golden SF2) and return state.data.
    """
    function decode_with_sf3(sf3_page::Vector{Int32}; toi0::Int = 120)
        payload = build_payload(sf2_info, sf3_page)
        state = decode(GPSL1C_DDecoderState(7), build_stream(toi0, 4, payload), 4 * 1800)
        return state
    end

    @testset "SF3 page 1 — UTC + iono + ISC" begin
        st = decode_with_sf3(build_sf3_page(7, 1) do b
            _setbits!(b, 15, 16, 1000)    # A0
            _setbits!(b, 31, 13, -500)    # A1
            _setbits!(b, 44, 7, 3)        # A2
            _setbits!(b, 51, 8, 18)       # ΔtLS
            _setbits!(b, 59, 16, 100)     # tot (scale 2^4)
            _setbits!(b, 75, 13, 2200)    # WNot
            _setbits!(b, 88, 13, 2201)    # WNLSF
            _setbits!(b, 101, 4, 6)       # DN
            _setbits!(b, 105, 8, 19)      # ΔtLSF
            _setbits!(b, 113, 8, 12)      # α0
            _setbits!(b, 121, 8, -3)      # α1
            _setbits!(b, 129, 8, 4)       # α2
            _setbits!(b, 137, 8, -1)      # α3
            _setbits!(b, 145, 8, 7)       # β0
            _setbits!(b, 153, 8, -2)      # β1
            _setbits!(b, 161, 8, 5)       # β2
            _setbits!(b, 169, 8, -4)      # β3
            _setbits!(b, 177, 13, 9)      # ISC_L1CA
            _setbits!(b, 190, 13, -9)     # ISC_L2C
            _setbits!(b, 203, 13, 11)     # ISC_L5I5
            _setbits!(b, 216, 13, -11)    # ISC_L5Q5
        end)
        d = st.data
        @test d.A0_UTC ≈ 1000 * 2.0^-35
        @test d.A1_UTC ≈ -500 * 2.0^-51
        @test d.A2_UTC ≈ 3 * 2.0^-68
        @test d.Δt_LS == 18
        @test d.t_ot == 100 * 16
        @test d.WN_ot == 2200
        @test d.WN_LSF == 2201
        @test d.DN == 6
        @test d.Δt_LSF == 19
        @test d.α_0 ≈ 12 * 2.0^-30
        @test d.α_1 ≈ -3 * 2.0^-27
        @test d.α_2 ≈ 4 * 2.0^-24
        @test d.α_3 ≈ -1 * 2.0^-24
        @test d.β_0 ≈ 7 * 2.0^11
        @test d.β_1 ≈ -2 * 2.0^14
        @test d.β_2 ≈ 5 * 2.0^16
        @test d.β_3 ≈ -4 * 2.0^16
        @test d.ISC_L1CA ≈ 9 * 2.0^-35
        @test d.ISC_L2C ≈ -9 * 2.0^-35
        @test d.ISC_L5I5 ≈ 11 * 2.0^-35
        @test d.ISC_L5Q5 ≈ -11 * 2.0^-35
        @test d.num_sf3_pages_received >= 1
    end

    @testset "SF3 page 2 — GGTO + EOP" begin
        # The EOP fields are anchored to Spirent's GSS L1C reference capture.
        # The raw integers fed in below are exactly what Spirent transmitted; the
        # assertions are the physical values Spirent's *own* reference decoder
        # reported.
        #
        # Spirent prints the 31-bit field as Delta_UT1 = UT1−UTC with the
        # pre-Rev-G scale 2⁻²⁴ (-0.25315952301025391 s). IS-GPS-800J Table 3.5-5
        # redefines it as ΔUT_GPS = UT1−GPS with scale 2⁻²³ — exactly twice that.
        st = decode_with_sf3(build_sf3_page(7, 2) do b
            # GGTO — synthetic non-zero values (the capture carries no GGTO data:
            # GNSS ID = 0, all coefficients 0), kept to exercise the GGTO scaling.
            _setbits!(b, 15, 3, 1)             # GGTO_ID = Galileo
            _setbits!(b, 18, 16, 50)           # tGGTO (scale 2^4)
            _setbits!(b, 34, 13, 2100)         # WNGGTO
            _setbits!(b, 47, 16, 800)          # A0GGTO
            _setbits!(b, 63, 13, -200)         # A1GGTO
            _setbits!(b, 76, 7, 2)             # A2GGTO
            # EOP — raw values as transmitted in the Spirent capture (message 2).
            _setbits!(b, 83, 16, 16200)        # tEOP         -> 259200 s
            _setbits!(b, 99, 21, 111281)       # PM_X         -> 0.10612583160400391
            _setbits!(b, 120, 15, 2779)        # PM_X_dot     -> 0.0013251304626464844
            _setbits!(b, 135, 21, 467600)      # PM_Y         -> 0.4459381103515625
            _setbits!(b, 156, 15, 843)         # PM_Y_dot     -> 0.00040197372436523438
            _setbits!(b, 171, 31, -4247312)    # ΔUT_GPS raw  -> Spirent UT1−UTC @2⁻²⁴
            _setbits!(b, 202, 19, -15502)      # ΔUT_GPS_dot  -> -0.00046199560165405273
        end)
        d = st.data
        @test d.GGTO_ID == 1
        @test d.t_GGTO == 50 * 16
        @test d.WN_GGTO == 2100
        @test d.A0_GGTO ≈ 800 * 2.0^-35
        @test d.A1_GGTO ≈ -200 * 2.0^-51
        @test d.A2_GGTO ≈ 2 * 2.0^-68
        # EOP — Spirent reference-decoder values (independent of our scale LSBs).
        @test d.t_EOP == 259200
        @test d.PM_X ≈ 0.10612583160400391
        @test d.PM_X_dot ≈ 0.0013251304626464844
        @test d.PM_Y ≈ 0.4459381103515625
        @test d.PM_Y_dot ≈ 0.00040197372436523438
        # Rev-J ΔUT_GPS = 2 × Spirent's UT1−UTC (2⁻²³ vs Spirent's 2⁻²⁴); a
        # regression to 2⁻²⁴ would reproduce Spirent's value and FAIL here.
        @test d.ΔUT_GPS ≈ 2 * -0.25315952301025391
        @test d.ΔUT_GPS_dot ≈ -0.00046199560165405273
    end

    @testset "SF3 page 3 — reduced almanac (multi-packet)" begin
        PI = GPSL1C_DDecoderState(7).constants.PI
        st = decode_with_sf3(build_sf3_page(7, 3) do b
            _setbits!(b, 15, 13, 2200)    # WNa
            _setbits!(b, 28, 8, 30)       # toa (scale 2^12)
            # Packet 1 (bit 36): PRN 11
            _setbits!(b, 36, 8, 11)
            _setbits!(b, 44, 8, 5)        # δA
            _setbits!(b, 52, 7, -3)       # Ω0
            _setbits!(b, 59, 7, 2)        # Φ0
            b[66] = false;
            b[67] = true;
            b[68] = false  # L1/L2/L5 health
            # Packet 2 (bit 69): PRN 22
            _setbits!(b, 69, 8, 22)
            _setbits!(b, 77, 8, -4)       # δA
            # Packet 3 (bit 102): PRN 0 ⇒ terminates list
            _setbits!(b, 102, 8, 0)
        end)
        d = st.data
        @test !isnothing(d.reduced_almanacs)
        @test haskey(d.reduced_almanacs, 11)
        @test haskey(d.reduced_almanacs, 22)
        @test !haskey(d.reduced_almanacs, 0)
        a = d.reduced_almanacs[11]
        @test a.WN_a == 2200
        @test a.t_oa == 30 * 4096
        @test a.δA ≈ 5 * 2.0^9
        @test a.Ω_0 ≈ -3 * 2.0^-6 * PI
        @test a.Φ_0 ≈ 2 * 2.0^-6 * PI
        @test a.l1_health == false
        @test a.l2_health == true
        @test a.l5_health == false
        @test d.reduced_almanacs[22].δA ≈ -4 * 2.0^9
    end

    @testset "SF3 page 4 — midi almanac" begin
        PI = GPSL1C_DDecoderState(7).constants.PI
        st = decode_with_sf3(build_sf3_page(7, 4) do b
            _setbits!(b, 15, 13, 2150)    # WNa
            _setbits!(b, 28, 8, 40)       # toa (scale 2^12)
            _setbits!(b, 36, 8, 19)       # PRNa
            b[44] = false;
            b[45] = true;
            b[46] = false  # L1/L2/L5 health
            _setbits!(b, 47, 11, 100)     # e
            _setbits!(b, 58, 11, -50)     # δi
            _setbits!(b, 69, 11, -7)      # Ω_dot
            _setbits!(b, 80, 17, 81920)   # √A
            _setbits!(b, 97, 16, 2000)    # Ω0
            _setbits!(b, 113, 16, -1500)  # ω
            _setbits!(b, 129, 16, 12345)  # M0
            _setbits!(b, 145, 11, 33)     # af0
            _setbits!(b, 156, 10, -11)    # af1
        end)
        d = st.data
        @test !isnothing(d.midi_almanacs)
        @test haskey(d.midi_almanacs, 19)
        a = d.midi_almanacs[19]
        @test a.WN_a == 2150
        @test a.t_oa == 40 * 4096
        @test a.l1_health == false
        @test a.l2_health == true
        @test a.l5_health == false
        @test a.e ≈ 100 * 2.0^-16
        @test a.δi ≈ -50 * 2.0^-14 * PI
        @test a.Ω_dot ≈ -7 * 2.0^-33 * PI
        @test a.sqrt_A ≈ 81920 * 2.0^-4
        @test a.Ω_0 ≈ 2000 * 2.0^-15 * PI
        @test a.ω ≈ -1500 * 2.0^-15 * PI
        @test a.M_0 ≈ 12345 * 2.0^-15 * PI
        @test a.a_f0 ≈ 33 * 2.0^-20
        @test a.a_f1 ≈ -11 * 2.0^-37
    end

    @testset "SF3 page 5 — differential correction" begin
        PI = GPSL1C_DDecoderState(7).constants.PI
        st = decode_with_sf3(build_sf3_page(7, 5) do b
            _setbits!(b, 15, 11, 12)      # t_op-D (scale 300)
            _setbits!(b, 26, 11, 24)      # t_OD (scale 300)
            b[37] = false                 # DC data type = CNAV-2
            # CDC segment at bit 38
            _setbits!(b, 38, 8, 19)       # PRN ID
            _setbits!(b, 46, 13, 100)     # δaf0
            _setbits!(b, 59, 8, -10)      # δaf1
            _setbits!(b, 67, 5, 3)        # UDRA
            # EDC segment at bit 72
            _setbits!(b, 72, 8, 19)       # PRN ID (same SV)
            _setbits!(b, 80, 14, 50)      # Δα
            _setbits!(b, 94, 14, -25)     # Δβ
            _setbits!(b, 108, 15, 7)      # Δγ
            _setbits!(b, 123, 12, -3)     # Δi
            _setbits!(b, 135, 12, 6)      # ΔΩ
            _setbits!(b, 147, 12, -8)     # ΔA
            _setbits!(b, 159, 5, -2)      # UDRA-dot
        end)
        d = st.data
        @test !isnothing(d.differential_corrections)
        @test haskey(d.differential_corrections, 19)
        c = d.differential_corrections[19]
        @test c.t_op_D == 12 * 300
        @test c.t_OD == 24 * 300
        @test c.dc_data_type == false
        @test c.δa_f0 ≈ 100 * 2.0^-35
        @test c.δa_f1 ≈ -10 * 2.0^-51
        @test c.UDRA_index == 3
        @test c.Δα ≈ 50 * 2.0^-34
        @test c.Δβ ≈ -25 * 2.0^-34
        @test c.Δγ ≈ 7 * 2.0^-32 * PI
        @test c.Δi ≈ -3 * 2.0^-32 * PI
        @test c.ΔΩ ≈ 6 * 2.0^-32 * PI
        @test c.ΔA ≈ -8 * 2.0^-9
        @test c.UDRA_dot_index == -2
    end

    @testset "SF3 page 6 — text message" begin
        msg = "HELLO L1C-D #39 CNAV2 TEST!!!"  # 29 ASCII characters
        @assert length(msg) == 29
        st = decode_with_sf3(build_sf3_page(7, 6) do b
            for (k, ch) in enumerate(collect(msg))
                _setbits!(b, 19 + 8 * (k - 1), 8, Int(ch))
            end
        end)
        @test st.data.text_message == msg
    end

    @testset "SF3 unknown/reserved page is ignored but still counted" begin
        # Page 7 is reserved (SV config); page 9+ is undefined. Both must be
        # silently ignored — no exception — while the page counter increments.
        st = decode_with_sf3(build_sf3_page(7, 7) do b
            # arbitrary payload bits
            _setbits!(b, 15, 16, 0xABCD)
        end)
        d = st.data
        @test d.num_sf3_pages_received >= 1
        @test isnothing(d.reduced_almanacs)
        @test isnothing(d.midi_almanacs)
        @test isnothing(d.differential_corrections)
        @test isnothing(d.text_message)
        @test isnothing(d.A0_UTC)
        @test isnothing(d.A0_GGTO)
    end

    # --- Spirent-derived recording fixtures (committed; always run) ----------
    #
    # `test/data/gps_l1c_d_prn{1,2}_symbols.bin` hold the post-FEC channel
    # symbols of two satellites (PRN 1 and PRN 2) from a Spirent GSS7000 L1C
    # simulation: 69 consecutive 1800-symbol CNAV-2 frames each (TOI 1..69),
    # packed 8 hard symbols per byte MSB-first (15,525 bytes per file). They
    # were extracted from the simulator's GSS-CNAVDATA output
    # (`nav_data_fec.L1_cnv`, Spirent FAQ14061: 16-byte header + 225-byte
    # post-FEC symbol blocks, one block per satellite channel per 18-second
    # epoch in round-robin order; the full recording carries 31 channels =
    # PRN 1..31) by demultiplexing satellite channels 0 and 1. These are the
    # *encoded* channel symbols, so the test exercises the complete receive
    # chain: BCH TOI sync → 38×46 deinterleave → LDPC BP decode → CRC-24Q →
    # field parsing. The golden field values below come from Spirent's own
    # decode of the matching pre-FEC bits file (`nav_data_bits.L1_cnv.txt`,
    # message 1, PRN 1).
    """
    Unpack a packed-bit fixture file to ±1 Float32 soft symbols (MSB of each byte first).
    """
    function load_packed_symbols(path::String)
        Float32[(b >> (7 - j)) & 0x01 == 0 ? 1.0f0 : -1.0f0 for b in read(path) for j = 0:7]
    end

    """
    Assert the decoded `state` matches Spirent's pre-FEC field dump (PRN 1).
    """
    function assert_spirent_golden(state)
        d = state.data
        # Subframe 2 (semicircle-valued fields are stored in radians, × π).
        @test d.WN == 2106
        @test d.ITOW == 36
        @test d.t_op == 259200
        @test d.l1c_health == false
        @test d.t_0e == 264600
        @test d.ΔA ≈ 922.048828125
        @test d.A_dot ≈ 0.0
        @test d.Δn_0 ≈ 0.0
        @test d.M_0 ≈ 0.95177634549327195 * π
        @test d.e ≈ 0.0
        @test d.ω ≈ -0.38260139431804419 * π
        @test d.Ω_0 ≈ 0.29631945816799998 * π
        @test d.i_0 ≈ 0.30555555550381541 * π
        @test d.ΔΩ_dot ≈ 2.6000179786933586e-9 * π
        @test d.a_f0 ≈ 0.0
        @test d.T_GD ≈ 0.0
        @test is_sat_healthy(state)

        # Subframe 3 pages observed in the recording.
        @test d.Δt_LS == 18
        @test d.A0_UTC ≈ 0.0
        @test d.t_GGTO == 259200
        @test d.text_message == "Test text message for page: 2"
        @test d.num_sf3_pages_received == 68

        # Reduced almanac (page 3): all 31 SVs broadcast; entries golden-
        # checked against Spirent's dump for PRN 1 and PRN 2 (semicircle
        # fields × π).
        @test !isnothing(d.reduced_almanacs) && length(d.reduced_almanacs) == 31
        ra1 = d.reduced_almanacs[1]
        @test ra1.WN_a == 2106
        @test ra1.t_oa == 507904
        @test ra1.δA ≈ 1024.0
        @test ra1.Ω_0 ≈ 0.296875 * π
        @test ra1.Φ_0 ≈ -0.140625 * π
        @test !ra1.l1_health && !ra1.l2_health && !ra1.l5_health
        ra2 = d.reduced_almanacs[2]
        @test ra2.δA ≈ 1024.0
        @test ra2.Ω_0 ≈ -0.0625 * π
        @test ra2.Φ_0 ≈ -0.5625 * π

        # Midi almanac (page 4): PRNs 1-11 broadcast; PRN 1 entry golden-
        # checked against Spirent's dump.
        @test !isnothing(d.midi_almanacs) && length(d.midi_almanacs) == 11
        ma1 = d.midi_almanacs[1]
        @test ma1.WN_a == 2106
        @test ma1.t_oa == 507904
        @test ma1.e ≈ 0.0
        @test ma1.δi ≈ 0.00555419921875 * π
        @test ma1.Ω_dot ≈ 0.0
        @test ma1.sqrt_A ≈ 5153.6875
        @test ma1.Ω_0 ≈ 0.29632568359375 * π
        @test ma1.ω ≈ -0.382598876953125 * π
        @test ma1.M_0 ≈ 0.247406005859375 * π
        @test ma1.a_f0 ≈ 0.0
        @test ma1.a_f1 ≈ 0.0
        @test !ma1.l1_health && !ma1.l2_health && !ma1.l5_health
    end

    @testset "Spirent recording fixture (PRN 1)" begin
        soft = load_packed_symbols(joinpath(@__DIR__, "data", "gps_l1c_d_prn1_symbols.bin"))
        @test length(soft) == 69 * 1800

        state = decode(GPSL1C_DDecoderState(1), soft, length(soft))
        assert_spirent_golden(state)
        @test !state.is_shifted_by_180_degrees
        # The recording starts at TOI 1; the 69th frame cannot complete (its
        # validation needs the next frame's 52 BCH symbols), so the last
        # completed subframe carries TOI 68.
        @test state.data.toi == 68

        # Polarity-inverted stream must decode identically with the
        # 180°-flip flag set.
        state_inv = decode(GPSL1C_DDecoderState(1), -soft, length(soft))
        @test state_inv.is_shifted_by_180_degrees
        @test state_inv.data == state.data
    end

    @testset "Spirent recording fixture (PRN 2)" begin
        soft1 =
            load_packed_symbols(joinpath(@__DIR__, "data", "gps_l1c_d_prn1_symbols.bin"))
        soft2 =
            load_packed_symbols(joinpath(@__DIR__, "data", "gps_l1c_d_prn2_symbols.bin"))
        d1 = decode(GPSL1C_DDecoderState(1), soft1, length(soft1)).data
        state2 = decode(GPSL1C_DDecoderState(2), soft2, length(soft2))
        d2 = state2.data

        # Same constellation epoch, different satellite.
        @test d2.WN == 2106
        @test d2.ITOW == 36
        @test d2.t_0e == 264600
        @test d2.toi == 68
        @test is_sat_healthy(state2)
        @test d2.Ω_0 != d1.Ω_0  # PRN 2 broadcasts its own ephemeris ...
        # ... consistent (at almanac quantisation) with PRN 2's midi-almanac
        # entry independently decoded from PRN 1's stream.
        ma2 = d1.midi_almanacs[2]
        @test isapprox(d2.Ω_0, ma2.Ω_0; atol = 1e-3)
        @test isapprox(d2.ω, ma2.ω; atol = 1e-3)

        # The subframe-3 broadcast (almanacs, UTC, GGTO, text) is common to
        # the whole constellation — both satellites must decode it equally.
        @test d2.reduced_almanacs == d1.reduced_almanacs
        @test d2.midi_almanacs == d1.midi_almanacs
        @test d2.Δt_LS == d1.Δt_LS
        @test d2.t_GGTO == d1.t_GGTO
        @test d2.text_message == d1.text_message
        @test d2.num_sf3_pages_received == 68
    end
end
