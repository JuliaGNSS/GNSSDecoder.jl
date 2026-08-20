# Shared LDPC decode helper used by every LDPC-coded signal decoder: GPS
# L1C-D (CNAV-2 subframes 2 and 3, IS-GPS-800G) and the BeiDou B-CNAV family
# (B-CNAV1 subframes 2 and 3, B-CNAV2, B-CNAV3 — all binary images of the
# ICDs' 64-ary LDPC codes, see `scripts/generate_beidou_alist.jl`). All of
# them pair a rate-1/2 systematic LDPC block with a trailing CRC-24Q inside
# the info block, so the decode → CRC gate → MSB-first packing pipeline is
# one shared function.
#
# Included after `crc.jl` (consumes `crc24q`) and before the per-signal
# decoders that call it.

import Aff3ct
using Aff3ct: LDPCMatrix, LDPCBPDecoder

# Run an Aff3ct LDPC BP decode, CRC-check the info block, and pack it MSB-first
# into a wide word for the shared `get_bits` helpers. CRC failure ⇒ `nothing`
# (the caller silently drops the subframe). `T` is the packed-word type holding
# the `info_bits`-long block (e.g. `UInt600` for a GPS L1C-D subframe 2).
"""
Decode, CRC-check, and pack one LDPC info block into a `T`-typed word; `nothing` on CRC failure.
"""
function ldpc_decode_word(
    decoder::LDPCBPDecoder,
    symbols,
    info_bits::Int,
    ::Type{T},
) where {T}
    # AFF3CT LLR convention matches ours: positive ⇒ bit 0, negative ⇒ bit 1.
    llr = collect(Float32, symbols)
    info = Aff3ct.decode(decoder, llr)
    bits = Vector{Bool}(undef, info_bits)
    @inbounds for i = 1:info_bits
        bits[i] = info[i] != 0
    end
    # CRC-24Q over the whole info block (message bits + trailing 24-bit CRC) is
    # 0 iff the checksum matches; check on the bit vector before packing.
    crc24q(bits) == 0 || return nothing
    # Pack MSB-first so bit 1 is the most-significant bit and bit `info_bits`
    # the least-significant (right-aligned), matching `word_length = info_bits`.
    word = T(0)
    @inbounds for b in bits
        word = (word << 1) | T(b ? 1 : 0)
    end
    return word
end

"""
    load_ldpc_decoder(alist_path; num_iterations = 50) -> LDPCBPDecoder

Load a committed `.alist` parity-check matrix and build an Aff3ct flooding
belief-propagation decoder for it. All the alist artefacts in `data/` lay the
codeword out as `[info | parity]` (the ICDs' systematic layout), but Aff3ct's
alist loader derives its own info-bit positions by Gaussian elimination and
does not generally pick the first K columns, so the ICD layout is forced
before the decoder captures the positions.
"""
function load_ldpc_decoder(alist_path::AbstractString; num_iterations::Integer = 50)
    H = LDPCMatrix(String(alist_path))
    H.info_bits_pos = collect(UInt32, 0:(H.K-1))
    LDPCBPDecoder(H; num_iterations = num_iterations)
end
