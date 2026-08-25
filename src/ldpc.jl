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

"""
$(TYPEDEF)

Long-lived scratch for one LDPC block decode: the AFF3CT belief-propagation
decoder plus the three buffers the decode would otherwise allocate per block.

The counterpart of `GalileoViterbiScratch` for the LDPC-coded signals, and held
the same way — one per coded block in the owning decoder's cache (GPS L1C-D and
BeiDou B1C each hold two, for their subframes 2 and 3; B2a and B2b one each).

# Fields

$(TYPEDFIELDS)
"""
struct LDPCScratch
    """
    AFF3CT flooding belief-propagation decoder over a committed `.alist` matrix
    """
    decoder::LDPCBPDecoder
    """
    `N` channel LLRs, so a caller passing a view or a lazily-typed window is
    materialised without allocating
    """
    llr::Vector{Float32}
    """
    `K` decoded information bits as AFF3CT returns them
    """
    info::Vector{Int32}
    """
    the same `K` bits as `Bool`, the form `crc24q` takes
    """
    bits::Vector{Bool}
end

"""
    LDPCScratch(alist_path; num_iterations = 50)

Load a committed `.alist` parity-check matrix, build the decoder for it (see
[`load_ldpc_decoder`](@ref)), and size the three buffers from its `N` and `K`.
"""
function LDPCScratch(alist_path::AbstractString; num_iterations::Integer = 50)
    decoder = load_ldpc_decoder(alist_path; num_iterations)
    LDPCScratch(
        decoder,
        Vector{Float32}(undef, decoder.N),
        Vector{Int32}(undef, decoder.K),
        Vector{Bool}(undef, decoder.K),
    )
end

# Run an Aff3ct LDPC BP decode, CRC-check the info block, and pack it MSB-first
# into a wide word for the shared `get_bits` helpers. CRC failure ⇒ `nothing`
# (the caller silently drops the subframe). `T` is the packed-word type holding
# the info block (e.g. `UInt600` for a GPS L1C-D subframe 2).
"""
Decode, CRC-check, and pack one LDPC info block into a `T`-typed word; `nothing` on CRC failure.

Every buffer comes from `scratch`, so a decode allocates nothing. The belief
propagation itself dominates by orders of magnitude — 630 µs against a few
kilobytes of garbage for a B2b frame — so this is consistency with the Galileo
FEC path rather than a throughput win.
"""
function ldpc_decode_word(scratch::LDPCScratch, symbols, ::Type{T}) where {T}
    # AFF3CT LLR convention matches ours: positive ⇒ bit 0, negative ⇒ bit 1.
    llr = scratch.llr
    length(symbols) == length(llr) || throw(
        DimensionMismatch("expected $(length(llr)) channel LLRs, got $(length(symbols))"),
    )
    copyto!(llr, symbols)
    info = Aff3ct.decode!(scratch.info, scratch.decoder, llr)
    bits = scratch.bits
    @inbounds for i in eachindex(bits)
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
