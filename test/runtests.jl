using Test, GNSSDecoder, BitIntegers, ViterbiDecoder, GNSSSignals
using Aqua
using Dictionaries

BitIntegers.@define_integers 4000

"""
    to_soft_symbols(bits, num_bits) -> Vector{Float32}

Convert a packed unsigned integer of hard-decision bits (MSB-first) into a
vector of ±1.0f0 soft symbols. The v2 decoder consumes soft symbols
(`Float32`, positive ⇒ bit 0, negative ⇒ bit 1, magnitude ⇒ confidence);
fixtures in this test suite were captured as hard bits, so the test
boundary uses this helper to feed them into the soft API.
"""
function to_soft_symbols(bits::T, num_bits::Int) where {T<:Unsigned}
    out = Vector{Float32}(undef, num_bits)
    for i in 1:num_bits
        shift = num_bits - i
        bit = (bits >> shift) & T(1)
        out[i] = bit == T(0) ? 1.0f0 : -1.0f0
    end
    out
end

@testset "GNSSDecoder.jl" begin
    @testset "Aqua" begin
        # Aff3ct is now referenced from `src/gps/l1c_d.jl` (issue #38), so it
        # is no longer a stale dependency.
        Aqua.test_all(GNSSDecoder)
    end

    include("bit_fiddling.jl")
    include("gpsl1.jl")
    include("gps_l1c_d.jl")
    include("galileo_e1b.jl")

    # v2 shared-utility deep-module tests (issue #36)
    include("crc.jl")
    include("bch_toi.jl")
    include("deinterleave.jl")
    include("alist.jl")
end
