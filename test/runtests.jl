using Test, GNSSDecoder, BitIntegers, GNSSSignals
using Aqua
using Dictionaries
import Aff3ct

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
    for i = 1:num_bits
        shift = num_bits - i
        bit = (bits >> shift) & T(1)
        out[i] = bit == T(0) ? 1.0f0 : -1.0f0
    end
    out
end

@testset "GNSSDecoder.jl" begin
    @testset "Aqua" begin
        # Aff3ct is referenced from `src/` (Galileo E1B's K=7 NSC Viterbi,
        # issue #37, and the GPS L1C-D LDPC decoders, issue #38), so it
        # needs no stale-dep exemption.
        Aqua.test_all(GNSSDecoder)
    end

    include("gnss_supertype.jl")
    include("bit_fiddling.jl")
    include("gpsl1.jl")
    include("gps_l1c_d.jl")
    # Shared CNAV transmit-chain helpers, consumed by both the L5I and L2C tests.
    include("cnav_test_utils.jl")
    include("gps_l5i.jl")
    include("gps_l2cm.jl")
    include("galileo_e1b.jl")
    include("galileo_e5a.jl")

    # v2 shared-utility deep-module tests (issue #36)
    include("crc.jl")
    include("bch_toi.jl")
    include("deinterleave.jl")
    include("alist.jl")

    # Opt-in real-data integration test (Fraunhofer Flexiband III-7a capture, the
    # same recording Tracking.jl uses). Self-gated on
    # ENV["GNSSDECODER_RUN_INTEGRATION_TEST"]; no-ops on a plain `]test`.
    include("flexiband_iii7a.jl")
end
