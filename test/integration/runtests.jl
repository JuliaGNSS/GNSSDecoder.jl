# The real-data integration test, in an environment of its own (the header of
# `flexiband_iii7a.jl` says why). Run it with
#
#     julia --project=test/integration test/integration/runtests.jl

using Test, GNSSDecoder

@testset "GNSSDecoder.jl integration" begin
    include("flexiband_iii7a.jl")
end
