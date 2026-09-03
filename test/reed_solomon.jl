# Reed-Solomon over GF(2^8) — the primitives behind Galileo HAS's HPVRS outer
# layer (`src/coding/reed_solomon.jl`).
#
# Ground truth is the Galileo HAS SIS ICD itself: Table 42 for the generator
# polynomial and the Annex B attachment for the 255 x 32 systematic generator
# matrix (spot rows in `has_test_vectors.jl`). Everything else is checked by
# field-algebra identities and encode/decode round-trips.

using Random: MersenneTwister, randperm
using GNSSDecoder:
    GALILEO_HAS_GF256,
    rs_generator_polynomial,
    rs_systematic_generator_matrix,
    rs_erasure_decode

@testset "GF(256) field arithmetic" begin
    field = GALILEO_HAS_GF256
    # α is a generator: the 255 powers must be exactly the nonzero elements.
    powers = [field.exponentials[i+1] for i = 0:254]
    @test sort(Int.(powers)) == collect(1:255)
    @test field.exponentials[1] == 0x01
    @test field.exponentials[2] == 0x02
    # p(α) = α^8 + α^4 + α^3 + α^2 + 1 ⟹ α^8 = α^4 + α^3 + α^2 + 1 = 0b00011101
    @test field.exponentials[9] == 0x1d
    # The table wraps at α^255 = α^0 so products can index a raw log sum.
    @test field.exponentials[256] == field.exponentials[1]
    @test field.exponentials[512] == field.exponentials[257]

    @test GNSSDecoder.gf_add(0xa5, 0xa5) == 0x00
    @test GNSSDecoder.gf_multiply(field, 0x00, 0x57) == 0x00
    @test GNSSDecoder.gf_multiply(field, 0x01, 0x57) == 0x57
    # Worked by hand in this field (note: *not* AES's 0x11b). The carry-less
    # product of 0x57 = x^6+x^4+x^2+x+1 and 0x83 = x^7+x+1 is 0x2b79; reducing
    # modulo 0x11d gives 0x2b79 ⊻ (0x11d<<5) = 0x08d9, then
    # 0x08d9 ⊻ (0x11d<<3) = 0x0031.
    @test GNSSDecoder.gf_multiply(field, 0x57, 0x83) == 0x31
    @test GNSSDecoder.gf_multiply(field, 0x83, 0x57) == 0x31
    @test all(
        GNSSDecoder.gf_multiply(field, x, GNSSDecoder.gf_inverse(field, x)) == 0x01 for
        x = 0x01:0xff
    )
    @test_throws DivideError GNSSDecoder.gf_inverse(field, 0x00)
    # Distributivity, on a handful of triples — cheap but catches table slips.
    for (a, b, c) in ((0x03, 0x1f, 0xf0), (0xff, 0x02, 0x7c), (0x11, 0x11, 0x11))
        @test GNSSDecoder.gf_multiply(field, a, GNSSDecoder.gf_add(b, c)) ==
              GNSSDecoder.gf_add(
            GNSSDecoder.gf_multiply(field, a, b),
            GNSSDecoder.gf_multiply(field, a, c),
        )
    end
end

@testset "RS generator polynomial matches HAS SIS ICD Table 42" begin
    field = GALILEO_HAS_GF256
    g = rs_generator_polynomial(field, 223)
    @test length(g) == 224
    @test g[1:14] == HAS_GENERATOR_POLYNOMIAL_HEAD
    # Monic: the leading coefficient g_223 is 1 by construction.
    @test g[224] == HAS_GENERATOR_POLYNOMIAL_TAIL
    # g(α^i) = 0 for every root α^1 … α^223 — the defining property, and an
    # independent check of the coefficients (Horner over GF(256)).
    for i in (1, 2, 57, 222, 223)
        root = field.exponentials[i+1]
        value = 0x00
        for j = length(g):-1:1
            value = GNSSDecoder.gf_add(GNSSDecoder.gf_multiply(field, value, root), g[j])
        end
        @test value == 0x00
    end
    # α^224 is not a root of a narrow-sense code with 223 parity symbols.
    let root = field.exponentials[225], value = 0x00
        for j = length(g):-1:1
            value = GNSSDecoder.gf_add(GNSSDecoder.gf_multiply(field, value, root), g[j])
        end
        @test value != 0x00
    end
    @test_throws ArgumentError rs_generator_polynomial(field, 0)
end

@testset "RS(255, 32) generator matrix matches HAS SIS ICD Annex B" begin
    G = rs_systematic_generator_matrix(GALILEO_HAS_GF256, 255, 32)
    @test size(G) == (255, 32)
    # Systematic: the first 32 rows are the identity (ICD Eq. 15).
    @test G[1:32, :] == UInt8[UInt8(i == j) for i = 1:32, j = 1:32]
    for (row, expected_hex) in HAS_ANNEX_B_ROWS_HEX
        @test G[row, :] == hex2bytes(expected_hex)
    end
    @test_throws ArgumentError rs_systematic_generator_matrix(GALILEO_HAS_GF256, 32, 255)
end

@testset "RS erasure decode round-trips" begin
    field = GALILEO_HAS_GF256
    G = rs_systematic_generator_matrix(field, 255, 32)
    rng = MersenneTwister(20260821)
    for k in (1, 2, 15, 32)
        # Random information block: k rows of J symbols, as HAS lays it out.
        J = 53
        information = rand(rng, UInt8, k, J)
        # Encode: every code-vector row is G[row, 1:k] · information (the
        # trailing 32 - k information symbols are zero, ICD §6.3).
        function encode_row(row)
            out = Vector{UInt8}(undef, J)
            for j = 1:J
                accumulator = 0x00
                for i = 1:k
                    accumulator = GNSSDecoder.gf_add(
                        accumulator,
                        GNSSDecoder.gf_multiply(field, G[row, i], information[i, j]),
                    )
                end
                out[j] = accumulator
            end
            out
        end
        # Any k distinct *transmittable* rows must recover the message — try the
        # systematic rows, an all-parity set, and a random mixture. Rows
        # k+1 … 32 are excluded: for a message shorter than the code dimension
        # they are identically zero in the first k columns, so they carry no
        # information and make the decoding matrix singular (ICD §6.3: "Pages
        # Ck+1, …, CK contain only zeroes and are excluded from transmission").
        transmittable = vcat(1:k, 33:255)
        row_sets = [
            collect(1:k),
            collect(200:(199+k)),
            sort(transmittable[randperm(rng, length(transmittable))[1:k]]),
        ]
        for rows in row_sets
            received = permutedims(reduce(hcat, encode_row.(rows)))
            decoded = rs_erasure_decode(field, G, rows, received, k)
            @test decoded == information
        end
    end
end

@testset "RS erasure decode input validation" begin
    field = GALILEO_HAS_GF256
    G = rs_systematic_generator_matrix(field, 255, 32)
    received = zeros(UInt8, 2, 53)
    @test_throws DimensionMismatch rs_erasure_decode(field, G, [1], received, 2)
    @test_throws DimensionMismatch rs_erasure_decode(field, G, [1, 2, 3], received, 3)
    @test_throws ArgumentError rs_erasure_decode(
        field,
        G,
        collect(1:33),
        zeros(UInt8, 33, 4),
        33,
    )
    @test_throws ArgumentError rs_erasure_decode(field, G, [1, 256], received, 2)
    # Repeating a row makes the decoding matrix singular, which is reported
    # rather than thrown — a receiver must never feed duplicate page IDs.
    @test isnothing(rs_erasure_decode(field, G, [7, 7], received, 2))
    # Nor may it feed a row in k+1 … 32: those code symbols are identically zero
    # for a message of k < 32 blocks, so they are never transmitted and cannot
    # contribute to the decode (ICD §6.3).
    @test isnothing(rs_erasure_decode(field, G, [1, 30], received, 2))
end

@testset "GF(256) matrix inverse" begin
    field = GALILEO_HAS_GF256
    A = UInt8[0x02 0x03; 0x01 0x01]
    Ainv = GNSSDecoder.gf_invert(field, A)
    @test !isnothing(Ainv)
    product = [
        reduce(
            GNSSDecoder.gf_add,
            (GNSSDecoder.gf_multiply(field, A[i, t], Ainv[t, j]) for t = 1:2),
        ) for i = 1:2, j = 1:2
    ]
    @test product == UInt8[0x01 0x00; 0x00 0x01]
    @test isnothing(GNSSDecoder.gf_invert(field, UInt8[0x01 0x01; 0x01 0x01]))
    @test_throws DimensionMismatch GNSSDecoder.gf_invert(field, zeros(UInt8, 2, 3))
end
