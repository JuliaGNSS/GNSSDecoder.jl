# Generator for the BeiDou B-CNAV 64-ary LDPC parity-check matrices in AFF3CT
# `.alist` format (binary image).
#
# Source of truth: the four non-binary parity-check matrices transcribed from
# the official ICDs into `beidou_ldpc_coords.jl` (BDS-SIS-ICD-B1C-1.0 §6.2.2,
# BDS-SIS-ICD-B2a-1.0 §6.2.2, BDS-SIS-ICD-B2b-1.0 §6.2.2). Each code is
# defined over GF(2^6) with primitive polynomial p(x) = 1 + x + x^6 and maps
# codeword symbols to bits MSB first (symbol value 1 ⇔ bit vector 000001,
# ICD §6.2.2 of each document).
#
# Multiplication by a fixed non-zero GF(2^6) element is a GF(2)-linear map on
# the 6 coefficient bits, so every GF(2^6) parity check expands into exactly 6
# binary parity checks on the bit stream: the *binary image* of the non-binary
# code. The binary image is the same set of bit sequences the satellite
# transmits, so no approximation enters the code definition. Only the decoder
# changes: binary belief propagation over the image instead of non-binary BP
# over GF(2^6), which buys reuse of the binary machinery (this mirrors how the
# GPS L1C-D LDPC codes are handled via `.alist` files and AFF3CT's binary BP
# decoder, see `generate_alist.jl`).
#
# That substitution costs sensitivity, and this package has not measured how
# much. Binary BP over the image is strictly weaker than non-binary BP over the
# GF(2^6) code, for two reasons that compound: the 6x6 companion-matrix
# expansion of each non-zero entry plants length-4 cycles densely through the
# Tanner graph, exactly where BP's independence assumption is least defensible;
# and the joint constraint binding a symbol's 6 bits is discarded. Both bite
# hardest for short, high-order codes, which is what these are — LDPC(200,100),
# (88,44), (96,48) and (162,81) over GF(64). Expect a loss on the order of a dB
# rather than a fraction of one, showing up as a raised frame-erasure rate at
# low C/N0 once the CRC gate turns a failed decode into a dropped frame.
#
# Quantifying it would need a reference non-binary (FFT-QSPA) decoder to compare
# against, which this repository does not have — `beidou_ldpc_encode` below is an
# encoder only. Until someone writes one and runs the curves, no number here is
# more than a guess.
#
# Layout of the binary image: GF symbol j (0-based) of the codeword occupies
# bit columns 6j+1 .. 6j+6 (1-based), MSB (coefficient of x^5) first —
# exactly the broadcast bit order. Check row i (0-based) of H expands to
# alist rows 6i+1 .. 6i+6. The first 6k columns are the information bits.
#
# Run from anywhere as:
#   julia --project=. --startup-file=no scripts/generate_beidou_alist.jl
#
# Output: `data/bcnv1_sf2.alist`, `data/bcnv1_sf3.alist`, `data/bcnv2.alist`,
# `data/bcnv3.alist`.
#
# The file also defines the GF(2^6) arithmetic and a reference *non-binary*
# systematic encoder (`beidou_ldpc_encode`) used by the test suite to
# synthesize transmit chains and to cross-check the binary expansion against
# an independent encoding path (see `test/beidou_ldpc.jl`).

# Guarded includes so repeated inclusion (several test files pull this script
# in for the GF(2^6) reference encoder) stays warning-free and cheap.
isdefined(@__MODULE__, :BEIDOU_LDPC_CODES) ||
    include(joinpath(@__DIR__, "beidou_ldpc_coords.jl"))
# Reuse the alist serializer shared with the GPS L1C-D generator. Guarded so
# a test run that already included `generate_alist.jl` (test/alist.jl) does
# not re-include it and trigger method-overwrite warnings.
isdefined(@__MODULE__, :write_alist) || include(joinpath(@__DIR__, "generate_alist.jl"))

# ---- GF(2^6) arithmetic, p(x) = 1 + x + x^6 (ICD §6.2.2) --------------------
#
# Elements are integers 0..63; bit b of the integer is the coefficient of x^b
# (so the ICD's MSB-first bit vector [b5 b4 b3 b2 b1 b0] reads the integer's
# bits from bit 5 down to bit 0).

const GF64_POLY = 0b1000011  # x^6 + x + 1

"Carry-less product of two GF(2^6) elements, reduced modulo x^6 + x + 1."
function gf64_mul(a::Integer, b::Integer)
    (0 <= a <= 63 && 0 <= b <= 63) || throw(ArgumentError("GF(64) elements are 0..63"))
    acc = 0
    x = Int(a)
    y = Int(b)
    while y != 0
        if y & 1 == 1
            acc ⊻= x
        end
        y >>= 1
        x <<= 1
        if x & 0b1000000 != 0
            x ⊻= GF64_POLY
        end
    end
    acc
end

"Multiplicative inverse in GF(2^6) (brute force; the field has 63 units)."
function gf64_inv(a::Integer)
    a == 0 && throw(ArgumentError("zero has no inverse"))
    for b in 1:63
        gf64_mul(a, b) == 1 && return b
    end
    error("unreachable")
end

"""
    gf64_mul_matrix(h) -> 6×6 BitMatrix

The GF(2)-linear map "multiply by `h`" on MSB-first coefficient vectors:
`M[r, c]` is the coefficient of `x^(5-r)` in `h · x^(5-c) mod p(x)`, so that
`bits(h * v) = M * bits(v)` over GF(2) with `bits(v)[c] = coeff of x^(5-c)`.
"""
function gf64_mul_matrix(h::Integer)
    M = falses(6, 6)
    for c in 1:6
        prod = gf64_mul(h, 1 << (6 - c))   # h · x^(5-(c-1))
        for r in 1:6
            M[r, c] = (prod >> (6 - r)) & 1 == 1
        end
    end
    M
end

# ---- Non-binary H assembly and reference encoder ----------------------------

"Assemble the dense (n-k) × n GF(2^6) parity-check matrix from its coords."
function assemble_gf64_H(code)
    m = code.n - code.k
    H = zeros(Int, m, code.n)
    for i in 1:m
        for (col, val) in zip(code.index[i], code.element[i])
            H[i, col+1] == 0 || error("duplicate entry in row $i")
            H[i, col+1] = val
        end
    end
    H
end

"""
    gf64_solve(A, B) -> X with A * X = B over GF(2^6)

Gaussian elimination over GF(2^6). `A` must be square and invertible — for
these codes that is the parity part `H[:, k+1:end]`, whose invertibility is
what makes the ICD's systematic form (Annex, "Step 2") well defined.
"""
function gf64_solve(A::Matrix{Int}, B::Matrix{Int})
    m = size(A, 1)
    size(A, 2) == m || throw(ArgumentError("A must be square"))
    size(B, 1) == m || throw(ArgumentError("size mismatch"))
    Ab = hcat(copy(A), copy(B))
    for col in 1:m
        pivot = findnext(!=(0), view(Ab, :, col), col)
        pivot === nothing && error("matrix is singular over GF(2^6)")
        if pivot != col
            Ab[[col, pivot], :] = Ab[[pivot, col], :]
        end
        inv_p = gf64_inv(Ab[col, col])
        for j in col:size(Ab, 2)
            Ab[col, j] = gf64_mul(inv_p, Ab[col, j])
        end
        for r in 1:m
            r == col && continue
            f = Ab[r, col]
            f == 0 && continue
            for j in col:size(Ab, 2)
                Ab[r, j] ⊻= gf64_mul(f, Ab[col, j])
            end
        end
    end
    Ab[:, m+1:end]
end

"""
    beidou_ldpc_encode(code, message::AbstractVector{<:Integer}) -> Vector{Int}

Systematically encode `k` GF(2^6) message symbols into the `n`-symbol
codeword `[message; parity]` such that `H * codeword = 0` over GF(2^6)
(BDS-SIS-ICD-B1C/B2a/B2b Annex: `c = m · G` with `G = [I_k, (H2^{-1} H1)^T]`).
"""
function beidou_ldpc_encode(code, message::AbstractVector{<:Integer})
    length(message) == code.k || throw(ArgumentError("expected $(code.k) symbols"))
    H = assemble_gf64_H(code)
    m = code.n - code.k
    syndrome = zeros(Int, m, 1)
    for i in 1:m
        acc = 0
        for j in 1:code.k
            acc ⊻= gf64_mul(H[i, j], message[j])
        end
        syndrome[i, 1] = acc
    end
    parity = gf64_solve(H[:, code.k+1:end], syndrome)
    vcat(Vector{Int}(message), vec(parity))
end

"GF(2^6) symbols → broadcast bits, MSB first (symbol 1 ⇔ 000001, ICD §6.2.2)."
gf64_symbols_to_bits(symbols) =
    [Bool((s >> (6 - b)) & 1) for s in symbols for b in 1:6]

"Broadcast bits (MSB first) → GF(2^6) symbols."
function gf64_bits_to_symbols(bits)
    length(bits) % 6 == 0 || throw(ArgumentError("bit count must be a multiple of 6"))
    [foldl((acc, b) -> acc << 1 | Int(b), view(bits, 6i-5:6i); init = 0) for
     i in 1:length(bits)÷6]
end

# ---- Binary image ------------------------------------------------------------

"""
    binary_image(code) -> BitMatrix

Expand the GF(2^6) parity-check matrix into its 6(n-k) × 6n binary image:
each non-zero entry `h` becomes the 6×6 GF(2) matrix of "multiply by `h`"
acting on the symbol's MSB-first bit vector. A bit sequence is a codeword of
the binary image iff its symbol sequence is a codeword of the non-binary code.
"""
function binary_image(code)
    m = code.n - code.k
    H = falses(6m, 6 * code.n)
    for i in 1:m
        for (col, val) in zip(code.index[i], code.element[i])
            H[6i-5:6i, 6col+1:6col+6] = gf64_mul_matrix(val)
        end
    end
    H
end

# Cross-check one code end-to-end: the binary image must annihilate the bit
# expansion of every GF(2^6) codeword produced by the independent non-binary
# systematic encoder.
function check_binary_image(code, H_bin; num_random = 4)
    for trial in 1:num_random
        message = rand(0:63, code.k)
        codeword = beidou_ldpc_encode(code, message)
        bits = gf64_symbols_to_bits(codeword)
        all(iszero, (count(H_bin[r, :] .& bits) % 2 for r in 1:size(H_bin, 1))) ||
            error("binary image violated by encoded codeword (trial $trial)")
    end
end

function generate_beidou_alist(data_dir = joinpath(@__DIR__, "..", "data"))
    isdir(data_dir) || mkpath(data_dir)
    for (label, code) in pairs(BEIDOU_LDPC_CODES)
        H_bin = binary_image(code)
        check_binary_image(code, H_bin)
        path = joinpath(data_dir, "$(label).alist")
        write_alist(path, Matrix{Bool}(H_bin))
        @info "Wrote $(label).alist" M = size(H_bin, 1) N = size(H_bin, 2) ones = count(H_bin) path = path
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    generate_beidou_alist()
end
