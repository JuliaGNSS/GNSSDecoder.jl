# Generator for GPS L1C-D CNAV-2 LDPC parity-check matrices in AFF3CT `.alist`
# format.
#
# Source of truth: IS-GPS-800G Tables 6.2-2 .. 6.2-13. Each table lists the
# (row, col) coordinates of the 1s in one of six submatrices that compose the
# parity-check matrix:
#
#     H = [ A  B  T ]
#         [ C  D  E ]
#
# Those coordinates are committed alongside this script in
# `cnv2_ldpc_coords.jl`, so the generator is fully self-contained — no PocketSDR
# (or any other external file) is needed to reproduce the `.alist` artefacts.
#
# Dimensions per IS-GPS-800G §3.2.3.4:
#   - Subframe 2: m = 600 rows, n = 1200 cols, gap g = 1
#   - Subframe 3: m = 274 rows, n =  548 cols, gap g = 1
#
# AFF3CT `.alist` format (from libldpc / AFF3CT spec):
#   line 1:  "N M"                (N = #cols, M = #rows)
#   line 2:  "dmax_col dmax_row"  (max degrees, by column and row)
#   line 3:  per-column degrees (N integers)
#   line 4:  per-row    degrees (M integers)
#   lines 5 .. 5+N-1:  per column, 1-based row indices of the 1s
#                      (padded with 0s up to dmax_col)
#   then  M lines:     per row, 1-based column indices of the 1s
#                      (padded with 0s up to dmax_row)
#
# Run from anywhere as:
#   julia --project=. --startup-file=no scripts/generate_alist.jl
#
# Output: `data/cnv2_sf2.alist`, `data/cnv2_sf3.alist`.

include(joinpath(@__DIR__, "cnv2_ldpc_coords.jl"))

# Assemble the binary H matrix from its (A, B, C, D, E, T) submatrices, matching
# PocketSDR's `gen_B_LDPC_H` block-insertion rule exactly (1-based here):
#     A: (r,        c)
#     B: (r,        m+c)        [g=1, c=1 only]
#     C: (m-g+r,    c)
#     D: (m-g+r,    m+c)        [g=1, c=1 only]
#     E: (m-g+r,    m+g+c)
#     T: (r,        m+g+c)
# The +1 ↔ +g distinction is carried faithfully so a future g>1 signal plugs in
# the same way.
function assemble_H(m::Int, n::Int, g::Int, A, B, C, D, E, T)
    H = falses(m, n)
    for (r, c) in A
        H[r, c] = true
    end
    for (r, c) in B
        H[r, m+c] = true
    end
    for (r, c) in C
        H[m-g+r, c] = true
    end
    for (r, c) in D
        H[m-g+r, m+c] = true
    end
    for (r, c) in E
        H[m-g+r, m+g+c] = true
    end
    for (r, c) in T
        H[r, m+g+c] = true
    end
    return H
end

# Serialise a binary parity-check matrix in AFF3CT alist format.
function write_alist(path::AbstractString, H::AbstractMatrix{Bool})
    M, N = size(H)
    cols = [findall(view(H, :, j)) for j in 1:N]
    rows = [findall(view(H, i, :)) for i in 1:M]
    col_deg = map(length, cols)
    row_deg = map(length, rows)
    dmax_col = maximum(col_deg)
    dmax_row = maximum(row_deg)
    open(path, "w") do io
        println(io, "$N $M")
        println(io, "$dmax_col $dmax_row")
        println(io, join(col_deg, " "))
        println(io, join(row_deg, " "))
        for col_idx in cols
            padded = vcat(col_idx, fill(0, dmax_col - length(col_idx)))
            println(io, join(padded, " "))
        end
        for row_idx in rows
            padded = vcat(row_idx, fill(0, dmax_row - length(row_idx)))
            println(io, join(padded, " "))
        end
    end
    return path
end

# Cross-check that the coordinates round-trip through `assemble_H` without
# dropping or duplicating any 1.
function check_self_consistency(H, expected_ones)
    actual = count(H)
    actual == expected_ones ||
        error("H matrix has $actual ones, expected $expected_ones (duplicate (r,c)?)")
end

function generate_alist(data_dir = joinpath(@__DIR__, "..", "data"))
    isdir(data_dir) || mkpath(data_dir)
    for (label, fname) in (("SF2", "cnv2_sf2.alist"), ("SF3", "cnv2_sf3.alist"))
        d = getproperty(CNV2_LDPC_COORDS, Symbol(label))
        expected_ones =
            length(d.A) + length(d.B) + length(d.C) + length(d.D) + length(d.E) + length(d.T)
        H = assemble_H(d.m, d.n, d.g, d.A, d.B, d.C, d.D, d.E, d.T)
        check_self_consistency(H, expected_ones)
        path = joinpath(data_dir, fname)
        write_alist(path, H)
        @info "Wrote $(fname)" M = d.m N = d.n ones = expected_ones path = path
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    generate_alist()
end
