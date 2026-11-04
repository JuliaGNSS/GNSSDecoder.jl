# Block deinterleaver for GPS L1C-D, plus the matching interleaver (round-trip tests,
# future encoders).
#
# IS-GPS-800G §3.2.3.5: a 38-row by 46-column block interleaver over the 1748-symbol
# payload of subframes 2 + 3 (1200 SF2 + 548 SF3 symbols). The transmitter writes by row
# and reads by column; the receiver writes the 1748 symbols by column and reads by row.
#
# Signal-agnostic: `rows` and `cols` are parameters.

"""
    deinterleave!(dst, src, rows, cols) -> dst

Reverse a `rows × cols` block interleaver. `src` is the received stream; `dst` receives
the original (row-major) stream. `length(dst) == length(src) == rows*cols` is required
and `dst` and `src` may not alias. Works for any element type (soft symbols, `Bool`, ...).
"""
function deinterleave!(
    dst::AbstractVector{T},
    src::AbstractVector{T},
    rows::Int,
    cols::Int,
) where {T}
    n = rows * cols
    # Lengths bound before the checks: Julia 1.10 leaves a throwing branch
    # uninferred, so a `length` call inside the message would dispatch dynamically.
    num_src, num_dst = length(src), length(dst)
    num_src == n || throw(DimensionMismatch("src has $num_src elements, expected $n"))
    num_dst == n || throw(DimensionMismatch("dst has $num_dst elements, expected $n"))
    # Column-major write, row-major read: src[r + (c-1)*rows] -> dst[c + (r-1)*cols].
    @inbounds for r = 1:rows
        for c = 1:cols
            dst[c+(r-1)*cols] = src[r+(c-1)*rows]
        end
    end
    return dst
end

"""
    deinterleave(src, rows, cols) -> Vector

Allocating variant of [`deinterleave!`](@ref).
"""
function deinterleave(src::AbstractVector{T}, rows::Int, cols::Int) where {T}
    dst = similar(src, rows * cols)
    return deinterleave!(dst, src, rows, cols)
end

"""
    interleave!(dst, src, rows, cols) -> dst

Transmit-side block interleaver, inverse of [`deinterleave!`](@ref): writes `src`
row-major into a `rows × cols` matrix and reads column-major into `dst`.
"""
function interleave!(
    dst::AbstractVector{T},
    src::AbstractVector{T},
    rows::Int,
    cols::Int,
) where {T}
    n = rows * cols
    # Lengths bound before the checks: Julia 1.10 leaves a throwing branch
    # uninferred, so a `length` call inside the message would dispatch dynamically.
    num_src, num_dst = length(src), length(dst)
    num_src == n || throw(DimensionMismatch("src has $num_src elements, expected $n"))
    num_dst == n || throw(DimensionMismatch("dst has $num_dst elements, expected $n"))
    @inbounds for r = 1:rows
        for c = 1:cols
            dst[r+(c-1)*rows] = src[c+(r-1)*cols]
        end
    end
    return dst
end

"""
    interleave(src, rows, cols) -> Vector

Allocating variant of [`interleave!`](@ref).
"""
function interleave(src::AbstractVector{T}, rows::Int, cols::Int) where {T}
    dst = similar(src, rows * cols)
    return interleave!(dst, src, rows, cols)
end
