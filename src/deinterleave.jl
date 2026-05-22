# Block deinterleaver for GPS L1C-D — and a matching interleaver, kept
# next to it for round-trip testing and for future encoders.
#
# IS-GPS-800G §3.2.3.5 specifies a 38-row by 46-column block interleaver
# applied to the 1748-symbol payload of subframes 2 + 3 (1200 SF2 + 548
# SF3 channel symbols). At the transmitter the row-major pre-interleaver
# stream is written *by row* into a 38×46 matrix and read out *by column*.
# At the receiver we reverse that: write 1748 received symbols into the
# same shaped matrix *by column*, read out *by row*.
#
# The functions below are signal-agnostic — `rows` and `cols` are
# parameters — so a future signal that uses a different block size (e.g.
# IRN L1 SPS) can reuse the same primitive.

"""
    deinterleave!(dst, src, rows, cols) -> dst

Reverse a `rows × cols` block interleaver. `src` is the received (column-
major-written, row-major-read) stream; `dst` receives the original
(row-major) stream. Both must satisfy `length(dst) == length(src) ==
rows*cols`. `dst` and `src` may not alias.

Element type `T` is preserved — works for `Float32` soft symbols, `Bool`
hard slices, `Int8`, anything.
"""
function deinterleave!(
    dst::AbstractVector{T},
    src::AbstractVector{T},
    rows::Int,
    cols::Int,
) where {T}
    n = rows * cols
    length(src) == n || throw(DimensionMismatch("src has $(length(src)) elements, expected $n"))
    length(dst) == n || throw(DimensionMismatch("dst has $(length(dst)) elements, expected $n"))
    # Column-major write, row-major read: src[r + (c-1)*rows] -> dst[c + (r-1)*cols].
    @inbounds for r in 1:rows
        for c in 1:cols
            dst[c + (r - 1) * cols] = src[r + (c - 1) * rows]
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

Forward (transmit-side) block interleaver. Inverse of
[`deinterleave!`](@ref): write `src` row-major into a `rows × cols`
matrix, read column-major into `dst`.
"""
function interleave!(
    dst::AbstractVector{T},
    src::AbstractVector{T},
    rows::Int,
    cols::Int,
) where {T}
    n = rows * cols
    length(src) == n || throw(DimensionMismatch("src has $(length(src)) elements, expected $n"))
    length(dst) == n || throw(DimensionMismatch("dst has $(length(dst)) elements, expected $n"))
    @inbounds for r in 1:rows
        for c in 1:cols
            dst[r + (c - 1) * rows] = src[c + (r - 1) * cols]
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
