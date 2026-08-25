# Shared transmit-chain helpers for the Galileo decoder tests.
#
# Every Galileo data channel — I/NAV on E1-B and E5b-I, F/NAV on E5a-I, C/NAV on
# E6-B — uses the *same* rate-1/2 constraint-length-7 non-systematic
# convolutional code with the G2 output inverted (Galileo OS SIS ICD, Issue 2.2
# §4.1.4; HAS SIS ICD, Issue 1.0 Table 3, which says so explicitly). Only the
# block-interleaver shape, the sync pattern and the payload length differ. So the
# encoder that synthesises an on-air symbol stream from decoded page bits is
# written once here and specialised per signal by its caller, the way
# `cnav_test_utils.jl` serves the GPS L5I and L2C tests.

"""
    galileo_conv_encode(info::Vector{Bool}) -> Vector{Float32}

Transmit-side rate-1/2 K=7 NSC convolutional encoder: G1 = 0o171, G2 = 0o133
with the G2 output inverted. Appends the 6 zero tail bits that terminate the
trellis and returns `2 * (length(info) + 6)` soft symbols in this package's
convention (bit 0 → +1, bit 1 → −1).
"""
function galileo_conv_encode(info::Vector{Bool})
    polys = (0o171, 0o133)
    reg = zeros(Bool, 6)
    out = Float32[]
    for b in vcat(info, zeros(Bool, 6))
        taps = vcat([b], reg)               # [current, last-6]; tap j matches poly MSB
        for (k, p) in enumerate(polys)
            acc = false
            for j = 1:7
                ((p >> (7 - j)) & 1 == 1) && (acc ⊻= taps[j])
            end
            sym = (k == 2) ? !acc : acc     # invert G2
            push!(out, sym ? -1.0f0 : 1.0f0)
        end
        reg = vcat([b], reg[1:(end-1)])
    end
    out
end

"""
    galileo_sync_symbols(pattern) -> Vector{Float32}

Turn a tuple/vector of sync-pattern bits into soft symbols (bit 0 → +1,
bit 1 → −1). The sync pattern is neither encoded nor interleaved on any Galileo
signal.
"""
galileo_sync_symbols(pattern) = Float32[b == 1 ? -1.0f0 : 1.0f0 for b in pattern]

"""
    galileo_push_field!(bits::Vector{Bool}, value, width) -> bits

Append the low `width` bits of `value`, MSB first — the order every Galileo ICD
field is laid out in. Used by the I/NAV and C/NAV page builders to compose page
parts field by field, exactly as the ICD tables read.
"""
function galileo_push_field!(bits::Vector{Bool}, value::Integer, width::Int)
    for i = 1:width
        push!(bits, (UInt64(value) >> (width - i)) & 1 == 1)
    end
    bits
end

"""
    galileo_page_symbols(info::Vector{Bool}, sync, columns) -> Vector{Float32}

One complete on-air page: the sync pattern, then the convolutionally encoded
information bits passed through the signal's block interleaver. Every Galileo
channel interleaves over `GNSSDecoder.GALILEO_INTERLEAVER_ROWS` (= 8) rows, so
`columns` is the only shape parameter: 30 for I/NAV, 61 for E5a, 123 for C/NAV.
"""
galileo_page_symbols(info::Vector{Bool}, pattern, columns::Int) = vcat(
    galileo_sync_symbols(pattern),
    GNSSDecoder.interleave(
        galileo_conv_encode(info),
        columns,
        GNSSDecoder.GALILEO_INTERLEAVER_ROWS,
    ),
)

"""
    galileo_symbol_stream(parts, sync, columns) -> Vector{Float32}

Concatenate page parts into one on-air soft-symbol stream, then append a trailing
sync pattern so the final part's sync window is closed by the next page's sync —
which is what every Galileo decoder's both-ends preamble match requires.
"""
function galileo_symbol_stream(parts::Vector{Vector{Bool}}, pattern, columns::Int)
    stream = Float32[]
    for part in parts
        append!(stream, galileo_page_symbols(part, pattern, columns))
    end
    append!(stream, galileo_sync_symbols(pattern))
    stream
end
