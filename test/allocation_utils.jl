"""
    decode_allocations(make_state, symbols) -> NamedTuple

Bytes `decode!` allocates when it runs `symbols` through a decoder state from
`make_state()`, measured three ways, once the call has been compiled:

  - `fresh`: on a newly constructed state, so every first write into the
    preallocated storage (first almanac, first promotion to `data`, ...) runs;
  - `warm`: the same stream a second time on the state the first pass returned,
    so every overwrite of an already-populated store runs;
  - `reset`: `reset_decoder_state!` on that state, then the stream once more.

All three must be zero for the decoder to be allocation-free. Each count is
taken inside this function (not at global scope), so it measures `decode!`
alone.
"""
function decode_allocations(make_state, symbols::Vector{Float32})
    n = length(symbols)
    # Compile every path first, including the reset.
    warmup = decode!(make_state(), symbols, n)
    warmup = decode!(reset_decoder_state!(decode!(warmup, symbols, n)), symbols, n)

    state = make_state()
    fresh = @allocated(state = decode!(state, symbols, n))
    warm = @allocated(state = decode!(state, symbols, n))
    reset = @allocated(state = decode!(reset_decoder_state!(state), symbols, n))
    return (; fresh, warm, reset, state)
end
