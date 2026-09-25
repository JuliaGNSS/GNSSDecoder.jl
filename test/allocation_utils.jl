"""
Whether the zero-allocation assertions are checked on this Julia version.

`decode!` is allocation-free from Julia 1.11 on. Julia 1.10 compiles some keyword
calls with `Union`-typed values through the keyword sorter's generic path, which
allocates; the decoders work around the cases found so far, but a stray
allocation on 1.10 is accepted, so there the assertions are skipped (the tests
still run the decode and check its result).
"""
const CHECK_ALLOCATIONS = VERSION >= v"1.11"

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

"""
    copy_decode_allocations(make_state, symbols) -> Int

Bytes `decode!` allocates on a `copy` of a state that has already decoded
`symbols` once: `copy` must keep every buffer's capacity (the voting tallies
grow within it), so the copy stays allocation-free too.
"""
function copy_decode_allocations(make_state, symbols::Vector{Float32})
    n = length(symbols)
    state = decode!(make_state(), symbols, n)
    decode!(reset_decoder_state!(copy(state)), symbols, n)
    copied = reset_decoder_state!(copy(state))
    return @allocated decode!(copied, symbols, n)
end
