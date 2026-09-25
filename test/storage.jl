using Test, GNSSDecoder, Dictionaries
using GNSSDecoder:
    writable_container,
    overwrite!,
    publish!,
    duplicate,
    push_vote_candidate!,
    replace_vote_candidates!,
    @split_nothing

@testset "Preallocated storage" begin
    @testset "SlotDictionary" begin
        d = SlotDictionary{Float64,8}()
        @test isempty(d)
        @test length(d) == 0
        set!(d, 5, 1.5)
        insert!(d, 2, 2.5)
        @test length(d) == 2
        @test collect(keys(d)) == [2, 5]  # ascending key order, not insertion order
        @test collect(values(d)) == [2.5, 1.5]
        @test d[5] == 1.5
        @test haskey(d, 2) && !haskey(d, 3) && !haskey(d, 42) && !haskey(d, -1)
        @test get(d, 3, 0.0) == 0.0
        d[5] = 3.5
        @test d[5] == 3.5
        set!(d, 5, 4.5)
        @test d[5] == 4.5 && length(d) == 2
        @test_throws Dictionaries.IndexError d[3]
        @test_throws Dictionaries.IndexError d[3] = 1.0
        @test_throws Dictionaries.IndexError insert!(d, 2, 0.0)
        @test_throws ArgumentError set!(d, 8, 0.0)
        @test_throws ArgumentError set!(d, -1, 0.0)
        unset!(d, 7)  # absent: no-op
        @test length(d) == 2
        delete!(d, 2)
        @test collect(keys(d)) == [5]
        @test_throws Dictionaries.IndexError delete!(d, 2)
        empty!(d)
        @test isempty(d)

        # Converts from any `AbstractDictionary`, and compares equal to a
        # `Dictionary` with the same keys in the same (ascending) order.
        dict = Dictionary([1, 4], [10, 40])
        slots = convert(SlotDictionary{Int,8}, dict)
        @test slots isa SlotDictionary{Int,8}
        @test slots == dict
        @test convert(SlotDictionary{Int,8}, slots) === slots

        # `copy` is independent.
        c = copy(slots)
        set!(c, 6, 60)
        @test length(c) == 3 && length(slots) == 2

        # Writes into a preallocated dictionary allocate nothing.
        fill_slots!(d) = (for k = 0:7
            set!(d, k, Float64(k))
        end;
        d)
        fill_slots!(d)
        @test @allocated(fill_slots!(d)) == 0
    end

    @testset "FixedText" begin
        text = FixedText{8}("abc")
        @test text == "abc"
        @test "abc" == text
        @test length(text) == 3
        @test ncodeunits(text) == 3
        @test String(text) == "abc"
        @test sprint(print, text) == "abc"
        @test collect(text) == ['a', 'b', 'c']
        @test convert(FixedText{8}, "xy") == "xy"
        @test convert(FixedText{8}, text) === text
        @test isbits(text)
        @test_throws ArgumentError FixedText{2}("abc")
        @test_throws ArgumentError FixedText{8}("é")
    end

    @testset "writable_container, overwrite! and publish!" begin
        spare = SlotDictionary{Int,4}()
        set!(spare, 1, 1)
        # First write: the spare, emptied.
        @test writable_container(nothing, spare) === spare
        @test isempty(spare)
        # Later writes: the container already in use.
        current = SlotDictionary{Int,4}()
        @test writable_container(current, spare) === current
        @test_throws ArgumentError writable_container(nothing, nothing)

        vec_spare = [1, 2, 3]
        @test writable_container(nothing, vec_spare) == [0, 0, 0]

        src = SlotDictionary{Int,4}()
        set!(src, 3, 30)
        dst = SlotDictionary{Int,4}()
        set!(dst, 0, 1)
        @test overwrite!(dst, src) === dst
        @test dst == src && dst !== src
        @test overwrite!(src, src) === src
        @test overwrite!([0, 0], [1, 2]) == [1, 2]

        @test isnothing(publish!(dst, nothing))
        @test publish!(dst, src) === dst
        @test_throws ArgumentError publish!(nothing, src)
    end

    @testset "voting tally" begin
        Candidate = NamedTuple{(:vote, :id),Tuple{Int,Int}}
        tally = sizehint!(Candidate[], 3)
        push_vote_candidate!(tally, (vote = 2, id = 1), 3)
        push_vote_candidate!(tally, (vote = 0, id = 2), 3)
        push_vote_candidate!(tally, (vote = 1, id = 3), 3)
        # Full: the weakest candidate is overwritten.
        push_vote_candidate!(tally, (vote = 0, id = 4), 3)
        @test [c.id for c in tally] == [1, 4, 3]
        replace_vote_candidates!(tally, (vote = 5, id = 5))
        @test tally == [(vote = 5, id = 5)]
    end

    @testset "@split_nothing" begin
        describe(x, y) = (x, y)
        for a in (nothing, 1), b in (nothing, 2.0)
            @test (@split_nothing (a, b) describe(a, b)) === (a, b)
        end
        c = nothing
        @test (@split_nothing c describe(c, c)) === (nothing, nothing)
    end

    @testset "duplicate and copy" begin
        d = Dictionary([1], [[1, 2]])
        dd = duplicate(d)
        @test dd == d
        dd[1][1] = 99
        @test d[1] == [1, 2]  # deep: element vectors are copied too
        m = [1 2; 3 4]
        @test duplicate(m) == m && duplicate(m) !== m
        nested = [[1], [2, 3]]
        nd = duplicate(nested)
        @test nd == nested && nd[1] !== nested[1]
        r = Ref{Union{Nothing,Vector{Int}}}([1])
        rd = duplicate(r)
        @test rd[] == [1] && rd[] !== r[]

        # `copy` of a decoder state shares nothing mutable with the original.
        state = GPSL1CADecoderState(1)
        GNSSDecoder.push_soft_symbol!(state, 1.0f0)
        copied = copy(state)
        @test copied == state
        GNSSDecoder.push_soft_symbol!(copied, -1.0f0)
        @test GNSSDecoder.num_bits_buffered(state) == 1
        @test GNSSDecoder.num_bits_buffered(copied) == 2
        @test copied.cache.storage.raw.almanacs !== state.cache.storage.raw.almanacs
    end
end
