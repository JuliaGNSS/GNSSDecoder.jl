# Preallocated storage for the allocation-free decode path (`decode!`).
#
# A decoder that allocates nothing per symbol must never grow a container
# while decoding: every container a decoder writes into is sized once, when the
# state is constructed, and later writes *overwrite* it in place. This file
# holds the three pieces that make that explicit:
#
#   - `SlotDictionary`, a dictionary with one preallocated slot per possible key,
#     used for every keyed store in the decoded data (almanacs, masks, ...);
#   - `FixedText`, an inline string for the broadcast text messages;
#   - `writable_container` / `overwrite!`, the only two ways decode code touches
#     such a container, plus `duplicate`, the trim-safe deep copy `decode` makes
#     so that it can keep value semantics on top of the overwriting `decode!`.

"""
    SlotIndices

The keys of a [`SlotDictionary`](@ref): the occupied slots, iterated in
ascending key order.
"""
mutable struct SlotIndices <: Dictionaries.AbstractIndices{Int}
    occupied::Vector{Bool}
    count::Int
end

Base.length(indices::SlotIndices) = indices.count
Base.in(key::Integer, indices::SlotIndices) =
    0 <= key < length(indices.occupied) && @inbounds indices.occupied[key+1]
Base.in(key, indices::SlotIndices) = false

function Base.iterate(indices::SlotIndices, slot::Int = 1)
    occupied = indices.occupied
    @inbounds while slot <= length(occupied)
        occupied[slot] && return (slot - 1, slot + 1)
        slot += 1
    end
    return nothing
end

"""
$(TYPEDEF)

A dictionary from `Int` keys in `0:N-1` to values of type `V`, with one slot per
possible key allocated up front — every `N`-sized buffer exists from
construction on, so inserting, overwriting or deleting an entry never
allocates. This is what makes [`decode!`](@ref) allocation-free for the keyed
stores of the decoded data (almanacs keyed by satellite, HAS masks keyed by
Mask ID, ...): the key range of each is small and fixed by its ICD.

It is an `AbstractDictionary` from Dictionaries.jl, so it reads like the
`Dictionary` it replaces (`d[key]`, `haskey`, `keys`, `pairs`, `length`, ...),
with one difference: it iterates in **ascending key order** rather than in
insertion order. Any `AbstractDictionary` converts to it, which is what lets a
decoded data container be constructed from a `Dictionary`.

Inserting a key outside `0:N-1` throws.
"""
mutable struct SlotDictionary{V,N} <: Dictionaries.AbstractDictionary{Int,V}
    # Mutable only for its identity: an immutable wrapper would be re-boxed —
    # allocated — every time it is stored into a `Union{Nothing,SlotDictionary}`
    # field of a rebuilt data container.
    """
    Slot `key + 1` holds the value for `key`; only occupied slots are defined
    """
    const values::Vector{V}
    """
    Occupancy of each slot, shared with (and iterated by) `keys`
    """
    const indices::SlotIndices
end

function SlotDictionary{V,N}() where {V,N}
    SlotDictionary{V,N}(Vector{V}(undef, N), SlotIndices(zeros(Bool, N), 0))
end

function SlotDictionary{V,N}(dict::Dictionaries.AbstractDictionary) where {V,N}
    out = SlotDictionary{V,N}()
    for (key, value) in pairs(dict)
        insert!(out, key, value)
    end
    return out
end

Base.convert(::Type{SlotDictionary{V,N}}, dict::SlotDictionary{V,N}) where {V,N} = dict
Base.convert(
    ::Type{SlotDictionary{V,N}},
    dict::Dictionaries.AbstractDictionary,
) where {V,N} = SlotDictionary{V,N}(dict)

Base.keys(dict::SlotDictionary) = dict.indices
Base.isassigned(dict::SlotDictionary, key::Int) = key in dict.indices

function Base.getindex(dict::SlotDictionary, key::Int)
    key in dict.indices ||
        throw(Dictionaries.IndexError("Dictionary does not contain index: $key"))
    return @inbounds dict.values[key+1]
end

Dictionaries.issettable(::SlotDictionary) = true
Dictionaries.isinsertable(::SlotDictionary) = true

function Base.setindex!(dict::SlotDictionary{V}, value::V, key::Int) where {V}
    key in dict.indices ||
        throw(Dictionaries.IndexError("Dictionary does not contain index: $key"))
    @inbounds dict.values[key+1] = value
    return dict
end

function check_slot(dict::SlotDictionary{V,N}, key::Int) where {V,N}
    0 <= key < N ||
        throw(ArgumentError("key $key is outside the $N slots 0:$(N-1) of this dictionary"))
end

"""
    set!(dict::SlotDictionary, key, value) -> dict

Insert `value` at `key`, overwriting whatever the slot held.
"""
function Dictionaries.set!(dict::SlotDictionary{V}, key::Int, value::V) where {V}
    check_slot(dict, key)
    indices = dict.indices
    @inbounds if !indices.occupied[key+1]
        indices.occupied[key+1] = true
        indices.count += 1
    end
    @inbounds dict.values[key+1] = value
    return dict
end

function Base.insert!(dict::SlotDictionary{V}, key::Int, value::V) where {V}
    key in dict.indices &&
        throw(Dictionaries.IndexError("Dictionary already contains index: $key"))
    return Dictionaries.set!(dict, key, value)
end

function Dictionaries.unset!(dict::SlotDictionary, key::Int)
    indices = dict.indices
    if key in indices
        @inbounds indices.occupied[key+1] = false
        indices.count -= 1
    end
    return dict
end

function Base.delete!(dict::SlotDictionary, key::Int)
    key in dict.indices ||
        throw(Dictionaries.IndexError("Dictionary does not contain index: $key"))
    return Dictionaries.unset!(dict, key)
end

function Base.empty!(dict::SlotDictionary)
    fill!(dict.indices.occupied, false)
    dict.indices.count = 0
    return dict
end

Base.copy(dict::SlotDictionary{V,N}) where {V,N} = SlotDictionary{V,N}(
    copy(dict.values),
    SlotIndices(copy(dict.indices.occupied), length(dict)),
)

"""
$(TYPEDEF)

An inline, immutable ASCII string of at most `N` characters: the broadcast
text messages (GPS CNAV message types 15 and 36, the GPS L1C-D text page),
whose length the ICD fixes. Unlike a `String` it lives inside the decoded data
container instead of on the heap, so decoding a text page allocates nothing.

It is an `AbstractString`, so it compares equal to a `String` of the same
characters and prints as one; `String(text)` converts.
"""
struct FixedText{N} <: AbstractString
    """
    The characters as ASCII code units, `length` of them meaningful
    """
    code_units::NTuple{N,UInt8}
    """
    Number of meaningful code units
    """
    length::Int
end

function FixedText{N}(text::AbstractString) where {N}
    units = codeunits(String(text))
    length(units) <= N || throw(ArgumentError("text longer than $N code units"))
    all(<(0x80), units) || throw(ArgumentError("text must be ASCII"))
    FixedText{N}(ntuple(i -> i <= length(units) ? units[i] : 0x00, Val(N)), length(units))
end

Base.convert(::Type{FixedText{N}}, text::FixedText{N}) where {N} = text
Base.convert(::Type{FixedText{N}}, text::AbstractString) where {N} = FixedText{N}(text)

Base.ncodeunits(text::FixedText) = text.length
Base.codeunit(::FixedText) = UInt8
Base.codeunit(text::FixedText, i::Integer) = text.code_units[i]
Base.isvalid(text::FixedText, i::Integer) = 1 <= i <= text.length

function Base.iterate(text::FixedText, i::Int = 1)
    i > text.length && return nothing
    return (Char(@inbounds text.code_units[i]), i + 1)
end

Base.length(text::FixedText) = text.length

"""
    writable_container(current, spare) -> container

The container a decode step writes a keyed or indexed field into, in place:
the one `raw_data` already holds (`current`), or — the first time the field is
written, while `current` is `nothing` — the preallocated `spare` from the
decoder's cache, emptied first so no entry from an earlier use of it survives.

The returned container is **overwritten**: the caller mutates it and stores it
back into `raw_data`. Any earlier state still referencing it sees the change,
which is the documented contract of [`decode!`](@ref).
"""
function writable_container(current, spare)
    current === nothing || return current
    # `spare` is read from a `DataStorage` field, typed `Union{Nothing,C}` like
    # the data field it backs; it is never `nothing`, and saying so keeps the
    # `clear!` below statically dispatched (Julia 1.10 would not infer it).
    spare === nothing && throw(ArgumentError("preallocated storage is missing"))
    return clear!(spare)
end

clear!(dict::SlotDictionary) = empty!(dict)
clear!(vector::Vector{T}) where {T} = fill!(vector, zero(T))

"""
    overwrite!(dst, src) -> dst

Make `dst` hold exactly what `src` holds, by overwriting `dst` in place — the
explicit copy a decoder makes when it promotes a raw container into the
validated `data`. `dst` and `src` must have the same capacity.
"""
function overwrite!(dst::SlotDictionary{V,N}, src::SlotDictionary{V,N}) where {V,N}
    dst === src && return dst
    copyto!(dst.indices.occupied, src.indices.occupied)
    dst.indices.count = src.indices.count
    @inbounds for slot = 1:N
        if src.indices.occupied[slot]
            dst.values[slot] = src.values[slot]
        end
    end
    return dst
end

overwrite!(dst::Vector{T}, src::Vector{T}) where {T} = copyto!(dst, src)

"""
    publish!(dst, src) -> Union{Nothing,typeof(dst)}

Overwrite the preallocated validated container `dst` with the raw container
`src`, returning `dst` — or `nothing` when `src` is `nothing` (the field was
never decoded). The value to store in the promoted `data`.
"""
function publish!(dst, src)
    src === nothing && return nothing
    dst === nothing && throw(ArgumentError("preallocated storage is missing"))
    return overwrite!(dst, src)
end


# Types whose instances are immutable values or are shared rather than copied.
is_shared_leaf(T) =
    isbitstype(T) ||
    T <: Union{Number,Nothing,Symbol,Enum,String,FixedText} ||
    T <: Union{Aff3ct.ConvViterbiDecoder,Aff3ct.LDPCBPDecoder}

"""
    duplicate_expr(T, ex, depth) -> Expr

Code deep-copying the value `ex` (of type `T`), built at generation time by
recursing over `T`'s field and element types.
"""
function duplicate_expr(T, ex, depth::Int)
    # A recursive type would never bottom out; fall back to a runtime call.
    depth > 16 && return :(duplicate($ex))
    if T isa Union
        value = gensym(:value)
        branches = foldr(Base.uniontypes(T); init = :(error("unreachable"))) do U, rest
            :($value isa $U ? $(duplicate_expr(U, value, depth + 1)) : $rest)
        end
        return :(let $value = $ex
            $branches
        end)
    end
    isconcretetype(T) || return :(duplicate($ex))
    is_shared_leaf(T) && return ex
    T <: SlotDictionary && return :(copy($ex))
    if T <: Base.RefValue
        S = T.parameters[1]
        return :($T($(duplicate_expr(S, :($ex[]), depth + 1))))
    end
    if T <: Array
        E = eltype(T)
        src, dst, i = gensym(:src), gensym(:dst), gensym(:i)
        copy_elements = if isbitstype(E)
            :(copyto!($dst, $src))
        else
            :(for $i in eachindex($src, $dst)
                $dst[$i] = $(duplicate_expr(E, :($src[$i]), depth + 1))
            end)
        end
        keep_capacity = T <: Vector ? :(sizehint!($dst, vector_capacity($src))) : nothing
        return :(let $src = $ex
            $dst = similar($src)
            $copy_elements
            $keep_capacity
            $dst
        end)
    end
    if T <: CircularDeque
        src, dst, v = gensym(:src), gensym(:dst), gensym(:v)
        return :(let $src = $ex
            $dst = $T(capacity($src))
            for $v in $src
                push!($dst, $v)
            end
            $dst
        end)
    end
    if T <: Dictionary
        V = T.parameters[2]
        src, dst, k, v = gensym(:src), gensym(:dst), gensym(:k), gensym(:v)
        return :(let $src = $ex
            $dst = $T()
            for ($k, $v) in pairs($src)
                insert!($dst, $k, $(duplicate_expr(V, v, depth + 1)))
            end
            $dst
        end)
    end
    # A mutable struct from another package may own memory outside Julia (the
    # AFF3CT handles do); copying its fields would alias that memory, so each
    # such type must be handled above explicitly.
    ismutabletype(T) && parentmodule(T) !== @__MODULE__() &&
        return :(throw(ArgumentError($(string("no `duplicate` rule for ", T)))))
    value = gensym(:value)
    fields = [
        duplicate_expr(fieldtype(T, i), :(getfield($value, $i)), depth + 1) for
        i = 1:fieldcount(T)
    ]
    return :(let $value = $ex
        $(Expr(:new, T, fields...))
    end)
end

"""
    duplicate(x)

Deep copy of decoder state for [`decode`](@ref), written out per type instead of
via `Base.deepcopy` (whose reflection `juliac --trim` rejects).

Every container is copied, keeping its capacity, so the copy decodes through
[`decode!`](@ref) without allocating either. The FEC decoder handles — AFF3CT
objects that live outside Julia — are shared rather than copied: they hold no
state between two decodes, but it does mean two copies of a state must not
decode concurrently from different threads.

The copy of a whole type tree is unrolled into one method body at compile time
(see `duplicate_expr`) rather than recursing through `duplicate` itself: a
method that calls itself on ever-different argument types — a HAS message, its
mask, the mask's satellite list — is exactly the recursion inference widens,
which `juliac --trim` then reports as an unresolved call.
"""
# Defined after `duplicate_expr`: a generator may only call functions that
# already exist when the generated function is defined.
@generated duplicate(x::T) where {T} = duplicate_expr(T, :x, 0)

@static if VERSION >= v"1.11"
    vector_capacity(x::Vector) = length(x.ref.mem)
else
    vector_capacity(x::Vector) = length(x)
end

"""
$(TYPEDEF)

The preallocated containers behind a decoder's `raw_data` and `data`: two
instances of the data type `D` with every container field (a
[`SlotDictionary`](@ref) or an `Array`) allocated, which the decode path only
ever writes into.

  - A decode step writes a container field of `raw_data` in place through
    [`writable_container`](@ref)`(raw_data.field, storage.raw.field)`.
  - Promotion to `data` goes through [`publish_data`](@ref), which overwrites
    `storage.validated`'s containers with the raw ones instead of sharing them,
    so `raw_data` keeps being written while `data` holds still.

Each per-signal cache holds one, built by `preallocated_data(D)`.

# Fields

$(TYPEDFIELDS)
"""
struct DataStorage{D}
    """
    Containers `raw_data` is decoded into
    """
    raw::D
    """
    Containers `data` is published into
    """
    validated::D
end

DataStorage{D}() where {D} = DataStorage{D}(preallocated_data(D), preallocated_data(D))

"""
    preallocated_data(D) -> D

An instance of the data type `D` with every container field allocated at the
size the ICD bounds it to — the backing store of a [`DataStorage`](@ref).
Defined per data type.
"""
function preallocated_data end

is_container_type(T) = T <: Union{SlotDictionary,Array}
is_container_field(T) = is_container_type(Base.nonnothingtype(T))

"""
    publish_data(storage::DataStorage{D}, raw::D) -> D

`raw` promoted to validated data: every container field is copied into
(overwriting) the matching container of `storage.validated` via
[`publish!`](@ref), every other field is taken as is. The result shares no
container with `raw`.
"""
@generated function publish_data(storage::DataStorage{D}, raw::D) where {D}
    fields = map(1:fieldcount(D)) do i
        if is_container_field(fieldtype(D, i))
            :(publish!(getfield(storage.validated, $i), getfield(raw, $i)))
        else
            :(getfield(raw, $i))
        end
    end
    return Expr(:new, D, fields...)
end

"""
    @split_nothing (a, b, ...) expr

Evaluate `expr` with each of the local variables `a`, `b`, ... known to be
either `nothing` or not: the expression is duplicated under
`a === nothing ? expr : expr`, so the compiler sees a concrete type for `a` in
each copy.

This is for keyword calls that rebuild a data container from an optional
(`Union{Nothing,T}`) value, e.g. `GPSL1CAData(raw_data; TOW)` with a `TOW` that
may be `nothing`. On Julia 1.10 the keyword `NamedTuple` of such a call has a
non-concrete type, and the keyword sorter then takes its generic, allocating
path; splitting the union at the call site keeps every copy concrete. `n`
variables produce `2^n` copies, so keep the list short.
"""
macro split_nothing(vars, ex)
    names = vars isa Symbol ? [vars] : vars.args
    body = esc(ex)
    for name in reverse(names)
        var = esc(name)
        body = :($var === nothing ? $body : $body)
    end
    return body
end

"""
    push_vote_candidate!(candidates, candidate, capacity) -> candidates

Append `candidate` to a repetition-voting tally in place. The tally never grows
past `capacity` (its preallocated size): when it is full, the candidate with
the fewest votes — the oldest of them on a tie — is overwritten instead.
"""
function push_vote_candidate!(candidates::Vector, candidate, capacity::Int)
    if length(candidates) < capacity
        push!(candidates, candidate)
    else
        weakest = 1
        for i = 2:length(candidates)
            if candidates[i].vote < candidates[weakest].vote
                weakest = i
            end
        end
        candidates[weakest] = candidate
    end
    return candidates
end

"""
    replace_vote_candidates!(candidates, candidate) -> candidates

Overwrite a repetition-voting tally in place so it holds `candidate` alone.
"""
function replace_vote_candidates!(candidates::Vector, candidate)
    empty!(candidates)
    push!(candidates, candidate)
    return candidates
end
