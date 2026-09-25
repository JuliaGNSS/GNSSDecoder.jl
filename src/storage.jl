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
Base.convert(::Type{SlotDictionary{V,N}}, dict::Dictionaries.AbstractDictionary) where {V,N} =
    SlotDictionary{V,N}(dict)

Base.keys(dict::SlotDictionary) = dict.indices
Base.isassigned(dict::SlotDictionary, key::Int) = key in dict.indices

function Base.getindex(dict::SlotDictionary, key::Int)
    key in dict.indices || throw(Dictionaries.IndexError("Dictionary does not contain index: $key"))
    return @inbounds dict.values[key+1]
end

Dictionaries.issettable(::SlotDictionary) = true
Dictionaries.isinsertable(::SlotDictionary) = true

function Base.setindex!(dict::SlotDictionary{V}, value::V, key::Int) where {V}
    key in dict.indices || throw(Dictionaries.IndexError("Dictionary does not contain index: $key"))
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
    key in dict.indices && throw(Dictionaries.IndexError("Dictionary already contains index: $key"))
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
    key in dict.indices || throw(Dictionaries.IndexError("Dictionary does not contain index: $key"))
    return Dictionaries.unset!(dict, key)
end

function Base.empty!(dict::SlotDictionary)
    fill!(dict.indices.occupied, false)
    dict.indices.count = 0
    return dict
end

Base.copy(dict::SlotDictionary{V,N}) where {V,N} =
    SlotDictionary{V,N}(copy(dict.values), SlotIndices(copy(dict.indices.occupied), length(dict)))

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
writable_container(current, spare) = isnothing(current) ? clear!(spare) : current

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
publish!(dst, src) = isnothing(src) ? nothing : overwrite!(dst, src)

"""
    duplicate(x)

Deep copy of decoder state for [`decode`](@ref), written out per type instead of
via `Base.deepcopy` (whose reflection `juliac --trim` rejects).

Every container is copied, keeping its capacity, so the copy decodes through
[`decode!`](@ref) without allocating either. The FEC decoder handles — AFF3CT
objects that live outside Julia — are shared rather than copied: they hold no
state between two decodes, but it does mean two copies of a state must not
decode concurrently from different threads.
"""
duplicate(x::Union{Number,Nothing,Symbol,Enum,String,FixedText}) = x
duplicate(x::SlotDictionary) = copy(x)
duplicate(x::Base.RefValue{T}) where {T} = Base.RefValue{T}(duplicate(x[]))
duplicate(x::Array{T}) where {T} = isbitstype(T) ? copy(x) : map(duplicate, x)
function duplicate(x::Vector{T}) where {T}
    y = isbitstype(T) ? copy(x) : map(duplicate, x)
    sizehint!(y, vector_capacity(x))
    return y
end
function duplicate(x::CircularDeque{T}) where {T}
    y = CircularDeque{T}(capacity(x))
    for v in x
        push!(y, v)
    end
    return y
end
duplicate(x::Dictionary) = Dictionary(copy(collect(keys(x))), map(duplicate, collect(values(x))))
duplicate(x::Aff3ct.ConvViterbiDecoder) = x
duplicate(x::Aff3ct.LDPCBPDecoder) = x

@generated function duplicate(x::T) where {T}
    isbitstype(T) && return :x
    # A mutable struct from another package may own memory outside Julia (the
    # AFF3CT handles do); copying its fields would alias that memory, so each
    # such type must be listed above explicitly.
    ismutabletype(T) && parentmodule(T) !== @__MODULE__() &&
        return :(throw(ArgumentError(string("no `duplicate` method for ", $T))))
    fields = [:(duplicate(getfield(x, $i))) for i = 1:fieldcount(T)]
    return Expr(:new, T, fields...)
end

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
