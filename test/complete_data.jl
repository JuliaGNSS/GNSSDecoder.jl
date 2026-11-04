"""
    complete_data(D; overrides...)

A `D` with every scalar field set (numbers to 1, flags to `false`, enums to
their first value), which passes `D`'s readiness gate, with `overrides` on top.
The group delay reads only the fields the overrides name, so the rest need only
be present.
"""
function complete_data(D; overrides...)
    fields = Dict{Symbol,Any}()
    for (name, T) in zip(fieldnames(D), fieldtypes(D))
        inner = T isa Union ? Base.typesplit(T, Nothing) : T
        if inner <: Bool
            fields[name] = false
        elseif inner <: Real
            fields[name] = inner(1)
        elseif inner <: Enum
            fields[name] = first(instances(inner))
        end
    end
    D(; merge(fields, Dict(overrides))...)
end

complete_state(state, D; overrides...) =
    GNSSDecoderState(state; data = complete_data(D; overrides...))
