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

# Each container spells the ephemeris epochs, the rate of inclination and the clock
# polynomial as its ICD does; the normalised `BroadcastEphemeris` uses one name for each.
ephemeris_icd_names(::Type{<:GNSSDecoder.AbstractGPSData}) =
    (; t_0e = :t_oe, i_dot = :IDOT, t_0c = :t_oc)
ephemeris_icd_names(::Type{<:GNSSDecoder.AbstractGalileoData}) = (;)
ephemeris_icd_names(::Type{BeiDouDNAVData}) =
    (; t_0e = :t_oe, i_dot = :IDOT, t_0c = :t_oc, a_f0 = :a_0, a_f1 = :a_1, a_f2 = :a_2)
ephemeris_icd_names(::Type{<:GNSSDecoder.AbstractBeiDouCNAVData}) = (;
    t_0e = :t_oe,
    i_dot = :i_0_dot,
    t_0c = :t_oc,
    a_f0 = :a_0,
    a_f1 = :a_1,
    a_f2 = :a_2,
)

"""
    icd_fields(D, fields) -> NamedTuple

`fields`, named as `BroadcastEphemeris` names them, under the names `D`'s ICD gives the
same quantities, ready to build a `D` from.
"""
function icd_fields(D, fields::NamedTuple)
    names = ephemeris_icd_names(D)
    (; (get(names, name, name) => value for (name, value) in pairs(fields))...)
end
