# The broadcast ionospheric coefficient sets, one type per model, read off the
# decoded data of every message that carries one. The models themselves (the delay
# a set predicts along a line of sight) are a positioning concern and live with the
# consumer.

"""
$(TYPEDEF)

The eight Klobuchar ionospheric coefficients of a GPS navigation message (LNAV,
CNAV or CNAV-2 — every civil signal broadcasts the same set), in IS-GPS-200 SI
units (seconds and seconds·semicircle⁻ⁿ). The field names mirror the GPS data
containers' (`α_0…α_3`, `β_0…β_3`). [`get_klobuchar_parameters`](@ref) returns
it.

# Fields

$(TYPEDFIELDS)
"""
struct KlobucharParameters
    α_0::Float64
    α_1::Float64
    α_2::Float64
    α_3::Float64
    β_0::Float64
    β_1::Float64
    β_2::Float64
    β_3::Float64
end

"""
$(TYPEDEF)

The eight coefficients of BeiDou's Klobuchar variant, from the legacy D1/D2
message (B1I/B3I), in BDS-SIS-ICD-B1I-3.0 §5.2.4.7 SI units (seconds and
seconds·π⁻ⁿ, per semicircle like the GPS set). A type of its own because the
ICD prescribes a different algorithm, with the delay defined at B1I
(1561.098 MHz) rather than L1. [`get_klobuchar_parameters`](@ref) returns it.

# Fields

$(TYPEDFIELDS)
"""
struct BeiDouKlobucharParameters
    α_0::Float64
    α_1::Float64
    α_2::Float64
    α_3::Float64
    β_0::Float64
    β_1::Float64
    β_2::Float64
    β_3::Float64
end

"""
$(TYPEDEF)

The Galileo Effective Ionisation Level coefficients `a_i0`/`a_i1`/`a_i2` (sfu,
sfu/deg, sfu/deg²). [`get_nequick_g_parameters`](@ref) returns it.

Galileo broadcasts them for NeQuick-G, its single-frequency ionospheric model
(Galileo OS SIS ICD §5.1.6).

Only what is broadcast: the set carries no reference epoch. The model is
evaluated at the user's epoch (its day of year and universal time), which the
caller supplies.

# Fields

$(TYPEDFIELDS)
"""
struct NeQuickGParameters
    a_i0::Float64
    a_i1::Float64
    a_i2::Float64
end

"""
$(TYPEDEF)

The nine BDGIM coefficients of a BDS-3 B-CNAV message (B1C, B2a or B2b), in
TECU. Only what is broadcast: the model is evaluated at the user's epoch (its
Modified Julian Date), which the caller supplies.

The names are the ICD's `α_1 … α_9` (BDS-SIS-ICD-B1C-1.0 Table 7-10), the data
containers' `α_bdgim_1 … α_bdgim_9`, which arrive with each scale factor and
sign applied (including `α_5`'s −2⁻³). [`get_bdgim_parameters`](@ref) returns
it.

# Fields

$(TYPEDFIELDS)
"""
struct BDGIMParameters
    α_1::Float64
    α_2::Float64
    α_3::Float64
    α_4::Float64
    α_5::Float64
    α_6::Float64
    α_7::Float64
    α_8::Float64
    α_9::Float64
end

# Every accessor answers on every decoder: a set where the message carries one,
# `nothing` where it carries none, each stated per constellation family rather
# than by a blanket fallback, for the reason `get_broadcast_time_offset` gives.

"""
    get_klobuchar_parameters(state::GNSSDecoderState)
        -> Union{Nothing,KlobucharParameters,BeiDouKlobucharParameters}

The Klobuchar coefficients of a GPS message (LNAV, CNAV or CNAV-2), or, as a
[`BeiDouKlobucharParameters`](@ref), of the BeiDou D1/D2 message; `nothing`
until all eight are decoded, and on Galileo and the BeiDou B-CNAV messages,
which broadcast none.

Read off `data`, like every accessor here. The set is not part of the
positioning readiness gate (see
[`is_decoding_completed_for_positioning`](@ref)): it is broadcast far less
often, so apply it when present rather than waiting for it.
"""
function get_klobuchar_parameters(state::GNSSDecoderState{<:AbstractGPSData})
    d = state.data
    has_klobuchar_coefficients(d) || return nothing
    KlobucharParameters(klobuchar_coefficients(d)...)
end
function get_klobuchar_parameters(state::GNSSDecoderState{<:BeiDouDNAVData})
    d = state.data
    has_klobuchar_coefficients(d) || return nothing
    BeiDouKlobucharParameters(klobuchar_coefficients(d)...)
end
get_klobuchar_parameters(
    ::GNSSDecoderState{<:Union{AbstractGalileoData,AbstractBeiDouCNAVData}},
) = nothing

# Checked field by field: `any(isnothing, …)` over a tuple of the `Union` fields
# boxes the tuple.
has_klobuchar_coefficients(d) = !(
    isnothing(d.α_0) ||
    isnothing(d.α_1) ||
    isnothing(d.α_2) ||
    isnothing(d.α_3) ||
    isnothing(d.β_0) ||
    isnothing(d.β_1) ||
    isnothing(d.β_2) ||
    isnothing(d.β_3)
)

klobuchar_coefficients(d) = (
    something(d.α_0),
    something(d.α_1),
    something(d.α_2),
    something(d.α_3),
    something(d.β_0),
    something(d.β_1),
    something(d.β_2),
    something(d.β_3),
)

"""
    get_nequick_g_parameters(state::GNSSDecoderState) -> Union{Nothing,NeQuickGParameters}

The Effective Ionisation Level coefficients of a Galileo message (I/NAV on E1-B
or E5b, F/NAV on E5a), the coefficients of NeQuick-G; `nothing` until all three
are decoded, on Galileo
E6-B, whose C/NAV broadcasts none, and on GPS and BeiDou.
"""
function get_nequick_g_parameters(state::GNSSDecoderState{<:AbstractGalileoEphemerisData})
    d = state.data
    (isnothing(d.a_i0) || isnothing(d.a_i1) || isnothing(d.a_i2)) && return nothing
    NeQuickGParameters(something(d.a_i0), something(d.a_i1), something(d.a_i2))
end
get_nequick_g_parameters(
    ::GNSSDecoderState{<:Union{GalileoE6BData,AbstractGPSData,AbstractBeiDouData}},
) = nothing

"""
    get_bdgim_parameters(state::GNSSDecoderState) -> Union{Nothing,BDGIMParameters}

The nine BDGIM coefficients of a BDS-3 B-CNAV message (B1C, B2a or B2b, all
broadcasting the same set: BDS-SIS-ICD-B1C-1.0 §7.8, B2a §7.8 and B2b §7.7);
`nothing` until all nine are decoded, on the legacy D1/D2 message, which
broadcasts Klobuchar, and on GPS and Galileo.
"""
function get_bdgim_parameters(state::GNSSDecoderState{<:AbstractBeiDouCNAVData})
    d = state.data
    (
        isnothing(d.α_bdgim_1) ||
        isnothing(d.α_bdgim_2) ||
        isnothing(d.α_bdgim_3) ||
        isnothing(d.α_bdgim_4) ||
        isnothing(d.α_bdgim_5) ||
        isnothing(d.α_bdgim_6) ||
        isnothing(d.α_bdgim_7) ||
        isnothing(d.α_bdgim_8) ||
        isnothing(d.α_bdgim_9)
    ) && return nothing
    BDGIMParameters(
        something(d.α_bdgim_1),
        something(d.α_bdgim_2),
        something(d.α_bdgim_3),
        something(d.α_bdgim_4),
        something(d.α_bdgim_5),
        something(d.α_bdgim_6),
        something(d.α_bdgim_7),
        something(d.α_bdgim_8),
        something(d.α_bdgim_9),
    )
end
get_bdgim_parameters(
    ::GNSSDecoderState{<:Union{BeiDouDNAVData,AbstractGPSData,AbstractGalileoData}},
) = nothing
