# Definitions shared across Galileo signals (I/NAV on E1-B and E5b-I, F/NAV on E5a-I,
# C/NAV on E6-B): enums, the almanac record and the FEC primitive. Framing, page layout
# and parsers live in `inav.jl`, `e5a.jl`, `e6b.jl`; per-signal layers in `e1b.jl`,
# `e5b.jl`. Definitions shared by all signals are in `src/gnss.jl`.

# All Galileo open-service signals share the rate-1/2, K=7 non-systematic convolutional
# FEC: G1 = 0o171, G2 = 0o133, G2 output inverted (Galileo OS SIS ICD, Issue 2.2,
# §4.1.4). Only interleaver dimensions and codeword length differ per signal.
const GALILEO_VITERBI_POLY = [0o171, 0o133]

# Every Galileo data channel interleaves over 8 rows: I/NAV 8×30, F/NAV 8×61, C/NAV
# 8×123 (OS SIS ICD, Issue 2.2 §4.1.4.2 and Table 27; HAS SIS ICD, Issue 1.0 Table 4).
# Each signal states only its column count.
const GALILEO_INTERLEAVER_ROWS = 8

# Ephemeris user-algorithm constants (Galileo OS SIS ICD, Issue 2.2, Table 68, and
# Eq. 16 for F). Only those differing from GPS WGS-84 live here; π and the Earth
# rotation rate are in `gnss.jl`. `GALILEO_F` = -2√μ/c², hence differs too.
const GALILEO_μ = 3.986004418e14   # geocentric gravitational constant (m³/s²)
const GALILEO_F = -4.442807309e-10 # relativistic correction constant (s/√m)

"""
    SignalHealth

Galileo signal health status enumeration.

Indicates the operational status of a Galileo signal component (I/NAV word type 5,
E5a F/NAV page type 1, and the per-satellite almanacs of both).

# Values

  - `signal_ok`: Signal is operating normally (value 0)
  - `signal_out_of_service`: Signal is out of service (value 1)
  - `signal_in_extended_operations_mode`: Signal is in Extended Operations Mode (value 2)
  - `signal_component_currently_in_test`: Signal component is currently in test (value 3)

# Reference

Galileo OS SIS ICD, Issue 2.2, Table 84
"""
@enum SignalHealth begin
    signal_ok
    signal_out_of_service
    signal_in_extended_operations_mode
    signal_component_currently_in_test
end

"""
    DataValidityStatus

Galileo navigation data validity status enumeration.

Indicates whether the broadcast navigation data should be trusted for positioning.

# Values

  - `navigation_data_valid`: Navigation data is valid (value 0)
  - `working_without_guarantee`: Navigation data is working without guarantee (value 1)

# Reference

Galileo OS SIS ICD, Issue 2.2, Table 81
"""
@enum DataValidityStatus begin
    navigation_data_valid
    working_without_guarantee
end

"""
    GalileoAlmanac

Almanac data for one Galileo satellite.

The almanac provides reduced-precision orbital and clock parameters for predicting
satellite positions and selecting satellites for tracking. `Δsqrt_A` and `δ_i` are
relative to the nominal constellation (`A_nominal = 29 600 000 m`, `i_nominal = 56°`;
OS SIS ICD Issue 2.2 Table 1). Both the I/NAV decoder (word types 7-10) and the F/NAV
decoder (page types 5-6) produce this record; they populate different signal-health
fields.

# Fields

  - `SVID::Int`: Satellite identifier (1-36 nominal range; 0 = unused entry)
  - `Δsqrt_A::Float64`: Difference of √(semi-major axis) from nominal (√m)
  - `e::Float64`: Eccentricity (dimensionless)
  - `ω::Float64`: Argument of perigee (rad)
  - `δ_i::Float64`: Inclination delta from nominal (rad)
  - `Ω_0::Float64`: Longitude of ascending node at weekly epoch (rad)
  - `Ω_dot::Float64`: Rate of change of right ascension (rad/s)
  - `M_0::Float64`: Mean anomaly at reference time (rad)
  - `a_f0::Float64`: Truncated SV clock bias (seconds)
  - `a_f1::Float64`: Truncated SV clock drift (s/s)
  - `E5b_SHS::SignalHealth`: Predicted E5b signal health status (Galileo I/NAV word types 7-10; `nothing` for F/NAV-decoded almanacs)
  - `E1B_SHS::SignalHealth`: Predicted E1-B/C signal health status (ICD E1-B_SHS, hyphen
    dropped; Galileo I/NAV word types 7-10; `nothing` for F/NAV-decoded almanacs)
  - `E5a_SHS::SignalHealth`: Predicted E5a signal health status (Galileo F/NAV page types 5-6; `nothing` for I/NAV-decoded almanacs)
  - `IOD_a::Int`: Almanac IOD
  - `WN_a::Int`: Almanac reference Week Number
  - `t_0a::Int`: Almanac reference time (seconds)

# Reference

Galileo OS SIS ICD, Issue 2.2, Table 86 (the almanac parameters, shared by
I/NAV and F/NAV); bit allocations in Tables 48-51 (I/NAV word types 7-10) and
Tables 34-35 (F/NAV page types 5-6)
"""
Base.@kwdef struct GalileoAlmanac
    SVID::Union{Nothing,Int} = nothing
    Δsqrt_A::Union{Nothing,Float64} = nothing
    e::Union{Nothing,Float64} = nothing
    ω::Union{Nothing,Float64} = nothing
    δ_i::Union{Nothing,Float64} = nothing
    Ω_0::Union{Nothing,Float64} = nothing
    Ω_dot::Union{Nothing,Float64} = nothing
    M_0::Union{Nothing,Float64} = nothing
    a_f0::Union{Nothing,Float64} = nothing
    a_f1::Union{Nothing,Float64} = nothing
    E5b_SHS::Union{Nothing,SignalHealth} = nothing
    E1B_SHS::Union{Nothing,SignalHealth} = nothing
    E5a_SHS::Union{Nothing,SignalHealth} = nothing
    IOD_a::Union{Nothing,Int} = nothing
    WN_a::Union{Nothing,Int} = nothing
    t_0a::Union{Nothing,Int} = nothing
end

@inline function GalileoAlmanac(
    a::GalileoAlmanac;
    SVID = a.SVID,
    Δsqrt_A = a.Δsqrt_A,
    e = a.e,
    ω = a.ω,
    δ_i = a.δ_i,
    Ω_0 = a.Ω_0,
    Ω_dot = a.Ω_dot,
    M_0 = a.M_0,
    a_f0 = a.a_f0,
    a_f1 = a.a_f1,
    E5b_SHS = a.E5b_SHS,
    E1B_SHS = a.E1B_SHS,
    E5a_SHS = a.E5a_SHS,
    IOD_a = a.IOD_a,
    WN_a = a.WN_a,
    t_0a = a.t_0a,
)
    GalileoAlmanac(
        SVID,
        Δsqrt_A,
        e,
        ω,
        δ_i,
        Ω_0,
        Ω_dot,
        M_0,
        a_f0,
        a_f1,
        E5b_SHS,
        E1B_SHS,
        E5a_SHS,
        IOD_a,
        WN_a,
        t_0a,
    )
end

# Slots of the per-SVID almanac store (`SlotDictionary{GalileoAlmanac,
# GALILEO_ALMANAC_SLOTS}`) of I/NAV and F/NAV data. Nominal SVIDs are 1-36, but the
# SVID field is 6 bits (OS SIS ICD Tables 48-51 and 34-35), so a CRC-valid page can
# carry any key 1-63. One slot per 6-bit value keeps every decodable almanac instead
# of throwing on a key past 36.
const GALILEO_ALMANAC_SLOTS = 64

# I/NAV and F/NAV broadcast the same orbital and clock parameters, so the "all
# present?" checks are stated once, on `AbstractGalileoEphemerisData` (supertype of
# `GalileoINAVData` and `GalileoFNAVData`).
#
# `GalileoE6BData` is deliberately *not* an `AbstractGalileoEphemerisData`: C/NAV
# carries HAS corrections to another signal's ephemeris and has no orbital fields
# (dispatching on the wider supertype would raise a `FieldError`), so it defines
# `has_positioning_messages` itself. Health-status and positioning-readiness checks
# differ per signal and stay in `inav.jl` / `e5a.jl` / `e6b.jl`.
function is_ephemeris_decoded(data::AbstractGalileoEphemerisData)
    !isnothing(data.t_0e) &&
        !isnothing(data.M_0) &&
        !isnothing(data.e) &&
        !isnothing(data.sqrt_A) &&
        !isnothing(data.Ω_0) &&
        !isnothing(data.i_0) &&
        !isnothing(data.ω) &&
        !isnothing(data.i_dot) &&
        !isnothing(data.Ω_dot) &&
        !isnothing(data.Δn) &&
        !isnothing(data.C_uc) &&
        !isnothing(data.C_us) &&
        !isnothing(data.C_rc) &&
        !isnothing(data.C_rs) &&
        !isnothing(data.C_ic) &&
        !isnothing(data.C_is)
end

function is_clock_correction_decoded(data::AbstractGalileoEphemerisData)
    !isnothing(data.t_0c) &&
        !isnothing(data.a_f0) &&
        !isnothing(data.a_f1) &&
        !isnothing(data.a_f2)
end

"""
$(TYPEDEF)

Long-lived scratch for one Galileo FEC decode: the AFF3CT Viterbi handle plus the
two buffers a decode would otherwise allocate per page. One type serves I/NAV, F/NAV
and C/NAV; each cache holds one sized for its signal. Bundling the buffers with the
handle keeps [`galileo_viterbi`](@ref) a one-argument call and rules out a mismatched
decoder and buffer.

# Fields

$(TYPEDFIELDS)
"""
struct GalileoViterbiScratch
    """
    AFF3CT K=7 NSC Viterbi decoder (`K` information bits, `N` channel symbols),
    built once and reused across pages
    """
    decoder::Aff3ct.ConvViterbiDecoder
    """
    `N` deinterleaved LLRs, written by `deinterleave!` on each decode
    """
    deinterleaved::Vector{Float32}
    """
    `K` decoded information bits, written by `Aff3ct.decode!` on each decode
    """
    info_bits::Vector{Int32}
end

"""
    GalileoViterbiScratch(K, N)

Build the decoder and its two buffers for a `K`-bit, `N`-symbol Galileo
codeword: I/NAV (114, 240), F/NAV (238, 488), C/NAV (486, 984).
"""
GalileoViterbiScratch(K::Int, N::Int) = GalileoViterbiScratch(
    Aff3ct.ConvViterbiDecoder(K, N, GALILEO_VITERBI_POLY),
    Vector{Float32}(undef, N),
    Vector{Int32}(undef, K),
)

# The GGTO reference week is broadcast in 6 bits (I/NAV word type 10 bits 123-128,
# F/NAV page type 4) against a 12-bit `WN`; see `resolve_reference_week`.
const GALILEO_GGTO_WN_MODULUS = 64

# `WN_0t` and `WN_LSF` are broadcast in 8 bits (I/NAV word type 6, F/NAV page type 4)
# against a 12-bit `WN`: see `get_utc_parameters`.
const GALILEO_UTC_WN_MODULUS = 256
const GALILEO_LSF_WN_MODULUS = 256

# BGD scaling `(f_E1/f)²` for a single-frequency user: 1.79 on E5a, 1.70 on E5b, 1 on
# E1 (OS SIS ICD §5.1.5). Omitting it leaves a metre-level error per satellite.
galileo_group_delay_scaling(signal) =
    (get_center_frequency(GalileoE1B) / get_center_frequency(signal))^2

"""
    galileo_ggto_offset(state::GNSSDecoderState, target::TimeSystem)
        -> Union{Nothing,InterSystemTimeOffsetParameters}

Normalise a decoded GGTO into an `InterSystemTimeOffsetParameters`, or `nothing` when the
satellite is not broadcasting one or `target` is not GPS time.

The only Galileo GNSS Time Offset is GST to GPST (OS SIS ICD Issue 2.2 §5.1.8), so
there is no broadcast identifier to check and no other `target` can be answered. The
quadratic term of the normalised record is zero (the polynomial has two terms).

`state.data` must carry the four GGTO fields as set by [`galileo_ggto`](@ref); `state`
supplies the week number that resolves the 6-bit `WN_0G`. Callers are I/NAV word type
10 and F/NAV page type 4. This is not a method on `AbstractGalileoEphemerisData`, so
a future ephemeris signal without a GGTO does not inherit a `FieldError`.
"""
function galileo_ggto_offset(state::GNSSDecoderState, target::TimeSystem)
    data = state.data
    target === GPST() || return nothing
    # All present or all `nothing`. Each is checked anyway to narrow the types: Julia
    # 1.10 cannot split the call below over that many `Union` arguments and would
    # dispatch at runtime and allocate.
    (; A_0G, A_1G, t_0G, WN_0G) = data
    (isnothing(A_0G) || isnothing(A_1G) || isnothing(t_0G) || isnothing(WN_0G)) &&
        return nothing
    assemble_inter_system_time_offset_parameters(
        target,
        A_0G,
        A_1G,
        0.0,
        t_0G,
        WN_0G,
        data.WN,
        GALILEO_GGTO_WN_MODULUS,
    )
end

"""
    galileo_ggto(A_0G_raw, A_1G_raw, t_0G_raw, WN_0G_raw)
        -> (; A_0G, A_1G, t_0G, WN_0G)

Scale the four broadcast GNSS-Time-Offset fields, or report the ICD's
"not valid" encoding as four `nothing`s.

Galileo OS SIS ICD, Issue 2.2, 5.1.8 (above Table 76): "When the GGTO is not
available the GGTO parameters disseminated are: A0G (all ones - 16 bits), A1G
(all ones - 12 bits), t0G (all ones - 8 bits), WN0G (all ones - 6 bits). When a
user receives all four parameters set to all ones the GGTO is considered as not
valid."

The *set* is the sentinel: `A_0G` alone being `0xffff` is the legal value -2^-35 s.
Scaling an unavailable GGTO would publish `A_0G ~= -2.9e-11`, `t_0G = 918000`,
`WN_0G = 63`, indistinguishable from a real offset. E6-B treats the HAS "data not
available" sentinels the same way.

Takes the *unsigned* raw fields; `A_0G` and `A_1G` are sign-extended here, after the
test. Callers: I/NAV word type 10 and F/NAV page type 4.
"""
function galileo_ggto(
    A_0G_raw::Integer,
    A_1G_raw::Integer,
    t_0G_raw::Integer,
    WN_0G_raw::Integer,
)
    if A_0G_raw == 0xffff && A_1G_raw == 0x0fff && t_0G_raw == 0xff && WN_0G_raw == 0x3f
        return (A_0G = nothing, A_1G = nothing, t_0G = nothing, WN_0G = nothing)
    end
    return (
        A_0G = get_twos_complement_num(UInt64(A_0G_raw), 16, 1, 16) * 2.0^-35,
        A_1G = get_twos_complement_num(UInt64(A_1G_raw), 12, 1, 12) * 2.0^-51,
        t_0G = Int(t_0G_raw) * 3600,
        WN_0G = Int(WN_0G_raw),
    )
end

"""
    with_ggto(data, ggto) -> typeof(data)

Rebuild the I/NAV or F/NAV `data` with the four GGTO fields of `ggto`, the
result of [`galileo_ggto`](@ref).

A function barrier: `galileo_ggto` returns one of two `NamedTuple` types, so the call
splits into one concrete method instance per type. Destructuring into keywords would
leave four `Union{Nothing,…}` values, which Julia 1.10 compiles through the allocating
generic keyword path.
"""
@inline with_ggto(data::D, ggto::NamedTuple) where {D<:AbstractGalileoEphemerisData} =
    D(data; ggto.A_0G, ggto.A_1G, ggto.t_0G, ggto.WN_0G)

"""
    galileo_viterbi(scratch, soft_page, interleaver_columns, ::Type{T}) -> T

Recover one page's information bits from `soft_page`, the polarity-corrected
`Float32` LLRs between the leading and trailing page-sync sequences. Shared by I/NAV,
F/NAV and C/NAV. `scratch` is the caller's cached [`GalileoViterbiScratch`](@ref).

The transmit FEC chain (Galileo OS SIS ICD, Issue 2.2, §4.1.4) is undone in order:

 1. **Block deinterleave** (`deinterleave!` from `src/coding/deinterleave.jl`, into
    `scratch.deinterleaved`). The ICD interleaver is 8 rows by `interleaver_columns`
    columns; undoing it fills the *transposed* `interleaver_columns × 8` matrix by
    column and reads by row. Column counts: I/NAV 30, F/NAV 61, C/NAV 123.
 2. **Invert every second symbol** (the ICD inverts the G2 output). On soft symbols
    this is a sign flip, preserving confidence.
 3. **K=7 NSC Viterbi** via AFF3CT.jl's `ConvViterbiDecoder`, in place with
    `Aff3ct.decode!` into `scratch.info_bits`. AFF3CT's LLR sign convention matches
    ours (positive ⇒ bit 0). Only information bits are returned (the 6 tail bits are
    consumed by trellis termination).

The bits are packed MSB-first into the low bits of `T<:Unsigned` (`UInt128` for I/NAV's
114 bits, `UInt256` for F/NAV's 238, `UInt512` for C/NAV's 486).

Both buffers come from `scratch`: the allocating `deinterleave` and `Aff3ct.decode` would
cost 1.6 kB per I/NAV page part (1 s), 3.2 kB per F/NAV page (10 s) and 6.2 kB per C/NAV
page (1 s), per tracked satellite. A call allocates nothing. Callers must keep arithmetic
off the wide result: BitIntegers converts `UInt256`/`UInt512` to `Float64` through a
`BigInt`, so narrow a field to `UInt64` before scaling (see `fnav_field` in `e5a.jl`).
"""
function galileo_viterbi(
    scratch::GalileoViterbiScratch,
    soft_page::AbstractVector{Float32},
    interleaver_columns::Int,
    ::Type{T},
) where {T<:Unsigned}
    deinterleaved = deinterleave!(
        scratch.deinterleaved,
        soft_page,
        interleaver_columns,
        GALILEO_INTERLEAVER_ROWS,
    )
    @inbounds for i = 2:2:length(deinterleaved)
        deinterleaved[i] = -deinterleaved[i]
    end
    info_bits = Aff3ct.decode!(scratch.info_bits, scratch.decoder, deinterleaved)
    bits = T(0)
    @inbounds for b in info_bits
        bits = (bits << 1) | T(b)
    end
    return bits
end
