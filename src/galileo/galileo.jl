# Definitions shared across Galileo signals (I/NAV on E1-B and E5b-I, F/NAV on
# E5a-I, C/NAV on E6-B): the signal-health / data-validity enums, the
# per-satellite almanac record, and the common
# forward-error-correction primitive. Each message's framing, page layout, and
# parser live in its own file (`inav.jl`, `e5a.jl`, `e6b.jl`), with the
# per-signal layers on top (`e1b.jl`, `e5b.jl`). This mirrors how `src/gnss.jl`
# holds the definitions shared across *all* signals.

# All Galileo open-service signals share the same rate-1/2, constraint-length-7
# (K=7) non-systematic convolutional (NSC) FEC: generator polynomials G1 = 0o171,
# G2 = 0o133, with the G2 output inverted (Galileo OS SIS ICD, Issue 2.2, §4.1.4).
# Only the block-interleaver dimensions and codeword length differ per signal.
const GALILEO_VITERBI_POLY = [0o171, 0o133]

# Every Galileo data channel interleaves over 8 rows — I/NAV 8×30, E5a F/NAV
# 8×61, E6-B C/NAV 8×123 (OS SIS ICD, Issue 2.2 §4.1.4.2 and Table 5; HAS SIS
# ICD, Issue 1.0 Table 4). Only the column count differs, so it is the one number
# each signal states, and `galileo_viterbi` takes that alone.
const GALILEO_INTERLEAVER_ROWS = 8

# GTRF constants specific to Galileo (Galileo OS SIS ICD, Issue 2.2, Table 68).
# Only the two that differ from the GPS WGS-84 values live here; π, the speed of
# light, and the Earth rotation rate are shared package-wide (`GNSS_PI`,
# `SPEED_OF_LIGHT`, `EARTH_ROTATION_RATE` in `gnss.jl`). `GALILEO_F` follows from
# μ via F = -2√μ/c², which is why it too differs from GPS.
const GALILEO_μ = 3.986004418e14   # geocentric gravitational constant (m³/s²)
const GALILEO_F = -4.442807309e-10 # relativistic correction constant (s/√m)

"""
    SignalHealth

Galileo signal health status enumeration.

Indicates the operational status of a Galileo signal component (I/NAV word type 5,
E5a F/NAV word type 1, and the per-satellite almanacs of both).

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
satellite positions and selecting satellites for tracking. Differences (`Δsqrt_A`,
`δi`) are relative to nominal Galileo constellation values (`A_nominal = 29600.318 km`,
`i_nominal = 56°`). The same record is produced by both the I/NAV decoder (word
types 7-10) and the F/NAV decoder (word types 5-6); they differ only in which
signal-health facet they populate (see below).

# Fields

  - `SVID::Int`: Satellite identifier (1-36 nominal range; 0 = unused entry)
  - `Δsqrt_A::Float64`: Difference of √(semi-major axis) from nominal (√m)
  - `e::Float64`: Eccentricity (dimensionless)
  - `ω::Float64`: Argument of perigee (semi-circles)
  - `δi::Float64`: Inclination delta from nominal (semi-circles)
  - `Ω_0::Float64`: Longitude of ascending node at weekly epoch (semi-circles)
  - `Ω_dot::Float64`: Rate of change of right ascension (semi-circles/s)
  - `M_0::Float64`: Mean anomaly at reference time (semi-circles)
  - `a_f0::Float64`: Truncated SV clock bias (seconds)
  - `a_f1::Float64`: Truncated SV clock drift (s/s)
  - `signal_health_e5b::SignalHealth`: Predicted E5b signal health status (Galileo I/NAV word types 7-10)
  - `signal_health_e1b::SignalHealth`: Predicted E1-B/C signal health status (Galileo I/NAV word types 7-10)
  - `signal_health_e5a::SignalHealth`: Predicted E5a signal health status (Galileo F/NAV word types 5-6; `nothing` for I/NAV-decoded almanacs)
  - `IOD_a::Int`: Almanac IOD
  - `WN_a::Int`: Almanac reference Week Number
  - `t_0a::Int`: Almanac reference time (seconds)

# Reference

Galileo OS SIS ICD, Issue 2.2, Table 86 (I/NAV) and Tables 75-76 (F/NAV)
"""
Base.@kwdef struct GalileoAlmanac
    SVID::Union{Nothing,Int} = nothing
    Δsqrt_A::Union{Nothing,Float64} = nothing
    e::Union{Nothing,Float64} = nothing
    ω::Union{Nothing,Float64} = nothing
    δi::Union{Nothing,Float64} = nothing
    Ω_0::Union{Nothing,Float64} = nothing
    Ω_dot::Union{Nothing,Float64} = nothing
    M_0::Union{Nothing,Float64} = nothing
    a_f0::Union{Nothing,Float64} = nothing
    a_f1::Union{Nothing,Float64} = nothing
    signal_health_e5b::Union{Nothing,SignalHealth} = nothing
    signal_health_e1b::Union{Nothing,SignalHealth} = nothing
    signal_health_e5a::Union{Nothing,SignalHealth} = nothing
    IOD_a::Union{Nothing,Int} = nothing
    WN_a::Union{Nothing,Int} = nothing
    t_0a::Union{Nothing,Int} = nothing
end

function GalileoAlmanac(
    a::GalileoAlmanac;
    SVID = a.SVID,
    Δsqrt_A = a.Δsqrt_A,
    e = a.e,
    ω = a.ω,
    δi = a.δi,
    Ω_0 = a.Ω_0,
    Ω_dot = a.Ω_dot,
    M_0 = a.M_0,
    a_f0 = a.a_f0,
    a_f1 = a.a_f1,
    signal_health_e5b = a.signal_health_e5b,
    signal_health_e1b = a.signal_health_e1b,
    signal_health_e5a = a.signal_health_e5a,
    IOD_a = a.IOD_a,
    WN_a = a.WN_a,
    t_0a = a.t_0a,
)
    GalileoAlmanac(
        SVID,
        Δsqrt_A,
        e,
        ω,
        δi,
        Ω_0,
        Ω_dot,
        M_0,
        a_f0,
        a_f1,
        signal_health_e5b,
        signal_health_e1b,
        signal_health_e5a,
        IOD_a,
        WN_a,
        t_0a,
    )
end

# Ephemeris/clock completeness are per-message facts: I/NAV (E1-B, E5b-I) and
# F/NAV (E5a-I) broadcast the same orbital and clock parameters, so the "all
# present?" checks are identical and are stated once — on
# `AbstractGalileoEphemerisData`, the supertype of exactly the Galileo data types
# that have these fields (`GalileoINAVData`, `GalileoE5aData`).
#
# `GalileoE6BData` is an `AbstractGalileoData` but deliberately *not* an
# `AbstractGalileoEphemerisData`: E6-B's C/NAV broadcasts HAS corrections *to*
# another signal's ephemeris and has no orbital fields, so it defines
# `is_decoding_completed_for_positioning` itself. Dispatching these on the wider
# supertype would give it methods that raise a `FieldError` on its own subtype.
# The health-status and positioning-readiness checks genuinely differ per signal
# and stay in `inav.jl` / `e5a.jl` / `e6b.jl`.
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
    galileo_ggto(A_0G_raw, A_1G_raw, t_0G_raw, WN_0G_raw)
        -> (; A_0G, A_1G, t_0G, WN_0G)

Scale the four broadcast GNSS-Time-Offset fields, or report the ICD's
"not valid" encoding as four `nothing`s.

Galileo OS SIS ICD, Issue 2.2, 5.1.8 (above Table 76): "When the GGTO is not
available the GGTO parameters disseminated are: A0G (all ones - 16 bits), A1G
(all ones - 12 bits), t0G (all ones - 8 bits), WN0G (all ones - 6 bits). When a
user receives all four parameters set to all ones the GGTO is considered as not
valid."

All four must be all-ones together — the ICD makes the *set* the sentinel, not
any one field, and `A_0G` alone being `0xffff` is the perfectly legal value
-2^-35 s. Scaling an unavailable GGTO regardless would publish
`A_0G ~= -2.9e-11`, `t_0G = 918000`, `WN_0G = 63`: numbers a consumer has no way
to tell from a real offset. Reporting `nothing` is the same treatment E6-B gives
the HAS "data not available" sentinels.

Takes the *unsigned* raw fields; the sign extension of `A_0G` and `A_1G` happens
here, after the test.

The two I/NAV/F/NAV call sites are the only ones: GGTO rides in I/NAV word
type 10 and F/NAV page type 4.
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
    galileo_viterbi(decoder, soft_page, interleaver_columns, ::Type{T}) -> T

Recover one Galileo page's information bits from `soft_page` — the
polarity-corrected `Float32` LLR soft symbols between the leading and trailing
page-sync sequences. Shared by I/NAV (E1-B, E5b-I), E5a F/NAV, and E6-B C/NAV,
which all use the same FEC and differ only in the block-interleaver shape and
codeword length. `decoder` is the caller's long-lived
`Aff3ct.ConvViterbiDecoder`, reused across pages.

The transmit FEC chain (Galileo OS SIS ICD, Issue 2.2, §4.1.4) is undone in order:

 1. **Block deinterleave** of the LLRs (`deinterleave` from
    `src/deinterleave.jl`). The ICD's interleaver is 8 rows by
    `interleaver_columns` columns; undoing it means filling the *transposed*
    matrix — `interleaver_columns × 8`, written by column and read by row — so
    only the column count varies per signal (I/NAV 30, E5a 61, E6-B C/NAV 123)
    and `GALILEO_INTERLEAVER_ROWS` supplies the other dimension.
 2. **Invert every second symbol** — the spec inverts the G2 output of the
    rate-1/2 encoder. On soft symbols an inversion is a sign flip (negation), so
    confidence magnitudes are preserved.
 3. **K=7 NSC Viterbi** via AFF3CT.jl's `ConvViterbiDecoder`. AFF3CT's LLR sign
    convention matches ours (positive ⇒ bit 0), so the LLRs feed in directly. The
    decoder returns only the information bits (the 6 tail bits are consumed by
    trellis termination).

The decoded bits are packed MSB-first into the low bits of `T<:Unsigned` (I/NAV
uses `UInt128` for its 114 bits, E5a `UInt256` for its 238, E6-B `UInt512` for
its 486).
"""
function galileo_viterbi(
    decoder::Aff3ct.ConvViterbiDecoder,
    soft_page::AbstractVector{Float32},
    interleaver_columns::Int,
    ::Type{T},
) where {T<:Unsigned}
    deinterleaved = deinterleave(soft_page, interleaver_columns, GALILEO_INTERLEAVER_ROWS)
    @inbounds for i = 2:2:length(deinterleaved)
        deinterleaved[i] = -deinterleaved[i]
    end
    info_bits = Aff3ct.decode(decoder, deinterleaved)
    bits = T(0)
    @inbounds for b in info_bits
        bits = (bits << 1) | T(b)
    end
    return bits
end
