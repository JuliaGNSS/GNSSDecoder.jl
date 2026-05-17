# UInt288 buffer for Galileo E1B
# which holds at least a complete Galileo E1B page
# plus 10 extra syncronization bits
BitIntegers.@define_integers 288

"""
    GalileoE1BConstants

GTRF constants and I/NAV message structure parameters for Galileo E1B signal decoding.

The physical constants are defined in the Galileo OS SIS ICD (Open Service Signal-In-Space
Interface Control Document) and are used for computing satellite positions and clock
corrections from broadcast ephemeris data.

# Fields
- `syncro_sequence_length::Int`: Length of synchronization sequence in bits (250 bits per page)
- `preamble::UInt16`: Page synchronization pattern (0101100000 binary)
- `preamble_length::Int`: Length of preamble in bits (10)
- `PI::Float64`: Mathematical constant π = 3.1415926535898 (Galileo OS SIS ICD Table 68)
- `Ω_dot_e::Float64`: Mean angular velocity of the Earth = 7.2921151467×10⁻⁵ rad/s
- `c::Float64`: Speed of light = 2.99792458×10⁸ m/s
- `μ::Float64`: Geocentric gravitational constant = 3.986004418×10¹⁴ m³/s²
- `F::Float64`: Relativistic correction constant = -4.442807309×10⁻¹⁰ s/√m

# Reference
Galileo OS SIS ICD, Issue 2.2, Table 68
"""
Base.@kwdef struct GalileoE1BConstants <: AbstractGNSSConstants
    syncro_sequence_length::Int = 250
    preamble::UInt16 = 0b0101100000
    preamble_length::Int = 10
    PI::Float64 = 3.1415926535898
    Ω_dot_e::Float64 = 7.2921151467e-5
    c::Float64 = 2.99792458e8
    μ::Float64 = 3.986004418e14
    F::Float64 = -4.442807309e-10
end

"""
    SignalHealth

Galileo signal health status enumeration.

Indicates the operational status of a Galileo signal component as broadcast in word type 5.

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

Almanac data for one Galileo satellite, decoded from word types 7-10.

The almanac provides reduced-precision orbital and clock parameters for predicting
satellite positions and selecting satellites for tracking. Differences (`Δsqrt_A`,
`δi`) are relative to nominal Galileo constellation values (`A_nominal = 29600.318 km`,
`i_nominal = 56°`).

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
- `signal_health_e5b::SignalHealth`: Predicted E5b signal health status
- `signal_health_e1b::SignalHealth`: Predicted E1-B/C signal health status
- `IOD_a::Int`: Almanac IOD
- `WN_a::Int`: Almanac reference Week Number
- `t_0a::Int`: Almanac reference time (seconds)

# Reference
Galileo OS SIS ICD, Issue 2.2, Table 86
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
        IOD_a,
        WN_a,
        t_0a,
    )
end

"""
    GalileoReducedCED

Reduced Clock and Ephemeris Data, decoded from word type 16.

A compact ephemeris/clock set transmitted within a single I/NAV word. Provides
reduced-accuracy parameters that allow a receiver to compute an initial position
fix before the full ephemeris (words 1-4) has been collected. Reduced CED
parameters must NOT be combined with full-precision parameters from words 1-4.

# Fields
- `ΔA_red::Float64`: Difference of semi-major axis from nominal (meters)
- `e_x_red::Float64`: Eccentricity vector x-component, e·cos(ω) (dimensionless)
- `e_y_red::Float64`: Eccentricity vector y-component, e·sin(ω) (dimensionless)
- `Δi_0_red::Float64`: Inclination delta from nominal (semi-circles)
- `Ω_0_red::Float64`: Longitude of ascending node at weekly epoch (semi-circles)
- `λ_0_red::Float64`: Mean argument of latitude, M0+ω (semi-circles)
- `a_f0_red::Float64`: SV clock bias correction coefficient (seconds)
- `a_f1_red::Float64`: SV clock drift correction coefficient (s/s)

# Reference
Galileo OS SIS ICD, Issue 2.2, Table 87
"""
Base.@kwdef struct GalileoReducedCED
    ΔA_red::Union{Nothing,Float64} = nothing
    e_x_red::Union{Nothing,Float64} = nothing
    e_y_red::Union{Nothing,Float64} = nothing
    Δi_0_red::Union{Nothing,Float64} = nothing
    Ω_0_red::Union{Nothing,Float64} = nothing
    λ_0_red::Union{Nothing,Float64} = nothing
    a_f0_red::Union{Nothing,Float64} = nothing
    a_f1_red::Union{Nothing,Float64} = nothing
end

# Page is splitted in even and odd parts
# Cache even part and decode after odd part
# Page contains 120 bits
#
# Almanac chain partials carry across word types 7-10 within one subframe:
#   WT7  fills SV-position-1 first half
#   WT8  fills SV-position-1 second half (flush to almanacs[SVID])
#         and SV-position-2 first half
#   WT9  fills SV-position-2 second half (flush)
#         and SV-position-3 first half
#   WT10 fills SV-position-3 second half (flush)
struct GalileoE1BCache <: AbstractGNSSCache
    even_page_part_bits::UInt128
    almanac_chain_pos1::GalileoAlmanac
    almanac_chain_pos2::GalileoAlmanac
end

GalileoE1BCache() =
    GalileoE1BCache(UInt128(0), GalileoAlmanac(), GalileoAlmanac())

function GalileoE1BCache(
    cache::GalileoE1BCache;
    even_page_part_bits = cache.even_page_part_bits,
    almanac_chain_pos1 = cache.almanac_chain_pos1,
    almanac_chain_pos2 = cache.almanac_chain_pos2,
)
    GalileoE1BCache(
        even_page_part_bits,
        almanac_chain_pos1,
        almanac_chain_pos2,
    )
end

"""
    GalileoE1BData

Decoded Galileo E1B I/NAV navigation message data.

Contains ephemeris, clock correction, signal health, group delay, ionospheric
correction, GST-UTC and GST-GPS conversion, almanac, and Reduced CED parameters
decoded from the Galileo I/NAV message. All parameters conform to the Galileo
OS SIS ICD, Issue 2.2.

# Galileo System Time (GST) Fields
- `WN::Int64`: Week Number (0-4095)
- `TOW::Int64`: Time of Week at message transmission (seconds, 0-604799)

# Satellite Identification (Word Type 4)
- `SVID::Int`: Satellite Identifier (1-36 nominal range)

# Ephemeris Parameters (Word Types 1-3)
- `t_0e::Float64`: Ephemeris reference time (seconds)
- `M_0::Float64`: Mean anomaly at reference time (semi-circles)
- `e::Float64`: Eccentricity (dimensionless)
- `sqrt_A::Float64`: Square root of semi-major axis (√m)
- `Ω_0::Float64`: Longitude of ascending node at weekly epoch (semi-circles)
- `i_0::Float64`: Inclination angle at reference time (semi-circles)
- `ω::Float64`: Argument of perigee (semi-circles)
- `i_dot::Float64`: Rate of change of inclination angle (semi-circles/s)
- `Ω_dot::Float64`: Rate of change of right ascension (semi-circles/s)
- `Δn::Float64`: Mean motion difference from computed value (semi-circles/s)
- `C_uc::Float64`: Cosine harmonic correction to argument of latitude (rad)
- `C_us::Float64`: Sine harmonic correction to argument of latitude (rad)
- `C_rc::Float64`: Cosine harmonic correction to orbit radius (meters)
- `C_rs::Float64`: Sine harmonic correction to orbit radius (meters)
- `C_ic::Float64`: Cosine harmonic correction to inclination (rad)
- `C_is::Float64`: Sine harmonic correction to inclination (rad)

# Signal-In-Space Accuracy (Word Type 3)
- `SISA_e1_e5b::Int`: SISA index for dual frequency E1-E5b (Table 91/92; 255 = NAPA)

# Clock Correction Parameters (Word Type 4)
- `t_0c::Float64`: Clock correction reference time (seconds)
- `a_f0::Float64`: SV clock bias correction coefficient (seconds)
- `a_f1::Float64`: SV clock drift correction coefficient (s/s)
- `a_f2::Float64`: SV clock drift rate correction coefficient (s/s²)

# Issue of Data (Word Types 1-4)
- `IOD_nav1::UInt`: Issue of Data from word type 1 (10-bit)
- `IOD_nav2::UInt`: Issue of Data from word type 2 (10-bit)
- `IOD_nav3::UInt`: Issue of Data from word type 3 (10-bit)
- `IOD_nav4::UInt`: Issue of Data from word type 4 (10-bit)
- `num_pages_after_last_TOW::Int`: Pages decoded since last TOW update
- `num_bits_after_valid_syncro_sequence_after_last_TOW::Int`: Bits since last TOW sync

# Signal Health and Data Validity (Word Type 5)
- `signal_health_e1b::SignalHealth`: E1-B/C signal health status (0=OK, 1=out of service, 2=Extended Operations Mode, 3=in test)
- `signal_health_e5b::SignalHealth`: E5b signal health status
- `data_validity_status_e1b::DataValidityStatus`: E1-B data validity (0=valid, 1=working without guarantee)
- `data_validity_status_e5b::DataValidityStatus`: E5b data validity

# Broadcast Group Delay (Word Type 5)
- `broadcast_group_delay_e1_e5a::Float64`: E1-E5a group delay correction (seconds)
- `broadcast_group_delay_e1_e5b::Float64`: E1-E5b group delay correction (seconds)

# Ionospheric Correction (Word Type 5)
- `a_i0::Float64`: Effective Ionisation Level 1st-order coefficient (sfu)
- `a_i1::Float64`: Effective Ionisation Level 2nd-order coefficient (sfu/degree)
- `a_i2::Float64`: Effective Ionisation Level 3rd-order coefficient (sfu/degree²)
- `iono_storm_flag_region1..5::Bool`: Ionospheric Disturbance (storm) flags for regions 1-5

# GST-UTC Conversion (Word Type 6)
- `A_0_utc::Float64`: Constant term of polynomial (s)
- `A_1_utc::Float64`: 1st-order term of polynomial (s/s)
- `Δt_LS::Int`: Leap Second count before leap second adjustment (s)
- `t_0t::Int`: UTC data reference Time of Week (s)
- `WN_0t::Int`: UTC data reference Week Number (8-bit, modulo 256)
- `WN_LSF::Int`: Week Number of leap second adjustment (8-bit, modulo 256)
- `DN::Int`: Day Number at end of which leap second becomes effective (1=Sunday … 7=Saturday)
- `Δt_LSF::Int`: Leap Second count after leap second adjustment (s)

# GST-GPS Conversion / GGTO (Word Type 10)
- `A_0G::Float64`: Constant term of GST-GPS offset polynomial (s)
- `A_1G::Float64`: Rate of change of GST-GPS offset (s/s)
- `t_0G::Int`: GGTO reference time (s)
- `WN_0G::Int`: GGTO reference Week Number (6-bit)

# Almanac (Word Types 7-10)
- `almanacs::Dictionary{Int,GalileoAlmanac}`: Decoded almanacs keyed by SVID.
  Entries are inserted as Galileo broadcasts the almanac chain across word
  types 7→10. SVIDs not yet seen are absent from the dictionary. In-flight
  chain partials live in the decoder cache and are flushed here only once a
  full almanac for an SVID has been assembled with a consistent IODa.

# Reduced Clock and Ephemeris Data (Word Type 16)
- `reduced_ced::GalileoReducedCED`: Reduced CED for fast initial fix

# Reference
Galileo OS SIS ICD, Issue 2.2, Tables 42-55, 67-87
"""
Base.@kwdef struct GalileoE1BData <: AbstractGNSSData
    WN::Union{Nothing,Int64} = nothing
    TOW::Union{Nothing,Int64} = nothing

    SVID::Union{Nothing,Int} = nothing

    t_0e::Union{Nothing,Float64} = nothing
    M_0::Union{Nothing,Float64} = nothing
    e::Union{Nothing,Float64} = nothing
    sqrt_A::Union{Nothing,Float64} = nothing
    Ω_0::Union{Nothing,Float64} = nothing
    i_0::Union{Nothing,Float64} = nothing
    ω::Union{Nothing,Float64} = nothing
    i_dot::Union{Nothing,Float64} = nothing
    Ω_dot::Union{Nothing,Float64} = nothing
    Δn::Union{Nothing,Float64} = nothing
    C_uc::Union{Nothing,Float64} = nothing
    C_us::Union{Nothing,Float64} = nothing
    C_rc::Union{Nothing,Float64} = nothing
    C_rs::Union{Nothing,Float64} = nothing
    C_ic::Union{Nothing,Float64} = nothing
    C_is::Union{Nothing,Float64} = nothing

    SISA_e1_e5b::Union{Nothing,Int} = nothing

    t_0c::Union{Nothing,Float64} = nothing
    a_f0::Union{Nothing,Float64} = nothing
    a_f1::Union{Nothing,Float64} = nothing
    a_f2::Union{Nothing,Float64} = nothing

    IOD_nav1::Union{Nothing,UInt} = nothing
    IOD_nav2::Union{Nothing,UInt} = nothing
    IOD_nav3::Union{Nothing,UInt} = nothing
    IOD_nav4::Union{Nothing,UInt} = nothing
    num_pages_after_last_TOW::Int = 0
    num_bits_after_valid_syncro_sequence_after_last_TOW::Union{Nothing,Int} = nothing

    signal_health_e1b::Union{Nothing,SignalHealth} = nothing
    signal_health_e5b::Union{Nothing,SignalHealth} = nothing
    data_validity_status_e1b::Union{Nothing,DataValidityStatus} = nothing
    data_validity_status_e5b::Union{Nothing,DataValidityStatus} = nothing

    broadcast_group_delay_e1_e5a::Union{Nothing,Float64} = nothing
    broadcast_group_delay_e1_e5b::Union{Nothing,Float64} = nothing

    a_i0::Union{Nothing,Float64} = nothing
    a_i1::Union{Nothing,Float64} = nothing
    a_i2::Union{Nothing,Float64} = nothing
    iono_storm_flag_region1::Union{Nothing,Bool} = nothing
    iono_storm_flag_region2::Union{Nothing,Bool} = nothing
    iono_storm_flag_region3::Union{Nothing,Bool} = nothing
    iono_storm_flag_region4::Union{Nothing,Bool} = nothing
    iono_storm_flag_region5::Union{Nothing,Bool} = nothing

    A_0_utc::Union{Nothing,Float64} = nothing
    A_1_utc::Union{Nothing,Float64} = nothing
    Δt_LS::Union{Nothing,Int} = nothing
    t_0t::Union{Nothing,Int} = nothing
    WN_0t::Union{Nothing,Int} = nothing
    WN_LSF::Union{Nothing,Int} = nothing
    DN::Union{Nothing,Int} = nothing
    Δt_LSF::Union{Nothing,Int} = nothing

    A_0G::Union{Nothing,Float64} = nothing
    A_1G::Union{Nothing,Float64} = nothing
    t_0G::Union{Nothing,Int} = nothing
    WN_0G::Union{Nothing,Int} = nothing

    almanacs::Union{Nothing,Dictionary{Int,GalileoAlmanac}} = nothing

    reduced_ced::GalileoReducedCED = GalileoReducedCED()
end

function GalileoE1BData(
    data::GalileoE1BData;
    WN = data.WN,
    TOW = data.TOW,
    SVID = data.SVID,
    t_0e = data.t_0e,
    M_0 = data.M_0,
    e = data.e,
    sqrt_A = data.sqrt_A,
    Ω_0 = data.Ω_0,
    i_0 = data.i_0,
    ω = data.ω,
    i_dot = data.i_dot,
    Ω_dot = data.Ω_dot,
    Δn = data.Δn,
    C_uc = data.C_uc,
    C_us = data.C_us,
    C_rc = data.C_rc,
    C_rs = data.C_rs,
    C_ic = data.C_ic,
    C_is = data.C_is,
    SISA_e1_e5b = data.SISA_e1_e5b,
    t_0c = data.t_0c,
    a_f0 = data.a_f0,
    a_f1 = data.a_f1,
    a_f2 = data.a_f2,
    IOD_nav1 = data.IOD_nav1,
    IOD_nav2 = data.IOD_nav2,
    IOD_nav3 = data.IOD_nav3,
    IOD_nav4 = data.IOD_nav4,
    num_pages_after_last_TOW = data.num_pages_after_last_TOW,
    num_bits_after_valid_syncro_sequence_after_last_TOW = data.num_bits_after_valid_syncro_sequence_after_last_TOW,
    signal_health_e1b = data.signal_health_e1b,
    signal_health_e5b = data.signal_health_e5b,
    data_validity_status_e1b = data.data_validity_status_e1b,
    data_validity_status_e5b = data.data_validity_status_e5b,
    broadcast_group_delay_e1_e5a = data.broadcast_group_delay_e1_e5a,
    broadcast_group_delay_e1_e5b = data.broadcast_group_delay_e1_e5b,
    a_i0 = data.a_i0,
    a_i1 = data.a_i1,
    a_i2 = data.a_i2,
    iono_storm_flag_region1 = data.iono_storm_flag_region1,
    iono_storm_flag_region2 = data.iono_storm_flag_region2,
    iono_storm_flag_region3 = data.iono_storm_flag_region3,
    iono_storm_flag_region4 = data.iono_storm_flag_region4,
    iono_storm_flag_region5 = data.iono_storm_flag_region5,
    A_0_utc = data.A_0_utc,
    A_1_utc = data.A_1_utc,
    Δt_LS = data.Δt_LS,
    t_0t = data.t_0t,
    WN_0t = data.WN_0t,
    WN_LSF = data.WN_LSF,
    DN = data.DN,
    Δt_LSF = data.Δt_LSF,
    A_0G = data.A_0G,
    A_1G = data.A_1G,
    t_0G = data.t_0G,
    WN_0G = data.WN_0G,
    almanacs = data.almanacs,
    reduced_ced = data.reduced_ced,
)
    GalileoE1BData(
        WN,
        TOW,
        SVID,
        t_0e,
        M_0,
        e,
        sqrt_A,
        Ω_0,
        i_0,
        ω,
        i_dot,
        Ω_dot,
        Δn,
        C_uc,
        C_us,
        C_rc,
        C_rs,
        C_ic,
        C_is,
        SISA_e1_e5b,
        t_0c,
        a_f0,
        a_f1,
        a_f2,
        IOD_nav1,
        IOD_nav2,
        IOD_nav3,
        IOD_nav4,
        num_pages_after_last_TOW,
        num_bits_after_valid_syncro_sequence_after_last_TOW,
        signal_health_e1b,
        signal_health_e5b,
        data_validity_status_e1b,
        data_validity_status_e5b,
        broadcast_group_delay_e1_e5a,
        broadcast_group_delay_e1_e5b,
        a_i0,
        a_i1,
        a_i2,
        iono_storm_flag_region1,
        iono_storm_flag_region2,
        iono_storm_flag_region3,
        iono_storm_flag_region4,
        iono_storm_flag_region5,
        A_0_utc,
        A_1_utc,
        Δt_LS,
        t_0t,
        WN_0t,
        WN_LSF,
        DN,
        Δt_LSF,
        A_0G,
        A_1G,
        t_0G,
        WN_0G,
        almanacs,
        reduced_ced,
    )
end

# The default `==` for structs falls back to `===`, which is reference equality
# and so fails for the mutable `almanacs::Dictionary{...}` field even when the
# contents match. Compare field-by-field.
function Base.:(==)(a::GalileoE1BData, b::GalileoE1BData)
    for f in fieldnames(GalileoE1BData)
        getfield(a, f) == getfield(b, f) || return false
    end
    return true
end

function is_ephemeris_decoded(data::GalileoE1BData)
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

function is_clock_correction_decoded(data::GalileoE1BData)
    !isnothing(data.t_0c) &&
        !isnothing(data.a_f0) &&
        !isnothing(data.a_f1) &&
        !isnothing(data.a_f2)
end

function is_health_status_decoded(data::GalileoE1BData)
    !isnothing(data.signal_health_e1b) &&
        !isnothing(data.signal_health_e5b) &&
        !isnothing(data.data_validity_status_e1b) &&
        !isnothing(data.data_validity_status_e5b)
end

function is_decoding_completed_for_positioning(data::GalileoE1BData)
    !isnothing(data.TOW) &&
        !isnothing(data.WN) &&
        !isnothing(data.broadcast_group_delay_e1_e5a) &&
        !isnothing(data.broadcast_group_delay_e1_e5b) &&
        is_ephemeris_decoded(data) &&
        is_clock_correction_decoded(data) &&
        is_health_status_decoded(data)
end

"""
$(TYPEDSIGNATURES)

Create a decoder state for Galileo E1B I/NAV navigation messages.

Initializes a [`GNSSDecoderState`](@ref) configured for decoding Galileo E1B
(Open Service) navigation messages. The decoder extracts ephemeris, clock
correction, ionospheric parameters, and health data from the 250 bps I/NAV
data stream using Viterbi decoding.

# Arguments
- `prn::Int`: Pseudo-Random Noise code identifier (1-36 for Galileo satellites)

# Returns
- `GNSSDecoderState{GalileoE1BData}`: Initialized decoder state for Galileo E1B

# Example
```julia
state = GalileoE1BDecoderState(1)  # Create decoder for PRN 1
state = decode(state, bits, num_bits)
if is_sat_healthy(state)
    # Use state.data for positioning
end
```

# See Also
- [`GNSSDecoderState`](@ref): The underlying state structure
- [`decode`](@ref): Decode bits using this state
- [`reset_decoder_state`](@ref): Reset after signal loss
- [`is_sat_healthy`](@ref): Check satellite health status
"""
function GalileoE1BDecoderState(prn)
    GNSSDecoderState(
        prn,
        UInt288(0),
        UInt288(0),
        GalileoE1BData(),
        GalileoE1BData(),
        GalileoE1BConstants(),
        GalileoE1BCache(),
        0,
        nothing,
        false,
    )
end

function GNSSDecoderState(system::GalileoE1B, prn)
    GNSSDecoderState(
        prn,
        UInt288(0),
        UInt288(0),
        GalileoE1BData(),
        GalileoE1BData(),
        GalileoE1BConstants(),
        GalileoE1BCache(),
        0,
        nothing,
        false,
    )
end

"""
$(TYPEDSIGNATURES)

Reset the Galileo E1B decoder state after a signal loss or reacquisition.

Clears the bit buffers and time-of-week (TOW) field while preserving other
decoded ephemeris and clock data in `raw_data`. This allows faster recovery
after brief signal outages without requiring a full re-decode of all pages.

!!! note
    The week number (`WN`) field is intentionally not reset as it is not
    broadcast as frequently as TOW. This may cause brief errors if a week
    rollover occurs during a signal outage.

# Arguments
- `state::GNSSDecoderState{<:GalileoE1BData}`: Current Galileo E1B decoder state

# Returns
- `GNSSDecoderState{<:GalileoE1BData}`: Reset decoder state with cleared buffers

# Example
```julia
# After detecting signal loss
state = reset_decoder_state(state)
# Continue decoding with preserved ephemeris
state = decode(state, new_bits, num_bits)
```

# See Also
- [`GalileoE1BDecoderState`](@ref): Create a fresh decoder state
- [`decode`](@ref): Continue decoding after reset
"""
function reset_decoder_state(state::GNSSDecoderState{<:GalileoE1BData})
    # Reset bit buffers and TOW data field, while keeping the
    # remaining parameters in raw_data. This allows a GNSSReceiver
    # to use a satellite after a reacquisition without waiting for
    # the decoding of all data fields.
    # Note: WN is currently not reset as it is broadcast not as
    # frequently as the TOW and thus may increase the time until
    # the decoder is available again after an outage. This will
    # lead to erroneous decoder information for a few seconds after
    # reacquisition when a new week started during a signal outage.
    GNSSDecoderState(
        state;
        raw_buffer = UInt288(0),
        buffer = UInt288(0),
        raw_data = GalileoE1BData(
            state.raw_data;
            TOW = nothing,
            num_bits_after_valid_syncro_sequence_after_last_TOW = nothing,
        ),
        data = GalileoE1BData(),
        num_bits_buffered = 0,
        num_bits_after_valid_syncro_sequence = nothing
    )
end

function decode_syncro_sequence(state::GNSSDecoderState{<:GalileoE1BData})
    encoded_bits = bitstring(state.buffer >> state.constants.preamble_length)[sizeof(
        state.buffer,
    )*8-state.constants.syncro_sequence_length+state.constants.preamble_length+1:end]
    deinterleaved_encoded_bits = deinterleave(encoded_bits, 30, 8)
    inv_deinterleaved_encoded_bits = invert_every_second_bit(deinterleaved_encoded_bits)
    decoded_bits = viterbi_decode(7, [79, 109], inv_deinterleaved_encoded_bits)
    bits = parse(UInt128, decoded_bits; base = 2)
    is_even = !get_bit(bits, 114, 1)
    is_nominal_page = !get_bit(bits, 114, 2)
    state = GNSSDecoderState(
        state;
        raw_data = GalileoE1BData(
            state.raw_data;
            num_pages_after_last_TOW = state.raw_data.num_pages_after_last_TOW + 1,
        ),
    )
    if is_even
        state = GNSSDecoderState(
            state;
            cache = GalileoE1BCache(
                state.cache;
                even_page_part_bits = is_nominal_page ? bits : UInt128(0),
            ),
        )
    elseif state.cache.even_page_part_bits != 0 && is_nominal_page
        data =
            get_bits(state.cache.even_page_part_bits, 114, 3, 112) << 16 +
            get_bits(bits, 114, 3, 16)
        bits_to_check_CRC =
            UInt288(state.cache.even_page_part_bits) << 106 + get_bits(bits, 114, 1, 106)
        if galCRC24(reverse(digits(UInt8, bits_to_check_CRC; base = 256))) == 0
            data_type = get_bits(data, 128, 1, 6)
            if data_type == 0
                if get_bits(data, 128, 7, 2) == 2 # '10'
                    WN = get_bits(data, 128, 97, 12)
                    TOW = get_bits(data, 128, 109, 20)
                    state = GNSSDecoderState(
                        state;
                        raw_data = GalileoE1BData(
                            state.raw_data;
                            WN,
                            TOW,
                            num_pages_after_last_TOW = 1,
                            num_bits_after_valid_syncro_sequence_after_last_TOW = state.num_bits_after_valid_syncro_sequence,
                        ),
                    )
                end
            elseif data_type == 1
                IOD_nav1 = get_bits(data, 128, 7, 10)
                t_0e = get_bits(data, 128, 17, 14) * 60
                M_0 =
                    get_twos_complement_num(data, 128, 31, 32) * state.constants.PI /
                    1 << 31
                e = get_bits(data, 128, 63, 32) / 1 << 33
                sqrt_A = get_bits(data, 128, 95, 32) / 1 << 19
                state = GNSSDecoderState(
                    state;
                    raw_data = GalileoE1BData(
                        state.raw_data;
                        IOD_nav1,
                        t_0e,
                        M_0,
                        e,
                        sqrt_A,
                    ),
                )
            elseif data_type == 2
                IOD_nav2 = get_bits(data, 128, 7, 10)
                Ω_0 =
                    get_twos_complement_num(data, 128, 17, 32) * state.constants.PI /
                    1 << 31
                i_0 =
                    get_twos_complement_num(data, 128, 49, 32) * state.constants.PI /
                    1 << 31
                ω =
                    get_twos_complement_num(data, 128, 81, 32) * state.constants.PI /
                    1 << 31
                i_dot =
                    get_twos_complement_num(data, 128, 113, 14) * state.constants.PI /
                    1 << 43
                state = GNSSDecoderState(
                    state;
                    raw_data = GalileoE1BData(state.raw_data; IOD_nav2, Ω_0, i_0, ω, i_dot),
                )
            elseif data_type == 3
                IOD_nav3 = get_bits(data, 128, 7, 10)
                Ω_dot =
                    get_twos_complement_num(data, 128, 17, 24) * state.constants.PI /
                    1 << 43
                Δn =
                    get_twos_complement_num(data, 128, 41, 16) * state.constants.PI /
                    1 << 43
                C_uc = get_twos_complement_num(data, 128, 57, 16) / 1 << 29
                C_us = get_twos_complement_num(data, 128, 73, 16) / 1 << 29
                C_rc = get_twos_complement_num(data, 128, 89, 16) / 1 << 5
                C_rs = get_twos_complement_num(data, 128, 105, 16) / 1 << 5
                SISA_e1_e5b = Int(get_bits(data, 128, 121, 8))
                state = GNSSDecoderState(
                    state;
                    raw_data = GalileoE1BData(
                        state.raw_data;
                        IOD_nav3,
                        Ω_dot,
                        Δn,
                        C_uc,
                        C_us,
                        C_rc,
                        C_rs,
                        SISA_e1_e5b,
                    ),
                )
            elseif data_type == 4
                IOD_nav4 = get_bits(data, 128, 7, 10)
                SVID = Int(get_bits(data, 128, 17, 6))
                C_ic = get_twos_complement_num(data, 128, 23, 16) / 1 << 29
                C_is = get_twos_complement_num(data, 128, 39, 16) / 1 << 29
                t_0c = get_bits(data, 128, 55, 14) * 60
                a_f0 = get_twos_complement_num(data, 128, 69, 31) / 1 << 34
                a_f1 = get_twos_complement_num(data, 128, 100, 21) / 1 << 46
                a_f2 = get_twos_complement_num(data, 128, 121, 6) / 1 << 59
                state = GNSSDecoderState(
                    state;
                    raw_data = GalileoE1BData(
                        state.raw_data;
                        IOD_nav4,
                        SVID,
                        C_ic,
                        C_is,
                        t_0c,
                        a_f0,
                        a_f1,
                        a_f2,
                    ),
                )
            elseif data_type == 5
                a_i0 = get_bits(data, 128, 7, 11) / 1 << 2
                a_i1 = get_twos_complement_num(data, 128, 18, 11) / 1 << 8
                a_i2 = get_twos_complement_num(data, 128, 29, 14) / 1 << 15
                iono_storm_flag_region1 = get_bit(data, 128, 43)
                iono_storm_flag_region2 = get_bit(data, 128, 44)
                iono_storm_flag_region3 = get_bit(data, 128, 45)
                iono_storm_flag_region4 = get_bit(data, 128, 46)
                iono_storm_flag_region5 = get_bit(data, 128, 47)
                broadcast_group_delay_e1_e5a =
                    get_twos_complement_num(data, 128, 48, 10) / 1 << 32
                broadcast_group_delay_e1_e5b =
                    get_twos_complement_num(data, 128, 58, 10) / 1 << 32
                signal_health_e5b = SignalHealth(get_bits(data, 128, 68, 2))
                signal_health_e1b = SignalHealth(get_bits(data, 128, 70, 2))
                data_validity_status_e5b = DataValidityStatus(get_bit(data, 128, 72))
                data_validity_status_e1b = DataValidityStatus(get_bit(data, 128, 73))
                WN = get_bits(data, 128, 74, 12)
                TOW = get_bits(data, 128, 86, 20)
                state = GNSSDecoderState(
                    state;
                    raw_data = GalileoE1BData(
                        state.raw_data;
                        a_i0,
                        a_i1,
                        a_i2,
                        iono_storm_flag_region1,
                        iono_storm_flag_region2,
                        iono_storm_flag_region3,
                        iono_storm_flag_region4,
                        iono_storm_flag_region5,
                        broadcast_group_delay_e1_e5a,
                        broadcast_group_delay_e1_e5b,
                        signal_health_e5b,
                        signal_health_e1b,
                        data_validity_status_e5b,
                        data_validity_status_e1b,
                        WN,
                        TOW,
                        num_pages_after_last_TOW = 1,
                        num_bits_after_valid_syncro_sequence_after_last_TOW = state.num_bits_after_valid_syncro_sequence,
                    ),
                )
            elseif data_type == 6
                A_0_utc = get_twos_complement_num(data, 128, 7, 32) / 1 << 30
                A_1_utc = get_twos_complement_num(data, 128, 39, 24) / 1 << 50
                Δt_LS = Int(get_twos_complement_num(data, 128, 63, 8))
                t_0t = Int(get_bits(data, 128, 71, 8) * 3600)
                WN_0t = Int(get_bits(data, 128, 79, 8))
                WN_LSF = Int(get_bits(data, 128, 87, 8))
                DN = Int(get_bits(data, 128, 95, 3))
                Δt_LSF = Int(get_twos_complement_num(data, 128, 98, 8))
                TOW = get_bits(data, 128, 106, 20)
                state = GNSSDecoderState(
                    state;
                    raw_data = GalileoE1BData(
                        state.raw_data;
                        A_0_utc,
                        A_1_utc,
                        Δt_LS,
                        t_0t,
                        WN_0t,
                        WN_LSF,
                        DN,
                        Δt_LSF,
                        TOW,
                        num_pages_after_last_TOW = 1,
                        num_bits_after_valid_syncro_sequence_after_last_TOW = state.num_bits_after_valid_syncro_sequence,
                    ),
                )
            elseif data_type == 7
                IOD_a = Int(get_bits(data, 128, 7, 4))
                WN_a = Int(get_bits(data, 128, 11, 2))
                t_0a = Int(get_bits(data, 128, 13, 10) * 600)
                SVID = Int(get_bits(data, 128, 23, 6))
                Δsqrt_A = get_twos_complement_num(data, 128, 29, 13) / 1 << 9
                e = get_bits(data, 128, 42, 11) / 1 << 16
                ω =
                    get_twos_complement_num(data, 128, 53, 16) * state.constants.PI /
                    1 << 15
                δi =
                    get_twos_complement_num(data, 128, 69, 11) * state.constants.PI /
                    1 << 14
                Ω_0 =
                    get_twos_complement_num(data, 128, 80, 16) * state.constants.PI /
                    1 << 15
                Ω_dot =
                    get_twos_complement_num(data, 128, 96, 11) * state.constants.PI /
                    1 << 33
                M_0 =
                    get_twos_complement_num(data, 128, 107, 16) * state.constants.PI /
                    1 << 15
                almanac_pos1 = GalileoAlmanac(;
                    SVID,
                    Δsqrt_A,
                    e,
                    ω,
                    δi,
                    Ω_0,
                    Ω_dot,
                    M_0,
                    IOD_a,
                    WN_a,
                    t_0a,
                )
                valid_SVID = SVID >= 1
                state = GNSSDecoderState(
                    state;
                    cache = GalileoE1BCache(
                        state.cache;
                        almanac_chain_pos1 = valid_SVID ? almanac_pos1 : GalileoAlmanac(),
                        almanac_chain_pos2 = GalileoAlmanac(),
                    ),
                )
            elseif data_type == 8
                IOD_a = Int(get_bits(data, 128, 7, 4))
                a_f0_pos1 = get_twos_complement_num(data, 128, 11, 16) / 1 << 19
                a_f1_pos1 = get_twos_complement_num(data, 128, 27, 13) / 1 << 38
                signal_health_e5b_pos1 = SignalHealth(get_bits(data, 128, 40, 2))
                signal_health_e1b_pos1 = SignalHealth(get_bits(data, 128, 42, 2))
                SVID = Int(get_bits(data, 128, 44, 6))
                Δsqrt_A = get_twos_complement_num(data, 128, 50, 13) / 1 << 9
                e = get_bits(data, 128, 63, 11) / 1 << 16
                ω =
                    get_twos_complement_num(data, 128, 74, 16) * state.constants.PI /
                    1 << 15
                δi =
                    get_twos_complement_num(data, 128, 90, 11) * state.constants.PI /
                    1 << 14
                Ω_0 =
                    get_twos_complement_num(data, 128, 101, 16) * state.constants.PI /
                    1 << 15
                Ω_dot =
                    get_twos_complement_num(data, 128, 117, 11) * state.constants.PI /
                    1 << 33
                # Flush position-1 almanac if its WT7 partial is intact and IODa matches
                completed_pos1 = GalileoAlmanac(
                    state.cache.almanac_chain_pos1;
                    a_f0 = a_f0_pos1,
                    a_f1 = a_f1_pos1,
                    signal_health_e5b = signal_health_e5b_pos1,
                    signal_health_e1b = signal_health_e1b_pos1,
                )
                almanacs = state.raw_data.almanacs
                if completed_pos1.IOD_a == IOD_a &&
                   !isnothing(completed_pos1.SVID)
                    almanacs =
                        isnothing(almanacs) ? Dictionary{Int,GalileoAlmanac}() :
                        copy(almanacs)
                    set!(almanacs, completed_pos1.SVID, completed_pos1)
                end
                almanac_pos2 = GalileoAlmanac(;
                    SVID,
                    Δsqrt_A,
                    e,
                    ω,
                    δi,
                    Ω_0,
                    Ω_dot,
                    IOD_a,
                )
                valid_SVID = SVID >= 1
                state = GNSSDecoderState(
                    state;
                    raw_data = GalileoE1BData(
                        state.raw_data;
                        almanacs,
                    ),
                    cache = GalileoE1BCache(
                        state.cache;
                        almanac_chain_pos1 = GalileoAlmanac(),
                        almanac_chain_pos2 = valid_SVID ? almanac_pos2 : GalileoAlmanac(),
                    ),
                )
            elseif data_type == 9
                IOD_a = Int(get_bits(data, 128, 7, 4))
                WN_a = Int(get_bits(data, 128, 11, 2))
                t_0a = Int(get_bits(data, 128, 13, 10) * 600)
                M_0_pos2 =
                    get_twos_complement_num(data, 128, 23, 16) * state.constants.PI /
                    1 << 15
                a_f0_pos2 = get_twos_complement_num(data, 128, 39, 16) / 1 << 19
                a_f1_pos2 = get_twos_complement_num(data, 128, 55, 13) / 1 << 38
                signal_health_e5b_pos2 = SignalHealth(get_bits(data, 128, 68, 2))
                signal_health_e1b_pos2 = SignalHealth(get_bits(data, 128, 70, 2))
                SVID = Int(get_bits(data, 128, 72, 6))
                Δsqrt_A = get_twos_complement_num(data, 128, 78, 13) / 1 << 9
                e = get_bits(data, 128, 91, 11) / 1 << 16
                ω =
                    get_twos_complement_num(data, 128, 102, 16) * state.constants.PI /
                    1 << 15
                δi =
                    get_twos_complement_num(data, 128, 118, 11) * state.constants.PI /
                    1 << 14
                completed_pos2 = GalileoAlmanac(
                    state.cache.almanac_chain_pos2;
                    M_0 = M_0_pos2,
                    a_f0 = a_f0_pos2,
                    a_f1 = a_f1_pos2,
                    signal_health_e5b = signal_health_e5b_pos2,
                    signal_health_e1b = signal_health_e1b_pos2,
                    WN_a,
                    t_0a,
                )
                almanacs = state.raw_data.almanacs
                if completed_pos2.IOD_a == IOD_a &&
                   !isnothing(completed_pos2.SVID)
                    almanacs =
                        isnothing(almanacs) ? Dictionary{Int,GalileoAlmanac}() :
                        copy(almanacs)
                    set!(almanacs, completed_pos2.SVID, completed_pos2)
                end
                almanac_pos3 = GalileoAlmanac(;
                    SVID,
                    Δsqrt_A,
                    e,
                    ω,
                    δi,
                    IOD_a,
                    WN_a,
                    t_0a,
                )
                valid_SVID = SVID >= 1
                state = GNSSDecoderState(
                    state;
                    raw_data = GalileoE1BData(
                        state.raw_data;
                        almanacs,
                    ),
                    cache = GalileoE1BCache(
                        state.cache;
                        almanac_chain_pos1 = valid_SVID ? almanac_pos3 : GalileoAlmanac(),
                        almanac_chain_pos2 = GalileoAlmanac(),
                    ),
                )
            elseif data_type == 10
                IOD_a = Int(get_bits(data, 128, 7, 4))
                Ω_0_pos3 =
                    get_twos_complement_num(data, 128, 11, 16) * state.constants.PI /
                    1 << 15
                Ω_dot_pos3 =
                    get_twos_complement_num(data, 128, 27, 11) * state.constants.PI /
                    1 << 33
                M_0_pos3 =
                    get_twos_complement_num(data, 128, 38, 16) * state.constants.PI /
                    1 << 15
                a_f0_pos3 = get_twos_complement_num(data, 128, 54, 16) / 1 << 19
                a_f1_pos3 = get_twos_complement_num(data, 128, 70, 13) / 1 << 38
                signal_health_e5b_pos3 = SignalHealth(get_bits(data, 128, 83, 2))
                signal_health_e1b_pos3 = SignalHealth(get_bits(data, 128, 85, 2))
                A_0G = get_twos_complement_num(data, 128, 87, 16) / 1 << 35
                A_1G = get_twos_complement_num(data, 128, 103, 12) / 1 << 51
                t_0G = Int(get_bits(data, 128, 115, 8) * 3600)
                WN_0G = Int(get_bits(data, 128, 123, 6))
                # Position 3 partial sits in cache pos1 (set by WT9)
                completed_pos3 = GalileoAlmanac(
                    state.cache.almanac_chain_pos1;
                    Ω_0 = Ω_0_pos3,
                    Ω_dot = Ω_dot_pos3,
                    M_0 = M_0_pos3,
                    a_f0 = a_f0_pos3,
                    a_f1 = a_f1_pos3,
                    signal_health_e5b = signal_health_e5b_pos3,
                    signal_health_e1b = signal_health_e1b_pos3,
                )
                almanacs = state.raw_data.almanacs
                if completed_pos3.IOD_a == IOD_a &&
                   !isnothing(completed_pos3.SVID)
                    almanacs =
                        isnothing(almanacs) ? Dictionary{Int,GalileoAlmanac}() :
                        copy(almanacs)
                    set!(almanacs, completed_pos3.SVID, completed_pos3)
                end
                state = GNSSDecoderState(
                    state;
                    raw_data = GalileoE1BData(
                        state.raw_data;
                        almanacs,
                        A_0G,
                        A_1G,
                        t_0G,
                        WN_0G,
                    ),
                    cache = GalileoE1BCache(
                        state.cache;
                        almanac_chain_pos1 = GalileoAlmanac(),
                        almanac_chain_pos2 = GalileoAlmanac(),
                    ),
                )
            elseif data_type == 16
                ΔA_red = get_twos_complement_num(data, 128, 7, 5) * 256.0
                e_x_red = get_twos_complement_num(data, 128, 12, 13) / 1 << 22
                e_y_red = get_twos_complement_num(data, 128, 25, 13) / 1 << 22
                Δi_0_red =
                    get_twos_complement_num(data, 128, 38, 17) * state.constants.PI /
                    1 << 22
                Ω_0_red =
                    get_twos_complement_num(data, 128, 55, 23) * state.constants.PI /
                    1 << 22
                λ_0_red =
                    get_twos_complement_num(data, 128, 78, 23) * state.constants.PI /
                    1 << 22
                a_f0_red = get_twos_complement_num(data, 128, 101, 22) / 1 << 26
                a_f1_red = get_twos_complement_num(data, 128, 123, 6) / 1 << 35
                reduced_ced = GalileoReducedCED(;
                    ΔA_red,
                    e_x_red,
                    e_y_red,
                    Δi_0_red,
                    Ω_0_red,
                    λ_0_red,
                    a_f0_red,
                    a_f1_red,
                )
                state = GNSSDecoderState(
                    state;
                    raw_data = GalileoE1BData(state.raw_data; reduced_ced),
                )
            end
        end
    end
    return state
end

function validate_data(state::GNSSDecoderState{<:GalileoE1BData})
    if is_decoding_completed_for_positioning(state.raw_data) &&
       state.raw_data.IOD_nav1 ==
       state.raw_data.IOD_nav2 ==
       state.raw_data.IOD_nav3 ==
       state.raw_data.IOD_nav4
        num_bits_after_valid_syncro_sequence = 0
        if state.data.TOW == state.raw_data.TOW
            num_bits_after_valid_syncro_sequence = state.num_bits_after_valid_syncro_sequence
        elseif !isnothing(state.raw_data.num_bits_after_valid_syncro_sequence_after_last_TOW)
            num_bits_after_valid_syncro_sequence = state.num_bits_after_valid_syncro_sequence -
                (state.raw_data.num_bits_after_valid_syncro_sequence_after_last_TOW -
                2 * state.constants.syncro_sequence_length - state.constants.preamble_length)
        else # first succesful decoding
            num_bits_after_valid_syncro_sequence = state.constants.preamble_length +
                (
                    state.raw_data.num_pages_after_last_TOW + 1
                ) * 250
        end
        state = GNSSDecoderState(
            state;
            data = state.raw_data,
            num_bits_after_valid_syncro_sequence,
        )
    end
    return state
end

"""
$(TYPEDSIGNATURES)

Check if the Galileo satellite is healthy and usable for positioning.

Examines both the signal health status (`signal_health_e1b`) and data validity
status (`data_validity_status_e1b`) from page type 5. A satellite is considered
healthy only if both conditions are met:
- Signal health is `signal_ok`
- Data validity is `navigation_data_valid`

!!! warning
    This function requires that page type 5 has been successfully decoded.
    Check that `state.data.signal_health_e1b` is not `nothing` before relying
    on this result.

# Arguments
- `state::GNSSDecoderState{<:GalileoE1BData}`: Galileo E1B decoder state with decoded data

# Returns
- `Bool`: `true` if satellite health and data validity indicate normal operation

# Example
```julia
state = GalileoE1BDecoderState(1)
state = decode(state, bits, num_bits)
if is_sat_healthy(state)
    # Safe to use for positioning
end
```

# See Also
- [`GalileoE1BDecoderState`](@ref): Create decoder state
- [`decode`](@ref): Decode navigation data
"""
function is_sat_healthy(state::GNSSDecoderState{<:GalileoE1BData})
    state.data.signal_health_e1b == signal_ok &&
        state.data.data_validity_status_e1b == navigation_data_valid
end
