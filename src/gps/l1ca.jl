# UInt320 buffer for GPS L1 (more efficient than UInt312)
# which holds at least a complete GPS L1 subframe plus
# 8 extra syncronization bíts
BitIntegers.@define_integers 320

"""
    GPSL1CAConstants

WGS 84 constants and LNAV message structure parameters for GPS L1 C/A signal decoding.

The physical constants are defined in IS-GPS-200 (Interface Specification) and are used
for computing satellite positions and clock corrections from broadcast ephemeris data.

# Fields

  - `syncro_sequence_length::Int`: Length of synchronization sequence in bits (300 bits = 10 words × 30 bits)
  - `preamble::UInt8`: TLM word preamble pattern (10001011 binary, 0x8B)
  - `preamble_length::Int`: Length of preamble in bits (8)
  - `word_length::Int`: Length of each LNAV word in bits (30)
  - `PI::Float64`: Mathematical constant π = 3.1415926535898 (IS-GPS-200 Table 20-IV)
  - `Ω_dot_e::Float64`: WGS 84 Earth rotation rate = 7.2921151467×10⁻⁵ rad/s
  - `c::Float64`: Speed of light = 2.99792458×10⁸ m/s
  - `μ::Float64`: WGS 84 Earth gravitational parameter = 3.986005×10¹⁴ m³/s²
  - `F::Float64`: Relativistic correction constant = -4.442807633×10⁻¹⁰ s/√m

# Reference

IS-GPS-200N, Section 20.3.3 and Table 20-IV
"""
Base.@kwdef struct GPSL1CAConstants <: AbstractGNSSConstants
    syncro_sequence_length::Int = 300
    preamble::UInt8 = 0b10001011
    preamble_length::Int = 8
    word_length::Int = 30
    PI::Float64 = GNSS_PI
    Ω_dot_e::Float64 = EARTH_ROTATION_RATE
    c::Float64 = SPEED_OF_LIGHT
    μ::Float64 = 3.986005e14
    F::Float64 = -4.442807633e-10
end

"""
    GPSL1CAAlmanac

Almanac data for one GPS satellite, decoded from a single LNAV almanac page
(subframe 4 pages 2-5 / 7-10 for SV 25-32, subframe 5 pages 1-24 for SV 1-24).

The almanac provides reduced-precision orbital and clock parameters for satellite
acquisition planning and bootstrapping position fixes. Inclination is broadcast as a
delta from the nominal GPS constellation value (`i_nominal = 0.3 semi-circles ≈ 54°`),
which the decoder stores in `δi`.

# Fields

  - `e::Float64`: Eccentricity (dimensionless)
  - `t_oa::Int`: Almanac reference time of week (seconds)
  - `δi::Float64`: Inclination delta from nominal 0.3 semi-circles (rad)
  - `Ω_dot::Float64`: Rate of right ascension (rad/s)
  - `sv_health::String`: 8-bit SV health status (binary string; "00000000" = healthy)
  - `sqrt_A::Float64`: Square root of semi-major axis (√m)
  - `Ω_0::Float64`: Longitude of ascending node at weekly epoch (rad)
  - `ω::Float64`: Argument of perigee (rad)
  - `M_0::Float64`: Mean anomaly at reference time (rad)
  - `a_f0::Float64`: SV clock bias correction coefficient (seconds)
  - `a_f1::Float64`: SV clock drift correction coefficient (s/s)

# Reference

IS-GPS-200N, Section 20.3.3.5.1.2, Table 20-VI
"""
Base.@kwdef struct GPSL1CAAlmanac
    e::Union{Nothing,Float64} = nothing
    t_oa::Union{Nothing,Int} = nothing
    δi::Union{Nothing,Float64} = nothing
    Ω_dot::Union{Nothing,Float64} = nothing
    sv_health::Union{Nothing,String} = nothing
    sqrt_A::Union{Nothing,Float64} = nothing
    Ω_0::Union{Nothing,Float64} = nothing
    ω::Union{Nothing,Float64} = nothing
    M_0::Union{Nothing,Float64} = nothing
    a_f0::Union{Nothing,Float64} = nothing
    a_f1::Union{Nothing,Float64} = nothing
end

"""
    GPSL1CAData

Decoded GPS L1 C/A LNAV navigation message data.

Contains ephemeris, clock correction, and satellite health parameters decoded from
subframes 1, 2, and 3 of the GPS LNAV message. All parameters conform to IS-GPS-200N.

# Telemetry and Handover Word (TLM/HOW) Fields

  - `last_subframe_id::Int`: ID of the last decoded subframe (1-5)
  - `integrity_status_flag::Bool`: LNAV data integrity status (0=OK, 1=bad)
  - `TOW::Int64`: Time of Week at message transmission (seconds, 0-604784)
  - `alert_flag::Bool`: URA may be worse than indicated (0=OK, 1=alert)
  - `anti_spoof_flag::Bool`: Anti-spoofing mode (0=off, 1=on)

# Subframe 1 - Clock Correction Parameters

  - `trans_week::Int64`: GPS week number (modulo 1024)
  - `codeonl2::Int64`: Code on L2 channel (0=invalid, 1=P-code, 2=C/A-code, 3=invalid)
  - `ura::Float64`: User Range Accuracy (meters), derived from URA index
  - `sv_health::String`: 6-bit satellite health status ("000000" = healthy)
  - `IODC::String`: Issue of Data, Clock (10-bit binary string)
  - `l2pcode::Bool`: L2 P-code data flag (1=LNAV OFF on P-code)
  - `T_GD::Float64`: L1-L2 group delay correction (seconds)
  - `t_0c::Int64`: Clock reference time (seconds)
  - `a_f0::Float64`: Clock bias correction coefficient (seconds)
  - `a_f1::Float64`: Clock drift correction coefficient (s/s)
  - `a_f2::Float64`: Clock drift rate correction coefficient (s/s²)

# Subframe 2 - Ephemeris Parameters (Part 1)

  - `IODE_Sub_2::String`: Issue of Data, Ephemeris from subframe 2 (8-bit binary string)
  - `C_rs::Float64`: Sine harmonic correction to orbit radius (meters)
  - `Δn::Float64`: Mean motion difference from computed value (semi-circles/s)
  - `M_0::Float64`: Mean anomaly at reference time (semi-circles)
  - `C_uc::Float64`: Cosine harmonic correction to argument of latitude (rad)
  - `e::Float64`: Eccentricity (dimensionless, range 0-0.03)
  - `C_us::Float64`: Sine harmonic correction to argument of latitude (rad)
  - `sqrt_A::Float64`: Square root of semi-major axis (√m)
  - `t_0e::Int64`: Ephemeris reference time (seconds)
  - `fit_interval::Bool`: Curve fit interval flag (0=4h, 1=>4h)
  - `AODO::Int64`: Age of Data Offset for NMCT (seconds)

# Subframe 3 - Ephemeris Parameters (Part 2)

  - `C_ic::Float64`: Cosine harmonic correction to inclination (rad)
  - `Ω_0::Float64`: Longitude of ascending node at weekly epoch (semi-circles)
  - `C_is::Float64`: Sine harmonic correction to inclination (rad)
  - `i_0::Float64`: Inclination angle at reference time (semi-circles)
  - `C_rc::Float64`: Cosine harmonic correction to orbit radius (meters)
  - `ω::Float64`: Argument of perigee (semi-circles)
  - `Ω_dot::Float64`: Rate of right ascension (semi-circles/s)
  - `IODE_Sub_3::String`: Issue of Data, Ephemeris from subframe 3 (8-bit binary string)
  - `i_dot::Float64`: Rate of inclination angle (semi-circles/s)

# Reference

IS-GPS-200N, Tables 20-I, 20-II, 20-III, Sections 20.3.3.3-20.3.3.4
"""
Base.@kwdef struct GPSL1CAData <: AbstractGNSSData
    last_subframe_id::Int = 0
    integrity_status_flag::Union{Nothing,Bool} = nothing
    TOW::Union{Nothing,Int64} = nothing
    alert_flag::Union{Nothing,Bool} = nothing
    anti_spoof_flag::Union{Nothing,Bool} = nothing

    trans_week::Union{Nothing,Int64} = nothing
    codeonl2::Union{Nothing,Int64} = nothing
    ura::Union{Nothing,Float64} = nothing
    sv_health::Union{Nothing,String} = nothing
    IODC::Union{Nothing,String} = nothing
    l2pcode::Union{Nothing,Bool} = nothing
    T_GD::Union{Nothing,Float64} = nothing
    t_0c::Union{Nothing,Int64} = nothing
    a_f2::Union{Nothing,Float64} = nothing
    a_f1::Union{Nothing,Float64} = nothing
    a_f0::Union{Nothing,Float64} = nothing

    IODE_Sub_2::Union{Nothing,String} = nothing
    C_rs::Union{Nothing,Float64} = nothing
    Δn::Union{Nothing,Float64} = nothing
    M_0::Union{Nothing,Float64} = nothing
    C_uc::Union{Nothing,Float64} = nothing
    e::Union{Nothing,Float64} = nothing
    C_us::Union{Nothing,Float64} = nothing
    sqrt_A::Union{Nothing,Float64} = nothing
    t_0e::Union{Nothing,Int64} = nothing
    fit_interval::Union{Nothing,Bool} = nothing
    AODO::Union{Nothing,Int64} = nothing

    C_ic::Union{Nothing,Float64} = nothing
    Ω_0::Union{Nothing,Float64} = nothing
    C_is::Union{Nothing,Float64} = nothing
    i_0::Union{Nothing,Float64} = nothing
    C_rc::Union{Nothing,Float64} = nothing
    ω::Union{Nothing,Float64} = nothing
    Ω_dot::Union{Nothing,Float64} = nothing
    IODE_Sub_3::Union{Nothing,String} = nothing
    i_dot::Union{Nothing,Float64} = nothing

    # Subframe 4 page 18: Ionospheric parameters
    α_0::Union{Nothing,Float64} = nothing
    α_1::Union{Nothing,Float64} = nothing
    α_2::Union{Nothing,Float64} = nothing
    α_3::Union{Nothing,Float64} = nothing
    β_0::Union{Nothing,Float64} = nothing
    β_1::Union{Nothing,Float64} = nothing
    β_2::Union{Nothing,Float64} = nothing
    β_3::Union{Nothing,Float64} = nothing

    # Subframe 4 page 18: UTC parameters
    A_0::Union{Nothing,Float64} = nothing
    A_1::Union{Nothing,Float64} = nothing
    Δt_LS::Union{Nothing,Int64} = nothing
    t_ot::Union{Nothing,Int64} = nothing
    WN_t::Union{Nothing,Int64} = nothing
    WN_LSF::Union{Nothing,Int64} = nothing
    DN::Union{Nothing,Int64} = nothing
    Δt_LSF::Union{Nothing,Int64} = nothing

    # Subframe 4 page 25: A-S flags and SV configurations (32 SVs)
    sv_config::Union{Nothing,Vector{Int64}} = nothing

    # Subframe 4 page 25: SV health for SV 25-32 (6-bit health words)
    sv_health_sf4_25::Union{Nothing,Vector{String}} = nothing

    # Subframe 5 pages 1-24: Almanac data (stored per SV ID)
    almanac::Union{Nothing,Dictionary{Int64,GPSL1CAAlmanac}} = nothing

    # Subframe 5 page 25: SV health for SV 1-24 (6-bit health words)
    sv_health_sf5_25::Union{Nothing,Vector{String}} = nothing

    # Subframe 5 page 25: Almanac reference time and week
    t_oa::Union{Nothing,Int64} = nothing
    WN_a::Union{Nothing,Int64} = nothing
end

function GPSL1CAData(
    data::GPSL1CAData;
    last_subframe_id = data.last_subframe_id,
    integrity_status_flag = data.integrity_status_flag,
    TOW = data.TOW,
    alert_flag = data.alert_flag,
    anti_spoof_flag = data.anti_spoof_flag,
    trans_week = data.trans_week,
    codeonl2 = data.codeonl2,
    ura = data.ura,
    sv_health = data.sv_health,
    IODC = data.IODC,
    l2pcode = data.l2pcode,
    T_GD = data.T_GD,
    t_0c = data.t_0c,
    a_f2 = data.a_f2,
    a_f1 = data.a_f1,
    a_f0 = data.a_f0,
    IODE_Sub_2 = data.IODE_Sub_2,
    C_rs = data.C_rs,
    Δn = data.Δn,
    M_0 = data.M_0,
    C_uc = data.C_uc,
    e = data.e,
    C_us = data.C_us,
    sqrt_A = data.sqrt_A,
    t_0e = data.t_0e,
    fit_interval = data.fit_interval,
    AODO = data.AODO,
    C_ic = data.C_ic,
    Ω_0 = data.Ω_0,
    C_is = data.C_is,
    i_0 = data.i_0,
    C_rc = data.C_rc,
    ω = data.ω,
    Ω_dot = data.Ω_dot,
    IODE_Sub_3 = data.IODE_Sub_3,
    i_dot = data.i_dot,
    α_0 = data.α_0,
    α_1 = data.α_1,
    α_2 = data.α_2,
    α_3 = data.α_3,
    β_0 = data.β_0,
    β_1 = data.β_1,
    β_2 = data.β_2,
    β_3 = data.β_3,
    A_0 = data.A_0,
    A_1 = data.A_1,
    Δt_LS = data.Δt_LS,
    t_ot = data.t_ot,
    WN_t = data.WN_t,
    WN_LSF = data.WN_LSF,
    DN = data.DN,
    Δt_LSF = data.Δt_LSF,
    sv_config = data.sv_config,
    sv_health_sf4_25 = data.sv_health_sf4_25,
    almanac = data.almanac,
    sv_health_sf5_25 = data.sv_health_sf5_25,
    t_oa = data.t_oa,
    WN_a = data.WN_a,
)
    GPSL1CAData(
        last_subframe_id,
        integrity_status_flag,
        TOW,
        alert_flag,
        anti_spoof_flag,
        trans_week,
        codeonl2,
        ura,
        sv_health,
        IODC,
        l2pcode,
        T_GD,
        t_0c,
        a_f2,
        a_f1,
        a_f0,
        IODE_Sub_2,
        C_rs,
        Δn,
        M_0,
        C_uc,
        e,
        C_us,
        sqrt_A,
        t_0e,
        fit_interval,
        AODO,
        C_ic,
        Ω_0,
        C_is,
        i_0,
        C_rc,
        ω,
        Ω_dot,
        IODE_Sub_3,
        i_dot,
        α_0,
        α_1,
        α_2,
        α_3,
        β_0,
        β_1,
        β_2,
        β_3,
        A_0,
        A_1,
        Δt_LS,
        t_ot,
        WN_t,
        WN_LSF,
        DN,
        Δt_LSF,
        sv_config,
        sv_health_sf4_25,
        almanac,
        sv_health_sf5_25,
        t_oa,
        WN_a,
    )
end

struct VotedGPSL1CAData
    vote::Int
    data::GPSL1CAData
end

"""
$(TYPEDEF)

Per-decoder cache for the GPS L1 C/A signal.

Holds the soft-symbol `CircularDeque{Float32}` (capacity = 300 + 8 = 308) and
the data-voting cache used by `confirm_data`. The struct itself is immutable.
Following the framework convention (see [`GNSSDecoderState`](@ref)), the only
field mutated in place is the shared soft-symbol buffer; the `old_data` voting
tally is rebuilt immutably and threaded through a new cache by `confirm_data`
(mirroring the Galileo E1B / GPS L1C-D caches). The packed-bit buffer used for
sync is *not* stored here — it is computed as a local value at sync time and
threaded through the decode path (see `pack_buffer`).

# Fields

$(TYPEDFIELDS)
"""
struct GPSL1CACache <: AbstractGNSSCache
    """
    Soft-symbol buffer (308 = 300 syncro + 8 preamble)
    """
    soft_buffer::CircularDeque{Float32}
    """
    Voting tally used by `confirm_data` for subframe-level data validation
    """
    old_data::Vector{VotedGPSL1CAData}
end

function GPSL1CACache()
    GPSL1CACache(CircularDeque{Float32}(308), Vector{VotedGPSL1CAData}())
end

function GPSL1CACache(old_data::Vector{VotedGPSL1CAData})
    # Convenience constructor used by the tests to seed the voting cache with a
    # known tally (fresh, empty soft buffer).
    GPSL1CACache(CircularDeque{Float32}(308), old_data)
end

# Keyword "rebuild" constructor, mirroring `GalileoE1BCache`. Reuses the shared
# soft-symbol buffer by reference and swaps in a freshly-built `old_data`
# tally, so `confirm_data` can thread a new cache through `GNSSDecoderState`
# instead of mutating the voting vector in place.
function GPSL1CACache(
    cache::GPSL1CACache;
    soft_buffer = cache.soft_buffer,
    old_data = cache.old_data,
)
    GPSL1CACache(soft_buffer, old_data)
end

function Base.:(==)(a::GPSL1CACache, b::GPSL1CACache)
    deques_equal(a.soft_buffer, b.soft_buffer) && a.old_data == b.old_data
end

packed_buffer_type(::GNSSDecoderState{<:GPSL1CAData}) = UInt320

function is_subframe1_decoded(data::GPSL1CAData)
    !isnothing(data.trans_week) &&
        !isnothing(data.codeonl2) &&
        !isnothing(data.ura) &&
        !isnothing(data.sv_health) &&
        !isnothing(data.IODC) &&
        !isnothing(data.l2pcode) &&
        !isnothing(data.T_GD) &&
        !isnothing(data.t_0c) &&
        !isnothing(data.a_f2) &&
        !isnothing(data.a_f1) &&
        !isnothing(data.a_f0)
end

function is_subframe2_decoded(data::GPSL1CAData)
    !isnothing(data.IODE_Sub_2) &&
        !isnothing(data.C_rs) &&
        !isnothing(data.Δn) &&
        !isnothing(data.M_0) &&
        !isnothing(data.C_uc) &&
        !isnothing(data.e) &&
        !isnothing(data.C_us) &&
        !isnothing(data.sqrt_A) &&
        !isnothing(data.t_0e) &&
        !isnothing(data.fit_interval) &&
        !isnothing(data.AODO)
end

function is_subframe3_decoded(data::GPSL1CAData)
    !isnothing(data.C_ic) &&
        !isnothing(data.Ω_0) &&
        !isnothing(data.C_is) &&
        !isnothing(data.i_0) &&
        !isnothing(data.C_rc) &&
        !isnothing(data.ω) &&
        !isnothing(data.Ω_dot) &&
        !isnothing(data.IODE_Sub_3) &&
        !isnothing(data.i_dot)
end

function is_subframe4_decoded(data::GPSL1CAData)
    # Subframe 4 is considered decoded when we have ionospheric and UTC parameters
    # (page 18) and SV configurations/health (page 25)
    !isnothing(data.α_0) &&
        !isnothing(data.α_1) &&
        !isnothing(data.α_2) &&
        !isnothing(data.α_3) &&
        !isnothing(data.β_0) &&
        !isnothing(data.β_1) &&
        !isnothing(data.β_2) &&
        !isnothing(data.β_3) &&
        !isnothing(data.A_0) &&
        !isnothing(data.A_1) &&
        !isnothing(data.Δt_LS) &&
        !isnothing(data.t_ot) &&
        !isnothing(data.WN_t) &&
        !isnothing(data.WN_LSF) &&
        !isnothing(data.DN) &&
        !isnothing(data.Δt_LSF) &&
        !isnothing(data.sv_config) &&
        !isnothing(data.sv_health_sf4_25)
end

function is_subframe5_decoded(data::GPSL1CAData)
    # Subframe 5 is considered decoded when we have SV health (page 25)
    # and almanac reference time/week
    !isnothing(data.sv_health_sf5_25) && !isnothing(data.t_oa) && !isnothing(data.WN_a)
end

function is_decoding_completed_for_positioning(data::GPSL1CAData)
    !isnothing(data.integrity_status_flag) &&
        !isnothing(data.TOW) &&
        !isnothing(data.alert_flag) &&
        !isnothing(data.anti_spoof_flag) &&
        is_subframe1_decoded(data) &&
        is_subframe2_decoded(data) &&
        is_subframe3_decoded(data)
end

"""
$(TYPEDSIGNATURES)

Create a decoder state for GPS L1 C/A navigation messages.

Initializes a [`GNSSDecoderState`](@ref) configured for decoding GPS L1 C/A
(Coarse/Acquisition) civil navigation messages. The decoder extracts ephemeris,
clock correction, and health data from the 50 bps LNAV data stream.

# Arguments

  - `prn::Int`: Pseudo-Random Noise code identifier (1-32 for GPS satellites)

# Returns

  - `GNSSDecoderState{GPSL1CAData}`: Initialized decoder state for GPS L1

# Example

```julia
state = GPSL1CADecoderState(1)  # Create decoder for PRN 1
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
function GPSL1CADecoderState(prn)
    GNSSDecoderState(
        prn,
        GPSL1CAData(),
        GPSL1CAData(),
        GPSL1CAConstants(),
        GPSL1CACache(),
        nothing,
        false,
    )
end

function GNSSDecoderState(system::GPSL1CA, prn)
    GPSL1CADecoderState(prn)
end

"""
$(TYPEDSIGNATURES)

Reset the GPS L1 decoder state after a signal loss or reacquisition.

Clears the bit buffers and time-of-week (TOW) field while preserving other
decoded ephemeris and clock data in `raw_data`. This allows faster recovery
after brief signal outages without requiring a full re-decode of all subframes.

!!! note

    The `trans_week` field is intentionally not reset as it is only broadcast
    in subframe 1. This may cause brief errors if a GPS week rollover occurs
    during a signal outage.

# Arguments

  - `state::GNSSDecoderState{<:GPSL1CAData}`: Current GPS L1 decoder state

# Returns

  - `GNSSDecoderState{<:GPSL1CAData}`: Reset decoder state with cleared buffers

# Example

```julia
# After detecting signal loss
state = reset_decoder_state(state)
# Continue decoding with preserved ephemeris
state = decode(state, new_bits, num_bits)
```

# See Also

    # Reset bit buffers and TOW data field, while keeping the

  - [`GPSL1CADecoderState`](@ref): Create a fresh decoder state    # remaining parameters in raw_data. This allows a GNSSReceiver
  - [`decode`](@ref): Continue decoding after reset    # to use a satellite after a reacquisition without waiting for
"""
function reset_decoder_state(state::GNSSDecoderState{<:GPSL1CAData})
    # Reset bit buffers and TOW data field, while keeping the
    # remaining parameters in raw_data. This allows a GNSSReceiver
    # to use a satellite after a reacquisition without waiting for
    # the decoding of all data fields.
    # Note: trans_week is currently not reset as it is only
    # broadcast in subframe 1 and thus may increase the time until
    # the decoder is available again after an outage. This will
    # lead to erroneous decoder information for a few seconds after
    # reacquisition when a new GPS week started during a signal outage.
    empty!(state.cache.soft_buffer)
    GNSSDecoderState(
        state;
        raw_data = GPSL1CAData(state.raw_data; TOW = nothing),
        data = GPSL1CAData(),
        num_bits_after_valid_syncro_sequence = nothing,
    )
end

function check_gpsl1_parity(word::Unsigned, prev_29 = false, prev_30 = false)
    function bit(bit_number)
        cbit = get_bit(word, 30, bit_number)
        prev_30 ? !cbit : cbit
    end
    # Parity check to verify the data integrity:
    D_25 =
        prev_29 ⊻ bit(1) ⊻ bit(2) ⊻ bit(3) ⊻ bit(5) ⊻ bit(6) ⊻ bit(10) ⊻ bit(11) ⊻ bit(12) ⊻
        bit(13) ⊻ bit(14) ⊻ bit(17) ⊻ bit(18) ⊻ bit(20) ⊻ bit(23)
    D_26 =
        prev_30 ⊻ bit(2) ⊻ bit(3) ⊻ bit(4) ⊻ bit(6) ⊻ bit(7) ⊻ bit(11) ⊻ bit(12) ⊻ bit(13) ⊻
        bit(14) ⊻ bit(15) ⊻ bit(18) ⊻ bit(19) ⊻ bit(21) ⊻ bit(24)
    D_27 =
        prev_29 ⊻ bit(1) ⊻ bit(3) ⊻ bit(4) ⊻ bit(5) ⊻ bit(7) ⊻ bit(8) ⊻ bit(12) ⊻ bit(13) ⊻
        bit(14) ⊻ bit(15) ⊻ bit(16) ⊻ bit(19) ⊻ bit(20) ⊻ bit(22)
    D_28 =
        prev_30 ⊻ bit(2) ⊻ bit(4) ⊻ bit(5) ⊻ bit(6) ⊻ bit(8) ⊻ bit(9) ⊻ bit(13) ⊻ bit(14) ⊻
        bit(15) ⊻ bit(16) ⊻ bit(17) ⊻ bit(20) ⊻ bit(21) ⊻ bit(23)
    D_29 =
        prev_30 ⊻ bit(1) ⊻ bit(3) ⊻ bit(5) ⊻ bit(6) ⊻ bit(7) ⊻ bit(9) ⊻ bit(10) ⊻ bit(14) ⊻
        bit(15) ⊻ bit(16) ⊻ bit(17) ⊻ bit(18) ⊻ bit(21) ⊻ bit(22) ⊻ bit(24)
    D_30 =
        prev_29 ⊻ bit(3) ⊻ bit(5) ⊻ bit(6) ⊻ bit(8) ⊻ bit(9) ⊻ bit(10) ⊻ bit(11) ⊻ bit(13) ⊻
        bit(15) ⊻ bit(19) ⊻ bit(22) ⊻ bit(23) ⊻ bit(24)
    computed_parity_bits =
        (
            (((D_25 << UInt(1) + D_26) << UInt(1) + D_27) << UInt(1) + D_28) << UInt(1) +
            D_29
        ) << UInt(1) + D_30
    computed_parity_bits == get_bits(word, 30, 25, 6)
end

function get_word(buffer, state::GNSSDecoderState{<:GPSL1CAData}, word_number::Int)
    num_words = Int(state.constants.syncro_sequence_length / state.constants.word_length)
    word =
        buffer >> UInt(
            state.constants.word_length * (num_words - word_number) +
            state.constants.preamble_length,
        )
    UInt(word & (UInt(1) << UInt(state.constants.word_length) - UInt(1)))
end

function can_decode_word(
    decode_bits::Function,
    state::GNSSDecoderState,
    buffer,
    word_number::Int,
)
    word = get_word(buffer, state, word_number)
    prev_word = get_word(buffer, state, word_number - 1)
    prev_word_bit_30 = get_bit(prev_word, 30, 30)
    is_checked_word = word_number != 1 && word_number != 3
    if check_gpsl1_parity(
        word,
        is_checked_word * get_bit(prev_word, 30, 29),
        is_checked_word * prev_word_bit_30,
    )
        word_comp = (is_checked_word * prev_word_bit_30) ? ~word : word
        data = decode_bits(word_comp, state)
        state = GNSSDecoderState(state; raw_data = data)
    end
    return state
end

function can_decode_two_words(
    decode_bits::Function,
    state::GNSSDecoderState,
    buffer,
    word_number1::Int,
    word_number2::Int,
)
    word1 = get_word(buffer, state, word_number1)
    prev_word1 = get_word(buffer, state, word_number1 - 1)
    prev_word1_bit_30 = get_bit(prev_word1, 30, 30)
    is_checked_word1 = word_number1 != 1 && word_number1 != 3
    word2 = get_word(buffer, state, word_number2)
    prev_word2 = get_word(buffer, state, word_number2 - 1)
    prev_word2_bit_30 = get_bit(prev_word2, 30, 30)
    is_checked_word2 = word_number2 != 1 && word_number2 != 3
    if check_gpsl1_parity(
        word1,
        is_checked_word1 * get_bit(prev_word1, 30, 29),
        is_checked_word1 * prev_word1_bit_30,
    ) && check_gpsl1_parity(
        word2,
        is_checked_word2 * get_bit(prev_word2, 30, 29),
        is_checked_word2 * prev_word2_bit_30,
    )
        word1_comp = (is_checked_word1 * prev_word1_bit_30) ? ~word1 : word1
        word2_comp = (is_checked_word2 * prev_word2_bit_30) ? ~word2 : word2
        data = decode_bits(word1_comp, word2_comp, state)
        state = GNSSDecoderState(state; raw_data = data)
    end
    return state
end

function read_tlm_and_how_words(state, buffer)
    state = can_decode_word(state, buffer, 1) do tlm_word, state
        integrity_status_flag = get_bit(tlm_word, 30, 23)
        GPSL1CAData(state.raw_data; integrity_status_flag)
    end
    prev_TOW = state.raw_data.TOW
    state = can_decode_word(state, buffer, 2) do how_word, state
        TOW = get_bits(how_word, 30, 1, 17) * 6
        alert_flag = get_bit(how_word, 30, 18)
        anti_spoof_flag = get_bit(how_word, 30, 19)
        last_subframe_id = get_bits(how_word, 30, 20, 3)
        GPSL1CAData(state.raw_data; last_subframe_id, TOW, alert_flag, anti_spoof_flag)
    end
    if !isnothing(prev_TOW) && prev_TOW + 1 != state.raw_data.TOW
        # Time of week must be decodable
        state =
            GNSSDecoderState(state; raw_data = GPSL1CAData(state.raw_data; TOW = nothing))
    end
    state
end

function decode_syncro_sequence(state::GNSSDecoderState{<:GPSL1CAData}, buffer)
    state = read_tlm_and_how_words(state, buffer)
    subframe_id = state.raw_data.last_subframe_id

    if subframe_id == 1
        state = can_decode_word(state, buffer, 3) do word3, state
            trans_week = get_bits(word3, 30, 1, 10)

            # Codes on L2 Channel
            codeonl2 = get_bits(word3, 30, 11, 2)

            # SV Accuracy, user range accuracy
            ura = get_bits(word3, 30, 13, 4)
            if ura <= 6
                ura = 2^(1 + (ura / 2))
            elseif 6 < ura <= 14
                ura = 2^(ura - 2)
            elseif ura == 15
                ura = nothing
            end

            # Satellite Health
            sv_health = bitstring(get_bits(word3, 30, 17, 6))[(end-5):end]
            if get_bit(word3, 30, 17)
                @warn "Bad LNAV Data, SV-Health critical", sv_health
            end

            GPSL1CAData(state.raw_data; trans_week, codeonl2, ura, sv_health)
        end

        state = can_decode_two_words(state, buffer, 3, 8) do word3, word8, state
            # Issue of Data Clock
            # 2 MSB in Word 2, LSB 8 in Word 8
            IODC =
                bitstring(get_bits(word3, 30, 23, 2))[(end-1):end] *
                bitstring(get_bits(word8, 30, 1, 8))[(end-7):end]
            GPSL1CAData(state.raw_data; IODC)
        end

        state = can_decode_word(state, buffer, 4) do word4, state
            # True: LNAV Datastream on PCode commanded OFF
            l2pcode = get_bit(word4, 30, 1)
            GPSL1CAData(state.raw_data; l2pcode)
        end

        state = can_decode_word(state, buffer, 7) do word7, state
            # group time differential
            T_GD = get_twos_complement_num(word7, 30, 17, 8) / 1 << 31
            GPSL1CAData(state.raw_data; T_GD)
        end

        state = can_decode_word(state, buffer, 8) do word8, state
            # Clock data reference
            t_0c = get_bits(word8, 30, 9, 16) << 4
            GPSL1CAData(state.raw_data; t_0c)
        end

        state = can_decode_word(state, buffer, 9) do word9, state
            # clock correction parameter a_f2
            a_f2 = get_twos_complement_num(word9, 30, 1, 8) / 1 << 55
            # clock correction parameter a_f1
            a_f1 = get_twos_complement_num(word9, 30, 9, 16) / 1 << 43
            GPSL1CAData(state.raw_data; a_f2, a_f1)
        end

        state = can_decode_word(state, buffer, 10) do word10, state
            # Clock data reference
            a_f0 = get_twos_complement_num(word10, 30, 1, 22) / 1 << 31
            GPSL1CAData(state.raw_data; a_f0)
        end
    elseif subframe_id == 2
        state = can_decode_word(state, buffer, 3) do word3, state
            # Issue of ephemeris data
            IODE_Sub_2 = bitstring(get_bits(word3, 30, 1, 8))[(end-7):end]
            C_rs = get_twos_complement_num(word3, 30, 9, 16) / 1 << 5
            GPSL1CAData(state.raw_data; IODE_Sub_2, C_rs)
        end

        state = can_decode_word(state, buffer, 4) do word4, state
            # Mean motion difference from computed value
            Δn = get_twos_complement_num(word4, 30, 1, 16) * state.constants.PI / 1 << 43
            GPSL1CAData(state.raw_data; Δn)
        end

        state = can_decode_two_words(state, buffer, 4, 5) do word4, word5, state
            # Mean motion difference from computed value
            combined_word =
                UInt(get_bits(word4, 30, 17, 8) << 24 + get_bits(word5, 30, 1, 24))
            M_0 =
                get_twos_complement_num(combined_word, 32, 1, 32) * state.constants.PI /
                1 << 31
            GPSL1CAData(state.raw_data; M_0)
        end

        state = can_decode_word(state, buffer, 6) do word6, state
            # Amplitude of the Cosine Harmonic Correction Term to the Argument Latitude
            C_uc = get_twos_complement_num(word6, 30, 1, 16) / 1 << 29
            GPSL1CAData(state.raw_data; C_uc)
        end

        state = can_decode_two_words(state, buffer, 6, 7) do word6, word7, state
            # Eccentricity
            e = (get_bits(word6, 30, 17, 8) << 24 + get_bits(word7, 30, 1, 24)) / 1 << 33
            GPSL1CAData(state.raw_data; e)
        end

        state = can_decode_word(state, buffer, 8) do word8, state
            # Amplitude of the Sine Harmonic Correction Term to the Argument of Latitude
            C_us = get_twos_complement_num(word8, 30, 1, 16) / 1 << 29
            GPSL1CAData(state.raw_data; C_us)
        end

        state = can_decode_two_words(state, buffer, 8, 9) do word8, word9, state
            # Square Root of Semi-Major Axis
            sqrt_A =
                (get_bits(word8, 30, 17, 8) << 24 + get_bits(word9, 30, 1, 24)) / 1 << 19
            GPSL1CAData(state.raw_data; sqrt_A)
        end

        state = can_decode_word(state, buffer, 10) do word10, state
            # Reference Time ephemeris
            t_0e = get_bits(word10, 30, 1, 16) << 4
            fit_interval = get_bit(word10, 30, 17)
            AODO = get_bits(word10, 30, 18, 5)
            GPSL1CAData(state.raw_data; t_0e, fit_interval, AODO)
        end
    elseif subframe_id == 3
        state = can_decode_word(state, buffer, 3) do word3, state
            # Amplitude of the Cosine Harmonic Correction to Angle of Inclination
            C_ic = get_twos_complement_num(word3, 30, 1, 16) / 1 << 29
            GPSL1CAData(state.raw_data; C_ic)
        end

        state = can_decode_two_words(state, buffer, 3, 4) do word3, word4, state
            # Longitude of Ascending Node of Orbit Plane at Weekly Epoch
            combined_word =
                UInt(get_bits(word3, 30, 17, 8) << 24 + get_bits(word4, 30, 1, 24))
            Ω_0 =
                get_twos_complement_num(combined_word, 32, 1, 32) * state.constants.PI /
                1 << 31
            GPSL1CAData(state.raw_data; Ω_0)
        end

        state = can_decode_word(state, buffer, 5) do word5, state
            # Amplitude of the sine harmonic correction term to angle of Inclination
            C_is = get_twos_complement_num(word5, 30, 1, 16) / 1 << 29
            GPSL1CAData(state.raw_data; C_is)
        end

        state = can_decode_two_words(state, buffer, 5, 6) do word5, word6, state
            # inclination Angle at reference time
            combined_word =
                UInt(get_bits(word5, 30, 17, 8) << 24 + get_bits(word6, 30, 1, 24))
            i_0 =
                get_twos_complement_num(combined_word, 32, 1, 32) * state.constants.PI /
                1 << 31
            GPSL1CAData(state.raw_data; i_0)
        end

        state = can_decode_word(state, buffer, 7) do word7, state
            # Amplitude of the cosine harmonic correction term to orbit Radius
            C_rc = get_twos_complement_num(word7, 30, 1, 16) / 1 << 5
            GPSL1CAData(state.raw_data; C_rc)
        end

        state = can_decode_two_words(state, buffer, 7, 8) do word7, word8, state
            # Argument of Perigee
            combined_word =
                UInt(get_bits(word7, 30, 17, 8) << 24 + get_bits(word8, 30, 1, 24))
            ω =
                get_twos_complement_num(combined_word, 32, 1, 32) * state.constants.PI /
                1 << 31
            GPSL1CAData(state.raw_data; ω)
        end

        state = can_decode_word(state, buffer, 9) do word9, state
            # Amplitude of the cosine harmonic correction term to orbit Radius
            Ω_dot = get_twos_complement_num(word9, 30, 1, 24) * state.constants.PI / 1 << 43
            GPSL1CAData(state.raw_data; Ω_dot)
        end

        state = can_decode_word(state, buffer, 10) do word10, state
            # Issue of Ephemeris Data
            IODE_Sub_3 = bitstring(get_bits(word10, 30, 1, 8))[(end-7):end]
            # Rate of Inclination Angle
            i_dot =
                get_twos_complement_num(word10, 30, 9, 14) * state.constants.PI / 1 << 43
            GPSL1CAData(state.raw_data; IODE_Sub_3, i_dot)
        end
    elseif subframe_id == 4
        # Get page ID (SV ID) from word 3 bits 3-8
        state = can_decode_word(state, buffer, 3) do word3, state
            sv_page_id = get_bits(word3, 30, 3, 6)
            GPSL1CAData(state.raw_data; last_subframe_id = 4 + sv_page_id * 100) # encode page in subframe_id temporarily
        end
        sv_page_id = (state.raw_data.last_subframe_id - 4) ÷ 100
        state = GNSSDecoderState(
            state;
            raw_data = GPSL1CAData(state.raw_data; last_subframe_id = 4),
        )

        if sv_page_id == 56 # Page 18: Ionospheric and UTC data
            state = decode_subframe4_page18(state, buffer)
        elseif sv_page_id == 63 # Page 25: A-S flags and SV health
            state = decode_subframe4_page25(state, buffer)
        elseif sv_page_id in 25:32 # Pages 2-5, 7-10: Almanac for SV 25-32
            state = decode_almanac_page(state, buffer, sv_page_id)
        end
    elseif subframe_id == 5
        # Get SV ID from word 3 bits 3-8
        state = can_decode_word(state, buffer, 3) do word3, state
            sv_id = get_bits(word3, 30, 3, 6)
            GPSL1CAData(state.raw_data; last_subframe_id = 5 + sv_id * 100) # encode SV ID temporarily
        end
        sv_id = (state.raw_data.last_subframe_id - 5) ÷ 100
        state = GNSSDecoderState(
            state;
            raw_data = GPSL1CAData(state.raw_data; last_subframe_id = 5),
        )

        if sv_id == 51 # Page 25: SV health and almanac reference
            state = decode_subframe5_page25(state, buffer)
        elseif sv_id in 1:24 # Pages 1-24: Almanac for SV 1-24
            state = decode_almanac_page(state, buffer, sv_id)
        end
    end

    return state
end

function decode_subframe4_page18(state::GNSSDecoderState{<:GPSL1CAData}, buffer)
    # Page 18 contains ionospheric parameters and UTC parameters
    # Word 3: bits 9-16 = α0, bits 17-24 = α1
    # Word 4: bits 1-8 = α2, bits 9-16 = α3, bits 17-24 = β0
    # Word 5: bits 1-8 = β1, bits 9-16 = β2, bits 17-24 = β3
    # Word 6: bits 1-24 = A1 (24 bits)
    # Word 7: bits 1-24 = A0 MSBs (24 bits)
    # Word 8: bits 1-8 = A0 LSBs (8 bits), bits 9-16 = tot, bits 17-24 = WNt
    # Word 9: bits 1-8 = ΔtLS, bits 9-16 = WNLSF, bits 17-24 = DN
    # Word 10: bits 1-8 = ΔtLSF

    state = can_decode_word(state, buffer, 3) do word3, state
        α_0 = get_twos_complement_num(word3, 30, 9, 8) / 1 << 30
        α_1 = get_twos_complement_num(word3, 30, 17, 8) / 1 << 27
        GPSL1CAData(state.raw_data; α_0, α_1)
    end

    state = can_decode_word(state, buffer, 4) do word4, state
        α_2 = get_twos_complement_num(word4, 30, 1, 8) / 1 << 24
        α_3 = get_twos_complement_num(word4, 30, 9, 8) / 1 << 24
        β_0 = get_twos_complement_num(word4, 30, 17, 8) * (1 << 11)
        GPSL1CAData(state.raw_data; α_2, α_3, β_0)
    end

    state = can_decode_word(state, buffer, 5) do word5, state
        β_1 = get_twos_complement_num(word5, 30, 1, 8) * (1 << 14)
        β_2 = get_twos_complement_num(word5, 30, 9, 8) * (1 << 16)
        β_3 = get_twos_complement_num(word5, 30, 17, 8) * (1 << 16)
        GPSL1CAData(state.raw_data; β_1, β_2, β_3)
    end

    state = can_decode_word(state, buffer, 6) do word6, state
        A_1 = get_twos_complement_num(word6, 30, 1, 24) / 1 << 50
        GPSL1CAData(state.raw_data; A_1)
    end

    state = can_decode_two_words(state, buffer, 7, 8) do word7, word8, state
        # A0 is 32 bits: 24 MSBs in word 7, 8 LSBs in word 8
        combined_word = UInt(get_bits(word7, 30, 1, 24) << 8 + get_bits(word8, 30, 1, 8))
        A_0 = get_twos_complement_num(combined_word, 32, 1, 32) / 1 << 30
        GPSL1CAData(state.raw_data; A_0)
    end

    state = can_decode_word(state, buffer, 8) do word8, state
        t_ot = get_bits(word8, 30, 9, 8) << 12
        WN_t = get_bits(word8, 30, 17, 8)
        GPSL1CAData(state.raw_data; t_ot, WN_t)
    end

    state = can_decode_word(state, buffer, 9) do word9, state
        Δt_LS = get_twos_complement_num(word9, 30, 1, 8)
        WN_LSF = get_bits(word9, 30, 9, 8)
        DN = get_bits(word9, 30, 17, 8)
        GPSL1CAData(state.raw_data; Δt_LS, WN_LSF, DN)
    end

    state = can_decode_word(state, buffer, 10) do word10, state
        Δt_LSF = get_twos_complement_num(word10, 30, 1, 8)
        GPSL1CAData(state.raw_data; Δt_LSF)
    end

    return state
end

function decode_subframe4_page25(state::GNSSDecoderState{<:GPSL1CAData}, buffer)
    # Page 25 contains A-S flags and SV configurations for 32 SVs
    # and SV health for SV 25-32
    # Word 3: bits 9-24 = 4 SVs config (4 bits each)
    # Words 4-7: 24 MSBs = 6 SVs config each (4 bits each)
    # Word 8: bits 1-16 = 4 SVs config, bits 19-24 = SV25 health (6 bits)
    # Word 9: bits 1-24 = SV26-29 health (6 bits each)
    # Word 10: bits 1-18 = SV30-32 health (6 bits each)

    # Decode SV configurations (32 x 4-bit values)
    sv_config = Vector{Int64}(undef, 32)

    state = can_decode_word(state, buffer, 3) do word3, state
        for i = 1:4
            sv_config[i] = get_bits(word3, 30, 9 + (i - 1) * 4, 4)
        end
        GPSL1CAData(state.raw_data; sv_config)
    end

    state = can_decode_word(state, buffer, 4) do word4, state
        cfg = something(state.raw_data.sv_config, Vector{Int64}(undef, 32))
        for i = 1:6
            cfg[4+i] = get_bits(word4, 30, 1 + (i - 1) * 4, 4)
        end
        GPSL1CAData(state.raw_data; sv_config = cfg)
    end

    state = can_decode_word(state, buffer, 5) do word5, state
        cfg = something(state.raw_data.sv_config, Vector{Int64}(undef, 32))
        for i = 1:6
            cfg[10+i] = get_bits(word5, 30, 1 + (i - 1) * 4, 4)
        end
        GPSL1CAData(state.raw_data; sv_config = cfg)
    end

    state = can_decode_word(state, buffer, 6) do word6, state
        cfg = something(state.raw_data.sv_config, Vector{Int64}(undef, 32))
        for i = 1:6
            cfg[16+i] = get_bits(word6, 30, 1 + (i - 1) * 4, 4)
        end
        GPSL1CAData(state.raw_data; sv_config = cfg)
    end

    state = can_decode_word(state, buffer, 7) do word7, state
        cfg = something(state.raw_data.sv_config, Vector{Int64}(undef, 32))
        for i = 1:6
            cfg[22+i] = get_bits(word7, 30, 1 + (i - 1) * 4, 4)
        end
        GPSL1CAData(state.raw_data; sv_config = cfg)
    end

    state = can_decode_word(state, buffer, 8) do word8, state
        cfg = something(state.raw_data.sv_config, Vector{Int64}(undef, 32))
        for i = 1:4
            cfg[28+i] = get_bits(word8, 30, 1 + (i - 1) * 4, 4)
        end
        # SV 25 health (6 bits) at bits 19-24
        sv_health_sf4_25 = Vector{String}(undef, 8)
        sv_health_sf4_25[1] = bitstring(get_bits(word8, 30, 19, 6))[(end-5):end]
        GPSL1CAData(state.raw_data; sv_config = cfg, sv_health_sf4_25)
    end

    state = can_decode_word(state, buffer, 9) do word9, state
        health = something(state.raw_data.sv_health_sf4_25, Vector{String}(undef, 8))
        for i = 1:4
            health[1+i] = bitstring(get_bits(word9, 30, 1 + (i - 1) * 6, 6))[(end-5):end]
        end
        GPSL1CAData(state.raw_data; sv_health_sf4_25 = health)
    end

    state = can_decode_word(state, buffer, 10) do word10, state
        health = something(state.raw_data.sv_health_sf4_25, Vector{String}(undef, 8))
        for i = 1:3
            health[5+i] = bitstring(get_bits(word10, 30, 1 + (i - 1) * 6, 6))[(end-5):end]
        end
        GPSL1CAData(state.raw_data; sv_health_sf4_25 = health)
    end

    return state
end

function decode_subframe5_page25(state::GNSSDecoderState{<:GPSL1CAData}, buffer)
    # Page 25 contains SV health for SV 1-24 and almanac reference time/week
    # Word 3: bits 9-16 = toa, bits 17-24 = WNa
    # Word 4: bits 1-24 = SV1-4 health (6 bits each)
    # Word 5-8: bits 1-24 = SV health (6 bits each, 4 SVs per word)
    # Word 9: bits 1-24 = SV21-24 health (6 bits each)

    state = can_decode_word(state, buffer, 3) do word3, state
        t_oa = get_bits(word3, 30, 9, 8) << 12
        WN_a = get_bits(word3, 30, 17, 8)
        GPSL1CAData(state.raw_data; t_oa, WN_a)
    end

    sv_health_sf5_25 = Vector{String}(undef, 24)

    state = can_decode_word(state, buffer, 4) do word4, state
        for i = 1:4
            sv_health_sf5_25[i] =
                bitstring(get_bits(word4, 30, 1 + (i - 1) * 6, 6))[(end-5):end]
        end
        GPSL1CAData(state.raw_data; sv_health_sf5_25)
    end

    state = can_decode_word(state, buffer, 5) do word5, state
        health = something(state.raw_data.sv_health_sf5_25, Vector{String}(undef, 24))
        for i = 1:4
            health[4+i] = bitstring(get_bits(word5, 30, 1 + (i - 1) * 6, 6))[(end-5):end]
        end
        GPSL1CAData(state.raw_data; sv_health_sf5_25 = health)
    end

    state = can_decode_word(state, buffer, 6) do word6, state
        health = something(state.raw_data.sv_health_sf5_25, Vector{String}(undef, 24))
        for i = 1:4
            health[8+i] = bitstring(get_bits(word6, 30, 1 + (i - 1) * 6, 6))[(end-5):end]
        end
        GPSL1CAData(state.raw_data; sv_health_sf5_25 = health)
    end

    state = can_decode_word(state, buffer, 7) do word7, state
        health = something(state.raw_data.sv_health_sf5_25, Vector{String}(undef, 24))
        for i = 1:4
            health[12+i] = bitstring(get_bits(word7, 30, 1 + (i - 1) * 6, 6))[(end-5):end]
        end
        GPSL1CAData(state.raw_data; sv_health_sf5_25 = health)
    end

    state = can_decode_word(state, buffer, 8) do word8, state
        health = something(state.raw_data.sv_health_sf5_25, Vector{String}(undef, 24))
        for i = 1:4
            health[16+i] = bitstring(get_bits(word8, 30, 1 + (i - 1) * 6, 6))[(end-5):end]
        end
        GPSL1CAData(state.raw_data; sv_health_sf5_25 = health)
    end

    state = can_decode_word(state, buffer, 9) do word9, state
        health = something(state.raw_data.sv_health_sf5_25, Vector{String}(undef, 24))
        for i = 1:4
            health[20+i] = bitstring(get_bits(word9, 30, 1 + (i - 1) * 6, 6))[(end-5):end]
        end
        GPSL1CAData(state.raw_data; sv_health_sf5_25 = health)
    end

    return state
end

function decode_almanac_page(state::GNSSDecoderState{<:GPSL1CAData}, buffer, sv_id::Int)
    # Almanac pages (subframe 4 pages 2-5, 7-10 for SV 25-32;
    #                subframe 5 pages 1-24 for SV 1-24)
    # Word 3: bits 9-24 = e (16 bits)
    # Word 4: bits 1-8 = toa, bits 9-24 = δi (16 bits)
    # Word 5: bits 1-16 = Ω_dot, bits 17-24 = SV health (8 bits)
    # Word 6: bits 1-24 = √A (24 bits)
    # Word 7: bits 1-24 = Ω0 (24 bits)
    # Word 8: bits 1-24 = ω (24 bits)
    # Word 9: bits 1-24 = M0 (24 bits)
    # Word 10: bits 1-8 = af0 MSBs (8 bits), bits 9-19 = af1 (11 bits), bits 20-22 = af0 LSBs (3 bits)

    if sv_id == 0
        # Dummy SV, skip decoding
        return state
    end

    # Decode all almanac words and build the almanac entry
    # We need to decode all words and check parity for each
    word3 = get_word(buffer, state, 3)
    word4 = get_word(buffer, state, 4)
    word5 = get_word(buffer, state, 5)
    word6 = get_word(buffer, state, 6)
    word7 = get_word(buffer, state, 7)
    word8 = get_word(buffer, state, 8)
    word9 = get_word(buffer, state, 9)
    word10 = get_word(buffer, state, 10)

    # Check parity for all words
    prev_word2 = get_word(buffer, state, 2)
    prev_word3 = get_word(buffer, state, 3)
    prev_word4 = get_word(buffer, state, 4)
    prev_word5 = get_word(buffer, state, 5)
    prev_word6 = get_word(buffer, state, 6)
    prev_word7 = get_word(buffer, state, 7)
    prev_word8 = get_word(buffer, state, 8)
    prev_word9 = get_word(buffer, state, 9)

    parity_ok =
        check_gpsl1_parity(
            word3,
            get_bit(prev_word2, 30, 29),
            get_bit(prev_word2, 30, 30),
        ) &&
        check_gpsl1_parity(
            word4,
            get_bit(prev_word3, 30, 29),
            get_bit(prev_word3, 30, 30),
        ) &&
        check_gpsl1_parity(
            word5,
            get_bit(prev_word4, 30, 29),
            get_bit(prev_word4, 30, 30),
        ) &&
        check_gpsl1_parity(
            word6,
            get_bit(prev_word5, 30, 29),
            get_bit(prev_word5, 30, 30),
        ) &&
        check_gpsl1_parity(
            word7,
            get_bit(prev_word6, 30, 29),
            get_bit(prev_word6, 30, 30),
        ) &&
        check_gpsl1_parity(
            word8,
            get_bit(prev_word7, 30, 29),
            get_bit(prev_word7, 30, 30),
        ) &&
        check_gpsl1_parity(
            word9,
            get_bit(prev_word8, 30, 29),
            get_bit(prev_word8, 30, 30),
        ) &&
        check_gpsl1_parity(word10, get_bit(prev_word9, 30, 29), get_bit(prev_word9, 30, 30))

    if !parity_ok
        return state
    end

    # Apply bit complement if needed based on prev word bit 30
    word3_comp = get_bit(prev_word2, 30, 30) ? ~word3 : word3
    word4_comp = get_bit(prev_word3, 30, 30) ? ~word4 : word4
    word5_comp = get_bit(prev_word4, 30, 30) ? ~word5 : word5
    word6_comp = get_bit(prev_word5, 30, 30) ? ~word6 : word6
    word7_comp = get_bit(prev_word6, 30, 30) ? ~word7 : word7
    word8_comp = get_bit(prev_word7, 30, 30) ? ~word8 : word8
    word9_comp = get_bit(prev_word8, 30, 30) ? ~word9 : word9
    word10_comp = get_bit(prev_word9, 30, 30) ? ~word10 : word10

    # Decode almanac parameters
    alm_e = get_bits(word3_comp, 30, 9, 16) / 1 << 21
    alm_toa = Int(get_bits(word4_comp, 30, 1, 8) << 12)
    alm_δi = get_twos_complement_num(word4_comp, 30, 9, 16) * state.constants.PI / 1 << 19
    alm_Ω_dot =
        get_twos_complement_num(word5_comp, 30, 1, 16) * state.constants.PI / 1 << 38
    alm_sv_health = bitstring(get_bits(word5_comp, 30, 17, 8))[(end-7):end]
    alm_sqrt_A = get_bits(word6_comp, 30, 1, 24) / 1 << 11
    alm_Ω_0 = get_twos_complement_num(word7_comp, 30, 1, 24) * state.constants.PI / 1 << 23
    alm_ω = get_twos_complement_num(word8_comp, 30, 1, 24) * state.constants.PI / 1 << 23
    alm_M_0 = get_twos_complement_num(word9_comp, 30, 1, 24) * state.constants.PI / 1 << 23

    # af0: 8 MSBs at bits 1-8, 3 LSBs at bits 20-22 = 11 bits total
    af0_msbs = get_bits(word10_comp, 30, 1, 8)
    af0_lsbs = get_bits(word10_comp, 30, 20, 3)
    af0_combined = UInt(af0_msbs << 3 + af0_lsbs)
    alm_af0 = get_twos_complement_num(af0_combined, 11, 1, 11) / 1 << 20
    alm_af1 = get_twos_complement_num(word10_comp, 30, 9, 11) / 1 << 38

    # Store almanac data for this SV
    almanac_entry = GPSL1CAAlmanac(;
        e = alm_e,
        t_oa = alm_toa,
        δi = alm_δi,
        Ω_dot = alm_Ω_dot,
        sv_health = alm_sv_health,
        sqrt_A = alm_sqrt_A,
        Ω_0 = alm_Ω_0,
        ω = alm_ω,
        M_0 = alm_M_0,
        a_f0 = alm_af0,
        a_f1 = alm_af1,
    )

    almanac = something(state.raw_data.almanac, Dictionary{Int64,GPSL1CAAlmanac}())
    set!(almanac, sv_id, almanac_entry)
    state = GNSSDecoderState(state; raw_data = GPSL1CAData(state.raw_data; almanac))

    return state
end

function compare_data(data::GPSL1CAData, new_data::GPSL1CAData)
    data.IODC == new_data.IODC && # IODE_Sub_2 and IODE_Sub_3 is already checked for in validate_data
        (
            data.TOW == new_data.TOW ||
            (data.TOW > new_data.TOW && data.trans_week < new_data.trans_week) ||
            (data.TOW < new_data.TOW && data.trans_week >= new_data.trans_week)
        ) &&
        data.T_GD == new_data.T_GD &&
        data.t_0c == new_data.t_0c &&
        data.a_f0 == new_data.a_f0 &&
        data.a_f1 == new_data.a_f1 &&
        data.a_f2 == new_data.a_f2 &&
        data.C_rs == new_data.C_rs &&
        data.Δn == new_data.Δn &&
        data.M_0 == new_data.M_0 &&
        data.C_uc == new_data.C_uc &&
        data.e == new_data.e &&
        data.C_us == new_data.C_us &&
        data.sqrt_A == new_data.sqrt_A &&
        data.t_0e == new_data.t_0e &&
        data.fit_interval == new_data.fit_interval &&
        data.AODO == new_data.AODO &&
        data.C_ic == new_data.C_ic &&
        data.Ω_0 == new_data.Ω_0 &&
        data.C_is == new_data.C_is &&
        data.i_0 == new_data.i_0 &&
        data.C_rc == new_data.C_rc &&
        data.ω == new_data.ω &&
        data.Ω_dot == new_data.Ω_dot &&
        data.i_dot == new_data.i_dot
end

function increment_voting(old_vote, max_vote)
    min(max_vote, old_vote + 1)
end

# Thread an updated voting tally through a new cache, reusing the shared
# soft-symbol buffer by reference. Keeps `confirm_data` free of in-place
# mutation of `cache.old_data`, matching the Galileo E1B / GPS L1C-D caches and
# the framework's immutable-threading convention (`gnss.jl`).
with_old_data(state, new_old_data; kwargs...) = GNSSDecoderState(
    state;
    cache = GPSL1CACache(state.cache; old_data = new_old_data),
    kwargs...,
)

function confirm_data(state, max_vote = 20)
    old_data = state.cache.old_data

    # Check if any entry has same IODC
    has_same_iodc = any(e.data.IODC == state.raw_data.IODC for e in old_data)

    # Find matching entry (same IODC and same data)
    matching_idx = findfirst(old_data) do entry
        entry.data.IODC == state.raw_data.IODC && compare_data(entry.data, state.raw_data)
    end

    if isnothing(matching_idx)
        if has_same_iodc
            # Same IODC exists but data differs - add as new entry, don't use data yet
            new_old_data = push!(copy(old_data), VotedGPSL1CAData(0, state.raw_data))
            return with_old_data(state, new_old_data; raw_data = GPSL1CAData())
        else
            # New IODC entirely
            if state.data == GPSL1CAData() # no data yet - add to cache and use data
                new_old_data = [VotedGPSL1CAData(0, state.raw_data)]
                return with_old_data(
                    state,
                    new_old_data;
                    data = state.raw_data,
                    num_bits_after_valid_syncro_sequence = state.constants.preamble_length,
                )
            else # add as new entry, don't use data yet
                new_old_data = push!(copy(old_data), VotedGPSL1CAData(0, state.raw_data))
                return with_old_data(state, new_old_data; raw_data = GPSL1CAData())
            end
        end
    end

    # Found matching entry - upvote it
    curr_score = old_data[matching_idx].vote
    new_vote = increment_voting(curr_score, max_vote)

    # Find best score among entries with same IODC
    best_score = maximum(e.vote for e in old_data if e.data.IODC == state.raw_data.IODC)

    if best_score > curr_score
        # Another entry has higher score - reject this data
        new_old_data = copy(old_data)
        new_old_data[matching_idx] = VotedGPSL1CAData(new_vote, state.raw_data)
        return with_old_data(state, new_old_data; raw_data = GPSL1CAData())
    end

    # This entry has the best (or tied best) score - use the data
    new_old_data = if new_vote == max_vote && length(old_data) > 1
        # Max votes reached - keep only this entry
        [VotedGPSL1CAData(new_vote, state.raw_data)]
    else
        updated = copy(old_data)
        updated[matching_idx] = VotedGPSL1CAData(new_vote, state.raw_data)
        updated
    end

    with_old_data(
        state,
        new_old_data;
        data = state.raw_data,
        num_bits_after_valid_syncro_sequence = state.constants.preamble_length,
    )
end

function validate_data(state::GNSSDecoderState{<:GPSL1CAData})
    if is_decoding_completed_for_positioning(state.raw_data) &&
       state.raw_data.IODC[3:10] == state.raw_data.IODE_Sub_2 == state.raw_data.IODE_Sub_3
        state = confirm_data(state)
    end
    return state
end

"""
$(TYPEDSIGNATURES)

Check if the GPS satellite is healthy and usable for positioning.

Examines the 6-bit satellite health field (`sv_health`) from subframe 1. A satellite
is considered healthy only if all health bits are zero (`"000000"`).

!!! warning

    This function requires that subframe 1 has been successfully decoded.
    Check that `state.data.sv_health` is not `nothing` before relying on this result.

# Arguments

  - `state::GNSSDecoderState{<:GPSL1CAData}`: GPS L1 decoder state with decoded data

# Returns

  - `Bool`: `true` if satellite health status indicates normal operation

# Example

```julia
state = GPSL1CADecoderState(1)
state = decode(state, bits, num_bits)
if is_sat_healthy(state)
    # Safe to use for positioning
end
```

# See Also

  - [`GPSL1CADecoderState`](@ref): Create decoder state
  - [`decode`](@ref): Decode navigation data
"""
function is_sat_healthy(state::GNSSDecoderState{<:GPSL1CAData})
    state.data.sv_health == "000000"
end
