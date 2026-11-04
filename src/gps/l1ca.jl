"""
    GPSL1CAConstants

LNAV message structure parameters for GPS L1 C/A signal decoding (IS-GPS-200).

# Fields

  - `sync_unit_length::Int`: Length of the sync unit in symbols (300 symbols = 10 words × 30 bits)
  - `preamble::UInt8`: TLM word preamble pattern (10001011 binary, 0x8B)
  - `preamble_length::Int`: Length of preamble in bits (8)
  - `word_length::Int`: Length of each LNAV word in bits (30)

# Reference

IS-GPS-200N, Section 20.3.3 and Table 20-IV
"""
Base.@kwdef struct GPSL1CAConstants <: AbstractGNSSConstants
    sync_unit_length::Int = 300
    preamble::UInt8 = 0b10001011
    preamble_length::Int = 8
    word_length::Int = 30
end

"""
    GPSL1CAAlmanac

Almanac data for one GPS satellite, decoded from a single LNAV almanac page
(subframe 4 pages 2-5 / 7-10 for SV 25-32, subframe 5 pages 1-24 for SV 1-24).

Reduced-precision orbital and clock parameters. Inclination is broadcast as a delta
from the nominal `i_nominal = 0.3 semi-circles ≈ 54°`, stored in `δi`.

# Fields

  - `e::Float64`: Eccentricity (dimensionless)
  - `t_0a::Int`: Almanac reference time of week (seconds)
  - `δi::Float64`: Inclination delta from nominal 0.3 semi-circles (rad)
  - `Ω_dot::Float64`: Rate of right ascension (rad/s)
  - `sv_health::Int64`: Raw 8-bit SV health word (0 = healthy; IS-GPS-200N
    20.3.3.5.1.3, Table 20-VII — the MSB summarises the navigation data and the
    five LSBs are a coded signal-component status, so the word is kept raw)
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
    t_0a::Union{Nothing,Int} = nothing
    δi::Union{Nothing,Float64} = nothing
    Ω_dot::Union{Nothing,Float64} = nothing
    sv_health::Union{Nothing,Int64} = nothing
    sqrt_A::Union{Nothing,Float64} = nothing
    Ω_0::Union{Nothing,Float64} = nothing
    ω::Union{Nothing,Float64} = nothing
    M_0::Union{Nothing,Float64} = nothing
    a_f0::Union{Nothing,Float64} = nothing
    a_f1::Union{Nothing,Float64} = nothing
end

"""
    GPSLNAVData

Decoded GPS L1 C/A LNAV navigation message data.

Ephemeris, clock correction, health, ionosphere, UTC and almanac parameters decoded
from the GPS LNAV message (IS-GPS-200N).

# Telemetry and Handover Word (TLM/HOW) Fields

  - `last_subframe_id::Int`: ID of the last decoded subframe (1-5)
  - `integrity_status_flag::Bool`: Integrity level of the URA (IS-GPS-200N 20.3.3.1).
    `false` is *legacy* (4.42 x URA bounds the instantaneous URE with probability
    1 - 1e-5 per hour); `true` is *enhanced* (5.73 x URA at 1 - 1e-8 per hour).
    `true` is not a fault; the fault indicator is `alert_flag`.
  - `TOW::Int64`: Time of Week at the start of the next subframe (seconds, 0-604794 in steps of 6)
  - `alert_flag::Bool`: URA may be worse than indicated (0=OK, 1=alert)
  - `anti_spoof_flag::Bool`: Anti-spoofing mode (0=off, 1=on)
  - `tow_anchor_bit_count::Int`: Bit-counter value when `TOW` was decoded

# Subframe 1 - Clock Correction Parameters

  - `WN::Int64`: GPS week number (modulo 1024)
  - `code_on_L2::Int64`: Code on L2 channel (0=invalid, 1=P-code, 2=C/A-code, 3=invalid)
  - `URA_index::Int64`: User Range Accuracy index, the raw 4-bit broadcast value
    (IS-GPS-200N 20.3.3.3.1.3, Table 20-I). Not converted to metres: the table maps
    it to a *bound*, and index 15 ("no accuracy prediction available - use at own
    risk") has no metre value.
  - `sv_health::Int64`: Raw 6-bit satellite health word (0 = healthy;
    IS-GPS-200N 20.3.3.3.1.4, Table 20-VIII)
  - `IODC::Int64`: Issue of Data, Clock (10 bits, 0-1023). An integer because the ICD
    defines arithmetic on it: `IODE == IODC & 0xff` (20.3.4.4), inequality against
    the preceding six hours, and the Table 20-XII range tests selecting the curve-fit
    interval (`IODC < 240`, `240-247`, `248-255, 496`)
  - `L2_P_data_flag::Bool`: L2 P-code data flag (1=LNAV OFF on P-code)
  - `T_GD::Float64`: L1-L2 group delay correction (seconds)
  - `t_0c::Int64`: Clock reference time (seconds; the ICD writes `toc`)
  - `a_f0::Float64`: Clock bias correction coefficient (seconds)
  - `a_f1::Float64`: Clock drift correction coefficient (s/s)
  - `a_f2::Float64`: Clock drift rate correction coefficient (s/s²)

# Subframe 2 - Ephemeris Parameters (Part 1)

  - `IODE_subframe2::Int64`: Issue of Data, Ephemeris from subframe 2 (8 bits)
  - `C_rs::Float64`: Sine harmonic correction to orbit radius (meters)
  - `Δn::Float64`: Mean motion difference from computed value (rad/s)
  - `M_0::Float64`: Mean anomaly at reference time (rad)
  - `C_uc::Float64`: Cosine harmonic correction to argument of latitude (rad)
  - `e::Float64`: Eccentricity (dimensionless, range 0-0.03)
  - `C_us::Float64`: Sine harmonic correction to argument of latitude (rad)
  - `sqrt_A::Float64`: Square root of semi-major axis (√m)
  - `t_0e::Int64`: Ephemeris reference time (seconds; the ICD writes `toe`)
  - `fit_interval::Bool`: Curve fit interval flag (0 = 4 h, 1 = more than 4 h)
  - `AODO::Int64`: Age of Data Offset for the NMCT (seconds; 5-bit broadcast
    count with an LSB of 900 s, so 27900 s means "NMCT unavailable")

# Subframe 3 - Ephemeris Parameters (Part 2)

  - `C_ic::Float64`: Cosine harmonic correction to inclination (rad)
  - `Ω_0::Float64`: Longitude of ascending node at weekly epoch (rad)
  - `C_is::Float64`: Sine harmonic correction to inclination (rad)
  - `i_0::Float64`: Inclination angle at reference time (rad)
  - `C_rc::Float64`: Cosine harmonic correction to orbit radius (meters)
  - `ω::Float64`: Argument of perigee (rad)
  - `Ω_dot::Float64`: Rate of right ascension (rad/s)
  - `IODE_subframe3::Int64`: Issue of Data, Ephemeris from subframe 3 (8 bits)
  - `i_dot::Float64`: Rate of inclination angle (rad/s)

# Subframe 4 Page 18 - Ionosphere and UTC

  - `α_0..α_3::Float64`: Klobuchar ionospheric delay coefficients (s, s/sc, s/sc², s/sc³)
  - `β_0..β_3::Float64`: Klobuchar ionospheric period coefficients (s, s/sc, s/sc², s/sc³)
  - `A_0UTC::Float64`: Constant term of the UTC polynomial (s; the ICD writes `A0`)
  - `A_1UTC::Float64`: 1st-order term of the UTC polynomial (s/s; the ICD writes `A1`)
  - `Δt_LS::Int64`: Leap-second count before the pending adjustment (s)
  - `t_0t::Int64`: UTC reference time of week (s; the ICD writes `tot`)
  - `WN_0t::Int64`: UTC reference week number (modulo 256; the ICD writes `WNt`)
  - `WN_LSF::Int64`: Week number of the leap-second adjustment (modulo 256)
  - `DN::Int64`: Day number at the end of which the leap second becomes effective (1-7)
  - `Δt_LSF::Int64`: Leap-second count after the pending adjustment (s)

# Subframe 4 Page 25 / Subframe 5 Page 25 - Configuration and health

  - `sv_config::Vector{Int64}`: 4-bit A-S flag and SV configuration per SV 1-32
  - `sv_health_sf4_25::Vector{Int64}`: Raw 6-bit health words for SV 25-32
  - `sv_health_sf5_25::Vector{Int64}`: Raw 6-bit health words for SV 1-24

# Subframe 5 Pages 1-24 - Almanac

  - `almanacs::SlotDictionary{GPSL1CAAlmanac,33}`: Per-SV almanacs, keyed by SV ID
  - `t_0a::Int64`: Almanac reference time of week (s; the ICD writes `toa`)
  - `WN_a::Int64`: Almanac reference week number (modulo 256)

# Reference

IS-GPS-200N, Tables 20-I, 20-II, 20-III, Sections 20.3.3.3-20.3.3.4
"""
Base.@kwdef struct GPSLNAVData <: AbstractGPSData
    last_subframe_id::Int = 0
    integrity_status_flag::Union{Nothing,Bool} = nothing
    TOW::Union{Nothing,Int64} = nothing
    alert_flag::Union{Nothing,Bool} = nothing
    anti_spoof_flag::Union{Nothing,Bool} = nothing
    # `num_bits_since_time_reference` when `TOW` was decoded (see `is_plausible_TOW`, #82).
    tow_anchor_bit_count::Union{Nothing,Int} = nothing

    WN::Union{Nothing,Int64} = nothing
    code_on_L2::Union{Nothing,Int64} = nothing
    URA_index::Union{Nothing,Int64} = nothing
    sv_health::Union{Nothing,Int64} = nothing
    IODC::Union{Nothing,Int64} = nothing
    L2_P_data_flag::Union{Nothing,Bool} = nothing
    T_GD::Union{Nothing,Float64} = nothing
    t_0c::Union{Nothing,Int64} = nothing
    a_f2::Union{Nothing,Float64} = nothing
    a_f1::Union{Nothing,Float64} = nothing
    a_f0::Union{Nothing,Float64} = nothing

    IODE_subframe2::Union{Nothing,Int64} = nothing
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
    IODE_subframe3::Union{Nothing,Int64} = nothing
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
    A_0UTC::Union{Nothing,Float64} = nothing
    A_1UTC::Union{Nothing,Float64} = nothing
    Δt_LS::Union{Nothing,Int64} = nothing
    t_0t::Union{Nothing,Int64} = nothing
    WN_0t::Union{Nothing,Int64} = nothing
    WN_LSF::Union{Nothing,Int64} = nothing
    DN::Union{Nothing,Int64} = nothing
    Δt_LSF::Union{Nothing,Int64} = nothing

    # Subframe 4 page 25: A-S flags and SV configurations (32 SVs)
    sv_config::Union{Nothing,Vector{Int64}} = nothing

    # Subframe 4 page 25: SV health for SV 25-32 (6-bit health words)
    sv_health_sf4_25::Union{Nothing,Vector{Int64}} = nothing

    # Subframe 5 pages 1-24: Almanac data (stored per SV ID)
    almanacs::Union{Nothing,SlotDictionary{GPSL1CAAlmanac,33}} = nothing

    # Subframe 5 page 25: SV health for SV 1-24 (6-bit health words)
    sv_health_sf5_25::Union{Nothing,Vector{Int64}} = nothing

    # Subframe 5 page 25: Almanac reference time and week
    t_0a::Union{Nothing,Int64} = nothing
    WN_a::Union{Nothing,Int64} = nothing
end

@inline function GPSLNAVData(
    data::GPSLNAVData;
    last_subframe_id = data.last_subframe_id,
    integrity_status_flag = data.integrity_status_flag,
    TOW = data.TOW,
    alert_flag = data.alert_flag,
    anti_spoof_flag = data.anti_spoof_flag,
    tow_anchor_bit_count = data.tow_anchor_bit_count,
    WN = data.WN,
    code_on_L2 = data.code_on_L2,
    URA_index = data.URA_index,
    sv_health = data.sv_health,
    IODC = data.IODC,
    L2_P_data_flag = data.L2_P_data_flag,
    T_GD = data.T_GD,
    t_0c = data.t_0c,
    a_f2 = data.a_f2,
    a_f1 = data.a_f1,
    a_f0 = data.a_f0,
    IODE_subframe2 = data.IODE_subframe2,
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
    IODE_subframe3 = data.IODE_subframe3,
    i_dot = data.i_dot,
    α_0 = data.α_0,
    α_1 = data.α_1,
    α_2 = data.α_2,
    α_3 = data.α_3,
    β_0 = data.β_0,
    β_1 = data.β_1,
    β_2 = data.β_2,
    β_3 = data.β_3,
    A_0UTC = data.A_0UTC,
    A_1UTC = data.A_1UTC,
    Δt_LS = data.Δt_LS,
    t_0t = data.t_0t,
    WN_0t = data.WN_0t,
    WN_LSF = data.WN_LSF,
    DN = data.DN,
    Δt_LSF = data.Δt_LSF,
    sv_config = data.sv_config,
    sv_health_sf4_25 = data.sv_health_sf4_25,
    almanacs = data.almanacs,
    sv_health_sf5_25 = data.sv_health_sf5_25,
    t_0a = data.t_0a,
    WN_a = data.WN_a,
)
    GPSLNAVData(
        last_subframe_id,
        integrity_status_flag,
        TOW,
        alert_flag,
        anti_spoof_flag,
        tow_anchor_bit_count,
        WN,
        code_on_L2,
        URA_index,
        sv_health,
        IODC,
        L2_P_data_flag,
        T_GD,
        t_0c,
        a_f2,
        a_f1,
        a_f0,
        IODE_subframe2,
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
        IODE_subframe3,
        i_dot,
        α_0,
        α_1,
        α_2,
        α_3,
        β_0,
        β_1,
        β_2,
        β_3,
        A_0UTC,
        A_1UTC,
        Δt_LS,
        t_0t,
        WN_0t,
        WN_LSF,
        DN,
        Δt_LSF,
        sv_config,
        sv_health_sf4_25,
        almanacs,
        sv_health_sf5_25,
        t_0a,
        WN_a,
    )
end

# The default `==` is `===`, which fails for the `Vector` and `Dictionary` fields
# (`sv_config`, `sv_health_sf*_25`, `almanacs`) even when contents match. Compare
# field-by-field (mirrors `GPSCNAVData`).
#
# This is equality of the *published* data only; subframe voting uses `compare_data`
# below, which weighs just the ephemeris fields.
Base.:(==)(a::GPSLNAVData, b::GPSLNAVData) = fields_equal(a, b)

struct VotedGPSLNAVData
    vote::Int
    data::GPSLNAVData
end

# Not implied by the `GPSLNAVData` method: the default `==` compares fields with
# `===`, so without this the wrapper, `GPSL1CACache.old_data` and hence the whole
# decoder state would compare by reference.
Base.:(==)(a::VotedGPSLNAVData, b::VotedGPSLNAVData) = fields_equal(a, b)

"""
$(TYPEDEF)

Per-decoder cache for the GPS L1 C/A signal.

Holds the soft-symbol `CircularDeque{Float32}` (capacity 300 + 8 = 308) and the
voting cache used by `confirm_data`. The struct is immutable; per the framework
convention (see [`GNSSDecoderState`](@ref)) its buffers are mutated in place: the
soft-symbol buffer, and the `old_data` tally by `confirm_data`. The packed-bit buffer
for sync is a local value, not stored here (see `pack_buffer`).

# Fields

$(TYPEDFIELDS)
"""
struct GPSL1CACache <: AbstractGNSSCache
    """
    Soft-symbol buffer (308 = 300 sync unit + 8 preamble)
    """
    soft_buffer::CircularDeque{Float32}
    """
    Voting tally used by `confirm_data` for subframe-level data validation,
    overwritten in place and capped at `GPSL1CA_MAX_VOTE_CANDIDATES` entries
    """
    old_data::Vector{VotedGPSLNAVData}
    """
    Preallocated containers `raw_data` and `data` are decoded into
    """
    storage::DataStorage{GPSLNAVData}
end

# Distinct IODC candidates the voting tally keeps. Data only changes at an upload, so
# more than two live candidates means bit errors that slipped past parity.
const GPSL1CA_MAX_VOTE_CANDIDATES = 8

function GPSL1CACache(old_data::Vector{VotedGPSLNAVData} = VotedGPSLNAVData[])
    # Seeding `old_data` lets tests start from a known tally.
    GPSL1CACache(
        CircularDeque{Float32}(308),
        sizehint!(copy(old_data), GPSL1CA_MAX_VOTE_CANDIDATES),
        DataStorage{GPSLNAVData}(),
    )
end

# Container fields at ICD size: 32 SV configs, 8 + 24 SV health words, an almanac
# slot per SV ID (1-32).
preallocated_data(::Type{GPSLNAVData}) = GPSLNAVData(;
    sv_config = zeros(Int64, 32),
    sv_health_sf4_25 = zeros(Int64, 8),
    almanacs = SlotDictionary{GPSL1CAAlmanac,33}(),
    sv_health_sf5_25 = zeros(Int64, 24),
)

function Base.:(==)(a::GPSL1CACache, b::GPSL1CACache)
    deques_equal(a.soft_buffer, b.soft_buffer) && a.old_data == b.old_data
end

packed_buffer_type(::GNSSDecoderState{<:GPSLNAVData}) = UInt320

function is_subframe1_decoded(data::GPSLNAVData)
    !isnothing(data.WN) &&
        !isnothing(data.code_on_L2) &&
        !isnothing(data.URA_index) &&
        !isnothing(data.sv_health) &&
        !isnothing(data.IODC) &&
        !isnothing(data.L2_P_data_flag) &&
        !isnothing(data.T_GD) &&
        !isnothing(data.t_0c) &&
        !isnothing(data.a_f2) &&
        !isnothing(data.a_f1) &&
        !isnothing(data.a_f0)
end

function is_subframe2_decoded(data::GPSLNAVData)
    !isnothing(data.IODE_subframe2) &&
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

function is_subframe3_decoded(data::GPSLNAVData)
    !isnothing(data.C_ic) &&
        !isnothing(data.Ω_0) &&
        !isnothing(data.C_is) &&
        !isnothing(data.i_0) &&
        !isnothing(data.C_rc) &&
        !isnothing(data.ω) &&
        !isnothing(data.Ω_dot) &&
        !isnothing(data.IODE_subframe3) &&
        !isnothing(data.i_dot)
end

# `TOW` is already in seconds (the 17-bit HOW count is scaled by 6 at parse).
get_broadcast_time_of_week(data::GPSLNAVData) = data.TOW

# LNAV carries UTC parameters but no inter-GNSS offset; GGTO is CNAV / CNAV-2 only
# (`gps/cnav.jl`, `gps/l1c_d.jl`).
get_time_offset(::GNSSDecoderState{<:GPSLNAVData}, ::TimeSystem) = nothing

function has_positioning_messages(data::GPSLNAVData)
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

Initializes a [`GNSSDecoderState`](@ref) for decoding GPS L1 C/A LNAV messages
(50 bps): ephemeris, clock correction, health, ionosphere, UTC and almanac data.

# Arguments

  - `prn::Int`: Pseudo-Random Noise code identifier (1-32 for GPS satellites)

# Returns

  - `GNSSDecoderState{GPSLNAVData}`: Initialized decoder state for GPS L1

# Example

```julia
state = GPSL1CADecoderState(1)  # Create decoder for PRN 1
state = decode!(state, bits, num_bits)
if is_sat_healthy(state)
    # Use state.data for positioning
end
```

# See Also

  - [`GNSSDecoderState`](@ref): The underlying state structure
  - [`decode!`](@ref): Decode bits using this state
  - [`reset_decoder_state!`](@ref): Reset after signal loss
  - [`is_sat_healthy`](@ref): Check satellite health status
"""
function GPSL1CADecoderState(prn)
    GNSSDecoderState(
        prn,
        GPSLNAVData(),
        GPSLNAVData(),
        GPSL1CAConstants(),
        GPSL1CACache(),
        nothing,
        false,
    )
end

function GNSSDecoderState(system::GPSL1CA, prn)
    GPSL1CADecoderState(prn)
end

# The signal this decoder demodulates. All GNSSSignals metadata (ids, band, nav-message
# rate, GPS Time) is forwarded through this mapping. Keyed on the constants type, which
# is 1:1 with the signal and tells apart decoders sharing a data type (see `src/gps/l5.jl`).
get_signal_type(::GPSL1CAConstants) = GPSL1CA

"""
$(TYPEDSIGNATURES)

Reset the GPS L1 decoder state after a signal loss or reacquisition.

Clears the bit buffers, the published `data` and `TOW`, while preserving other decoded
ephemeris and clock data in `raw_data`, so a satellite is usable again after
reacquisition without re-decoding all subframes; `is_sat_healthy` is `false` until the
data is validated again.

!!! note

    `WN` is deliberately not reset (it is only broadcast in subframe 1), so it may be
    wrong briefly if a GPS week crossover occurs during an outage.

# Arguments

  - `state::GNSSDecoderState{<:GPSLNAVData}`: Current GPS L1 decoder state

# Returns

  - `GNSSDecoderState{<:GPSLNAVData}`: Reset decoder state with cleared buffers

# Example

```julia
# After detecting signal loss
state = reset_decoder_state!(state)
# Continue decoding with preserved ephemeris
state = decode!(state, new_bits, num_bits)
```

# See Also

  - [`GPSL1CADecoderState`](@ref): Create a fresh decoder state
  - [`decode!`](@ref): Continue decoding after reset
"""
function reset_decoder_state!(state::GNSSDecoderState{<:GPSLNAVData})
    # Keep the other `raw_data` fields so a GNSSReceiver can use the satellite after
    # reacquisition without waiting for all data. `WN` is kept too: it is only in
    # subframe 1, and clearing it would delay availability (see docstring).
    empty!(state.cache.soft_buffer)
    GNSSDecoderState(
        state;
        raw_data = GPSLNAVData(
            state.raw_data;
            TOW = nothing,
            tow_anchor_bit_count = nothing,
        ),
        data = GPSLNAVData(),
        num_bits_since_time_reference = nothing,
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

function get_word(buffer, state::GNSSDecoderState{<:GPSLNAVData}, word_number::Int)
    num_words = Int(state.constants.sync_unit_length / state.constants.word_length)
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

# Largest legal TOW count in a HOW: steps once per 6 s subframe, rolls over after
# 100799 (IS-GPS-200N, Section 20.3.3.2), though the 17-bit field holds up to 131071.
const LNAV_MAX_TOW_COUNT = SECONDS_PER_WEEK ÷ 6 - 1

# An LNAV subframe is 300 symbols (6 s at 50 bps).
const LNAV_SYMBOLS_PER_SUBFRAME = 300

# Fallback bound for the TOW continuity screen while the symbol counter is not running
# (before the first data promotion, after a reset): largest forward step in seconds
# accepted between two decoded HOWs. Generous against real gaps (blockage, parity
# failures) yet small (~0.1% of the week) enough to reject nearly every false-lock
# TOW. A genuine gap beyond it costs one discarded TOW; the next is accepted fresh.
const LNAV_MAX_TOW_GAP = 600

"""
    is_plausible_TOW(TOW_count, prev_TOW, prev_TOW_anchor, num_bits_since_time_reference)

Screen a TOW count freshly decoded from a HOW against what a genuine frame lock can
produce, before it is stored (in seconds) in `raw_data.TOW`. `prev_TOW` is the TOW
(seconds) from an earlier HOW (`nothing` if none), `prev_TOW_anchor` the symbol-counter
value recorded then, and the last argument the current counter value
(`state.num_bits_since_time_reference`).

Frame sync is only a 16-bit preamble match, so a false lock on noise can reach this
point with an arbitrary 17-bit count behind valid parity (issue #82). Returns `true` if
the count passes the screens, strongest first:

  - The count must lie inside the week (at most `LNAV_MAX_TOW_COUNT`).
  - If the symbol counter ran across both HOW decodes, the elapsed symbols predict the
    TOW exactly: `elapsed` must be a positive multiple of `LNAV_SYMBOLS_PER_SUBFRAME`
    (genuine locks sit on the same 300-symbol grid) and the TOW must have advanced by
    6 s per subframe (modulo `SECONDS_PER_WEEK`). This holds across arbitrarily long
    gaps as long as symbols keep flowing, with no tuning constant. After a symbol slip,
    one genuine TOW is discarded and the next accepted fresh.
  - Without the counter, accept a forward step of at most `LNAV_MAX_TOW_GAP` seconds.
    With no `prev_TOW`, any in-week count is accepted.
"""
is_plausible_TOW(TOW_count, prev_TOW, prev_TOW_anchor, num_bits_since_time_reference) =
    is_plausible_TOW(
        TOW_count,
        prev_TOW,
        elapsed_symbols(prev_TOW_anchor, num_bits_since_time_reference),
    )

# The decoder calls this form: three `Union{Nothing,Int64}` arguments give eight
# union-split combinations, past inference's limit, so the four-argument call would
# dispatch dynamically, which `juliac --trim` rejects. Folding the two counter readings
# into `elapsed` leaves four.
function is_plausible_TOW(TOW_count, prev_TOW, elapsed)
    TOW_count <= LNAV_MAX_TOW_COUNT || return false
    isnothing(prev_TOW) && return true
    if !isnothing(elapsed)
        (elapsed > 0 && elapsed % LNAV_SYMBOLS_PER_SUBFRAME == 0) || return false
        return Int64(TOW_count) * 6 ==
               mod(prev_TOW + 6 * (elapsed ÷ LNAV_SYMBOLS_PER_SUBFRAME), SECONDS_PER_WEEK)
    end
    ΔTOW = mod(Int64(TOW_count) * 6 - prev_TOW, SECONDS_PER_WEEK)
    return 0 < ΔTOW <= LNAV_MAX_TOW_GAP
end

function read_tlm_and_how_words(state, buffer)
    state = can_decode_word(state, buffer, 1) do tlm_word, state
        integrity_status_flag = get_bit(tlm_word, 30, 23)
        GPSLNAVData(state.raw_data; integrity_status_flag)
    end
    # `raw_data.TOW` must only hold a TOW from *this* subframe's HOW: `confirm_data`
    # re-anchors `num_bits_since_time_reference` to it on promotion, so a stale TOW would
    # shift the transmit time by a multiple of 6 s. Clear it (and its anchor) up front;
    # a HOW passing parity and the plausibility screen writes them back.
    # Capture the concrete state, not three optional scalars: those give a different
    # closure type per Int/Nothing combination, and a missing HOW after the counter starts
    # would compile a new `can_decode_word` specialization during live tracking (292 ms
    # on an Orin, GNSSReceiver.jl#107).
    previous_state = state
    state = GNSSDecoderState(
        state;
        raw_data = GPSLNAVData(
            state.raw_data;
            TOW = nothing,
            tow_anchor_bit_count = nothing,
        ),
    )
    state = can_decode_word(state, buffer, 2) do how_word, state
        TOW_count = get_bits(how_word, 30, 1, 17)
        alert_flag = get_bit(how_word, 30, 18)
        anti_spoof_flag = get_bit(how_word, 30, 19)
        last_subframe_id = get_bits(how_word, 30, 20, 3)
        is_plausible = is_plausible_TOW(
            TOW_count,
            previous_state.raw_data.TOW,
            elapsed_symbols(
                previous_state.raw_data.tow_anchor_bit_count,
                previous_state.num_bits_since_time_reference,
            ),
        )
        TOW = is_plausible ? Int64(TOW_count) * 6 : nothing
        TOW_anchor =
            is_plausible ? previous_state.num_bits_since_time_reference : nothing
        @split_nothing (TOW, TOW_anchor) GPSLNAVData(
            state.raw_data;
            last_subframe_id,
            TOW,
            tow_anchor_bit_count = TOW_anchor,
            alert_flag,
            anti_spoof_flag,
        )
    end
    state
end

function decode_sync_unit(state::GNSSDecoderState{<:GPSLNAVData}, buffer)
    state = read_tlm_and_how_words(state, buffer)
    subframe_id = state.raw_data.last_subframe_id

    if subframe_id == 1
        state = can_decode_word(state, buffer, 3) do word3, state
            WN = get_bits(word3, 30, 1, 10)

            # Codes on L2 Channel
            code_on_L2 = get_bits(word3, 30, 11, 2)

            # Raw URA index, unconverted (Table 20-I maps it to a bound; 15 means "no
            # accuracy prediction"). Mapping 15 to `nothing` was indistinguishable
            # from "subframe 1 not decoded" and kept the SV out of
            # `has_positioning_data` forever.
            URA_index = Int64(get_bits(word3, 30, 13, 4))

            # Satellite Health
            sv_health = Int64(get_bits(word3, 30, 17, 6))
            if get_bit(word3, 30, 17)
                @warn "Bad LNAV Data, SV-Health critical", sv_health
            end

            GPSLNAVData(state.raw_data; WN, code_on_L2, URA_index, sv_health)
        end

        state = can_decode_two_words(state, buffer, 3, 8) do word3, word8, state
            # Issue of Data Clock
            # IODC: 2 MSBs in word 3 bits 23-24, 8 LSBs in word 8 (Figure 20-1)
            IODC = Int64(get_bits(word3, 30, 23, 2)) << 8 | Int64(get_bits(word8, 30, 1, 8))
            GPSLNAVData(state.raw_data; IODC)
        end

        state = can_decode_word(state, buffer, 4) do word4, state
            # True: LNAV Datastream on PCode commanded OFF
            L2_P_data_flag = get_bit(word4, 30, 1)
            GPSLNAVData(state.raw_data; L2_P_data_flag)
        end

        state = can_decode_word(state, buffer, 7) do word7, state
            # group time differential
            T_GD = get_twos_complement_num(word7, 30, 17, 8) * 2.0^-31
            GPSLNAVData(state.raw_data; T_GD)
        end

        state = can_decode_word(state, buffer, 8) do word8, state
            # Clock data reference
            t_0c = get_bits(word8, 30, 9, 16) << 4
            GPSLNAVData(state.raw_data; t_0c)
        end

        state = can_decode_word(state, buffer, 9) do word9, state
            # clock correction parameter a_f2
            a_f2 = get_twos_complement_num(word9, 30, 1, 8) * 2.0^-55
            # clock correction parameter a_f1
            a_f1 = get_twos_complement_num(word9, 30, 9, 16) * 2.0^-43
            GPSLNAVData(state.raw_data; a_f2, a_f1)
        end

        state = can_decode_word(state, buffer, 10) do word10, state
            # Clock data reference
            a_f0 = get_twos_complement_num(word10, 30, 1, 22) * 2.0^-31
            GPSLNAVData(state.raw_data; a_f0)
        end
    elseif subframe_id == 2
        state = can_decode_word(state, buffer, 3) do word3, state
            # Issue of ephemeris data
            IODE_subframe2 = Int64(get_bits(word3, 30, 1, 8))
            C_rs = get_twos_complement_num(word3, 30, 9, 16) / 1 << 5
            GPSLNAVData(state.raw_data; IODE_subframe2, C_rs)
        end

        state = can_decode_word(state, buffer, 4) do word4, state
            # Mean motion difference from computed value
            Δn = get_twos_complement_num(word4, 30, 1, 16) * GNSS_PI * 2.0^-43
            GPSLNAVData(state.raw_data; Δn)
        end

        state = can_decode_two_words(state, buffer, 4, 5) do word4, word5, state
            # Mean anomaly at reference time
            combined_word =
                UInt(get_bits(word4, 30, 17, 8) << 24 + get_bits(word5, 30, 1, 24))
            M_0 = get_twos_complement_num(combined_word, 32, 1, 32) * GNSS_PI * 2.0^-31
            GPSLNAVData(state.raw_data; M_0)
        end

        state = can_decode_word(state, buffer, 6) do word6, state
            # Amplitude of the Cosine Harmonic Correction Term to the Argument Latitude
            C_uc = get_twos_complement_num(word6, 30, 1, 16) / 1 << 29
            GPSLNAVData(state.raw_data; C_uc)
        end

        state = can_decode_two_words(state, buffer, 6, 7) do word6, word7, state
            # Eccentricity
            e = (get_bits(word6, 30, 17, 8) << 24 + get_bits(word7, 30, 1, 24)) * 2.0^-33
            GPSLNAVData(state.raw_data; e)
        end

        state = can_decode_word(state, buffer, 8) do word8, state
            # Amplitude of the Sine Harmonic Correction Term to the Argument of Latitude
            C_us = get_twos_complement_num(word8, 30, 1, 16) / 1 << 29
            GPSLNAVData(state.raw_data; C_us)
        end

        state = can_decode_two_words(state, buffer, 8, 9) do word8, word9, state
            # Square Root of Semi-Major Axis
            sqrt_A =
                (get_bits(word8, 30, 17, 8) << 24 + get_bits(word9, 30, 1, 24)) / 1 << 19
            GPSLNAVData(state.raw_data; sqrt_A)
        end

        state = can_decode_word(state, buffer, 10) do word10, state
            # Reference Time ephemeris
            t_0e = get_bits(word10, 30, 1, 16) << 4
            fit_interval = get_bit(word10, 30, 17)
            # Five-bit unsigned count with LSB 900 s, 0-31 (IS-GPS-200N 20.3.3.4); 11111
            # is 27900 s (20.3.3.4.4). The raw count is not the value.
            AODO = Int64(get_bits(word10, 30, 18, 5)) * 900
            GPSLNAVData(state.raw_data; t_0e, fit_interval, AODO)
        end
    elseif subframe_id == 3
        state = can_decode_word(state, buffer, 3) do word3, state
            # Amplitude of the Cosine Harmonic Correction to Angle of Inclination
            C_ic = get_twos_complement_num(word3, 30, 1, 16) / 1 << 29
            GPSLNAVData(state.raw_data; C_ic)
        end

        state = can_decode_two_words(state, buffer, 3, 4) do word3, word4, state
            # Longitude of Ascending Node of Orbit Plane at Weekly Epoch
            combined_word =
                UInt(get_bits(word3, 30, 17, 8) << 24 + get_bits(word4, 30, 1, 24))
            Ω_0 = get_twos_complement_num(combined_word, 32, 1, 32) * GNSS_PI * 2.0^-31
            GPSLNAVData(state.raw_data; Ω_0)
        end

        state = can_decode_word(state, buffer, 5) do word5, state
            # Amplitude of the sine harmonic correction term to angle of Inclination
            C_is = get_twos_complement_num(word5, 30, 1, 16) / 1 << 29
            GPSLNAVData(state.raw_data; C_is)
        end

        state = can_decode_two_words(state, buffer, 5, 6) do word5, word6, state
            # inclination Angle at reference time
            combined_word =
                UInt(get_bits(word5, 30, 17, 8) << 24 + get_bits(word6, 30, 1, 24))
            i_0 = get_twos_complement_num(combined_word, 32, 1, 32) * GNSS_PI * 2.0^-31
            GPSLNAVData(state.raw_data; i_0)
        end

        state = can_decode_word(state, buffer, 7) do word7, state
            # Amplitude of the cosine harmonic correction term to orbit Radius
            C_rc = get_twos_complement_num(word7, 30, 1, 16) / 1 << 5
            GPSLNAVData(state.raw_data; C_rc)
        end

        state = can_decode_two_words(state, buffer, 7, 8) do word7, word8, state
            # Argument of Perigee
            combined_word =
                UInt(get_bits(word7, 30, 17, 8) << 24 + get_bits(word8, 30, 1, 24))
            ω = get_twos_complement_num(combined_word, 32, 1, 32) * GNSS_PI * 2.0^-31
            GPSLNAVData(state.raw_data; ω)
        end

        state = can_decode_word(state, buffer, 9) do word9, state
            # Rate of right ascension
            Ω_dot = get_twos_complement_num(word9, 30, 1, 24) * GNSS_PI * 2.0^-43
            GPSLNAVData(state.raw_data; Ω_dot)
        end

        state = can_decode_word(state, buffer, 10) do word10, state
            # Issue of Ephemeris Data
            IODE_subframe3 = Int64(get_bits(word10, 30, 1, 8))
            # Rate of Inclination Angle
            i_dot = get_twos_complement_num(word10, 30, 9, 14) * GNSS_PI * 2.0^-43
            GPSLNAVData(state.raw_data; IODE_subframe3, i_dot)
        end
    elseif subframe_id == 4
        # Page ID (SV ID) is word 3 bits 3-8. The closure's local must not share a name
        # with `sv_page_id` below: a captured, reassigned variable is boxed as `Any`,
        # which `juliac --trim` rejects.
        state = can_decode_word(state, buffer, 3) do word3, state
            page_id = get_bits(word3, 30, 3, 6)
            GPSLNAVData(state.raw_data; last_subframe_id = 4 + page_id * 100) # encode page in subframe_id temporarily
        end
        sv_page_id = (state.raw_data.last_subframe_id - 4) ÷ 100
        state = GNSSDecoderState(
            state;
            raw_data = GPSLNAVData(state.raw_data; last_subframe_id = 4),
        )

        if sv_page_id == 56 # Page 18: Ionospheric and UTC data
            state = decode_subframe4_page18(state, buffer)
        elseif sv_page_id == 63 # Page 25: A-S flags and SV health
            state = decode_subframe4_page25(state, buffer)
        elseif sv_page_id in 25:32 # Pages 2-5, 7-10: Almanac for SV 25-32
            state = decode_almanac_page(state, buffer, sv_page_id)
        end
    elseif subframe_id == 5
        # SV ID is word 3 bits 3-8 (named apart from `sv_id`, as above)
        state = can_decode_word(state, buffer, 3) do word3, state
            word3_sv_id = get_bits(word3, 30, 3, 6)
            GPSLNAVData(state.raw_data; last_subframe_id = 5 + word3_sv_id * 100) # encode SV ID temporarily
        end
        sv_id = (state.raw_data.last_subframe_id - 5) ÷ 100
        state = GNSSDecoderState(
            state;
            raw_data = GPSLNAVData(state.raw_data; last_subframe_id = 5),
        )

        if sv_id == 51 # Page 25: SV health and almanac reference
            state = decode_subframe5_page25(state, buffer)
        elseif sv_id in 1:24 # Pages 1-24: Almanac for SV 1-24
            state = decode_almanac_page(state, buffer, sv_id)
        end
    end

    return state
end

function decode_subframe4_page18(state::GNSSDecoderState{<:GPSLNAVData}, buffer)
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
        GPSLNAVData(state.raw_data; α_0, α_1)
    end

    state = can_decode_word(state, buffer, 4) do word4, state
        α_2 = get_twos_complement_num(word4, 30, 1, 8) / 1 << 24
        α_3 = get_twos_complement_num(word4, 30, 9, 8) / 1 << 24
        β_0 = get_twos_complement_num(word4, 30, 17, 8) * (1 << 11)
        GPSLNAVData(state.raw_data; α_2, α_3, β_0)
    end

    state = can_decode_word(state, buffer, 5) do word5, state
        β_1 = get_twos_complement_num(word5, 30, 1, 8) * (1 << 14)
        β_2 = get_twos_complement_num(word5, 30, 9, 8) * (1 << 16)
        β_3 = get_twos_complement_num(word5, 30, 17, 8) * (1 << 16)
        GPSLNAVData(state.raw_data; β_1, β_2, β_3)
    end

    state = can_decode_word(state, buffer, 6) do word6, state
        A_1UTC = get_twos_complement_num(word6, 30, 1, 24) * 2.0^-50
        GPSLNAVData(state.raw_data; A_1UTC)
    end

    state = can_decode_two_words(state, buffer, 7, 8) do word7, word8, state
        # A0 is 32 bits: 24 MSBs in word 7, 8 LSBs in word 8
        combined_word = UInt(get_bits(word7, 30, 1, 24) << 8 + get_bits(word8, 30, 1, 8))
        A_0UTC = get_twos_complement_num(combined_word, 32, 1, 32) / 1 << 30
        GPSLNAVData(state.raw_data; A_0UTC)
    end

    state = can_decode_word(state, buffer, 8) do word8, state
        t_0t = get_bits(word8, 30, 9, 8) << 12
        WN_0t = get_bits(word8, 30, 17, 8)
        GPSLNAVData(state.raw_data; t_0t, WN_0t)
    end

    state = can_decode_word(state, buffer, 9) do word9, state
        Δt_LS = get_twos_complement_num(word9, 30, 1, 8)
        WN_LSF = get_bits(word9, 30, 9, 8)
        DN = get_bits(word9, 30, 17, 8)
        GPSLNAVData(state.raw_data; Δt_LS, WN_LSF, DN)
    end

    state = can_decode_word(state, buffer, 10) do word10, state
        Δt_LSF = get_twos_complement_num(word10, 30, 1, 8)
        GPSLNAVData(state.raw_data; Δt_LSF)
    end

    return state
end

function decode_subframe4_page25(state::GNSSDecoderState{<:GPSLNAVData}, buffer)
    # Page 25 contains A-S flags and SV configurations for 32 SVs
    # and SV health for SV 25-32
    # Word 3: bits 9-24 = 4 SVs config (4 bits each)
    # Words 4-7: 24 MSBs = 6 SVs config each (4 bits each)
    # Word 8: bits 1-16 = 4 SVs config, bits 19-24 = SV25 health (6 bits)
    # Word 9: bits 1-24 = SV26-29 health (6 bits each)
    # Word 10: bits 1-18 = SV30-32 health (6 bits each)

    # Decode SV configurations (32 x 4-bit values) and SV 25-32 health into the
    # preallocated tables, overwriting what an earlier page 25 left there.
    storage = state.cache.storage.raw
    sv_config = writable_container(state.raw_data.sv_config, storage.sv_config)
    sv_health_sf4_25 =
        writable_container(state.raw_data.sv_health_sf4_25, storage.sv_health_sf4_25)

    state = can_decode_word(state, buffer, 3) do word3, state
        for i = 1:4
            sv_config[i] = get_bits(word3, 30, 9 + (i - 1) * 4, 4)
        end
        GPSLNAVData(state.raw_data; sv_config)
    end

    state = can_decode_word(state, buffer, 4) do word4, state
        for i = 1:6
            sv_config[4+i] = get_bits(word4, 30, 1 + (i - 1) * 4, 4)
        end
        GPSLNAVData(state.raw_data; sv_config)
    end

    state = can_decode_word(state, buffer, 5) do word5, state
        for i = 1:6
            sv_config[10+i] = get_bits(word5, 30, 1 + (i - 1) * 4, 4)
        end
        GPSLNAVData(state.raw_data; sv_config)
    end

    state = can_decode_word(state, buffer, 6) do word6, state
        for i = 1:6
            sv_config[16+i] = get_bits(word6, 30, 1 + (i - 1) * 4, 4)
        end
        GPSLNAVData(state.raw_data; sv_config)
    end

    state = can_decode_word(state, buffer, 7) do word7, state
        for i = 1:6
            sv_config[22+i] = get_bits(word7, 30, 1 + (i - 1) * 4, 4)
        end
        GPSLNAVData(state.raw_data; sv_config)
    end

    state = can_decode_word(state, buffer, 8) do word8, state
        for i = 1:4
            sv_config[28+i] = get_bits(word8, 30, 1 + (i - 1) * 4, 4)
        end
        # SV 25 health (6 bits) at bits 19-24
        sv_health_sf4_25[1] = Int64(get_bits(word8, 30, 19, 6))
        GPSLNAVData(state.raw_data; sv_config, sv_health_sf4_25)
    end

    state = can_decode_word(state, buffer, 9) do word9, state
        for i = 1:4
            sv_health_sf4_25[1+i] = Int64(get_bits(word9, 30, 1 + (i - 1) * 6, 6))
        end
        GPSLNAVData(state.raw_data; sv_health_sf4_25)
    end

    state = can_decode_word(state, buffer, 10) do word10, state
        for i = 1:3
            sv_health_sf4_25[5+i] = Int64(get_bits(word10, 30, 1 + (i - 1) * 6, 6))
        end
        GPSLNAVData(state.raw_data; sv_health_sf4_25)
    end

    return state
end

function decode_subframe5_page25(state::GNSSDecoderState{<:GPSLNAVData}, buffer)
    # Page 25 contains SV health for SV 1-24 and almanac reference time/week
    # Word 3: bits 9-16 = toa, bits 17-24 = WNa
    # Word 4: bits 1-24 = SV1-4 health (6 bits each)
    # Word 5-8: bits 1-24 = SV health (6 bits each, 4 SVs per word)
    # Word 9: bits 1-24 = SV21-24 health (6 bits each)

    state = can_decode_word(state, buffer, 3) do word3, state
        t_0a = get_bits(word3, 30, 9, 8) << 12
        WN_a = get_bits(word3, 30, 17, 8)
        GPSLNAVData(state.raw_data; t_0a, WN_a)
    end

    sv_health_sf5_25 = writable_container(
        state.raw_data.sv_health_sf5_25,
        state.cache.storage.raw.sv_health_sf5_25,
    )

    state = can_decode_word(state, buffer, 4) do word4, state
        for i = 1:4
            sv_health_sf5_25[i] = Int64(get_bits(word4, 30, 1 + (i - 1) * 6, 6))
        end
        GPSLNAVData(state.raw_data; sv_health_sf5_25)
    end

    state = can_decode_word(state, buffer, 5) do word5, state
        for i = 1:4
            sv_health_sf5_25[4+i] = Int64(get_bits(word5, 30, 1 + (i - 1) * 6, 6))
        end
        GPSLNAVData(state.raw_data; sv_health_sf5_25)
    end

    state = can_decode_word(state, buffer, 6) do word6, state
        for i = 1:4
            sv_health_sf5_25[8+i] = Int64(get_bits(word6, 30, 1 + (i - 1) * 6, 6))
        end
        GPSLNAVData(state.raw_data; sv_health_sf5_25)
    end

    state = can_decode_word(state, buffer, 7) do word7, state
        for i = 1:4
            sv_health_sf5_25[12+i] = Int64(get_bits(word7, 30, 1 + (i - 1) * 6, 6))
        end
        GPSLNAVData(state.raw_data; sv_health_sf5_25)
    end

    state = can_decode_word(state, buffer, 8) do word8, state
        for i = 1:4
            sv_health_sf5_25[16+i] = Int64(get_bits(word8, 30, 1 + (i - 1) * 6, 6))
        end
        GPSLNAVData(state.raw_data; sv_health_sf5_25)
    end

    state = can_decode_word(state, buffer, 9) do word9, state
        for i = 1:4
            sv_health_sf5_25[20+i] = Int64(get_bits(word9, 30, 1 + (i - 1) * 6, 6))
        end
        GPSLNAVData(state.raw_data; sv_health_sf5_25)
    end

    return state
end

function decode_almanac_page(state::GNSSDecoderState{<:GPSLNAVData}, buffer, sv_id::Int)
    # Almanac pages (subframe 4 pages 2-5, 7-10 for SV 25-32;
    #                subframe 5 pages 1-24 for SV 1-24)
    # Word 3: bits 9-24 = e (16 bits)
    # Word 4: bits 1-8 = toa, bits 9-24 = δi (16 bits)
    # Word 5: bits 1-16 = Ω_dot, bits 17-24 = SV health (8 bits)
    # Word 6: bits 1-24 = √A (24 bits)
    # Word 7: bits 1-24 = Ω0 (24 bits)
    # Word 8: bits 1-24 = ω (24 bits)
    # Word 9: bits 1-24 = M0 (24 bits)
    # Word 10: bits 1-8 = af0 MSBs (8 bits), bits 9-19 = af1 (11 bits),
    #          bits 20-22 = af0 LSBs (3 bits)

    if sv_id == 0
        # Dummy SV, skip decoding
        return state
    end

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

    # Complement per the previous word's bit 30
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
    alm_δi = get_twos_complement_num(word4_comp, 30, 9, 16) * GNSS_PI / 1 << 19
    alm_Ω_dot = get_twos_complement_num(word5_comp, 30, 1, 16) * GNSS_PI * 2.0^-38
    alm_sv_health = Int64(get_bits(word5_comp, 30, 17, 8))
    alm_sqrt_A = get_bits(word6_comp, 30, 1, 24) / 1 << 11
    alm_Ω_0 = get_twos_complement_num(word7_comp, 30, 1, 24) * GNSS_PI / 1 << 23
    alm_ω = get_twos_complement_num(word8_comp, 30, 1, 24) * GNSS_PI / 1 << 23
    alm_M_0 = get_twos_complement_num(word9_comp, 30, 1, 24) * GNSS_PI / 1 << 23

    # af0: 8 MSBs at bits 1-8, 3 LSBs at bits 20-22 = 11 bits total
    af0_msbs = get_bits(word10_comp, 30, 1, 8)
    af0_lsbs = get_bits(word10_comp, 30, 20, 3)
    af0_combined = UInt(af0_msbs << 3 + af0_lsbs)
    alm_af0 = get_twos_complement_num(af0_combined, 11, 1, 11) / 1 << 20
    alm_af1 = get_twos_complement_num(word10_comp, 30, 9, 11) * 2.0^-38

    # Store almanac data for this SV
    almanac_entry = GPSL1CAAlmanac(;
        e = alm_e,
        t_0a = alm_toa,
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

    # Overwrites this SV's slot of the preallocated almanac store in place.
    almanacs = writable_container(state.raw_data.almanacs, state.cache.storage.raw.almanacs)
    set!(almanacs, sv_id, almanac_entry)
    state = GNSSDecoderState(state; raw_data = GPSLNAVData(state.raw_data; almanacs))

    return state
end

function compare_data(data::GPSLNAVData, new_data::GPSLNAVData)
    data.IODC == new_data.IODC && # IODE_subframe2 and IODE_subframe3 is already checked for in validate_data
        (
            data.TOW == new_data.TOW ||
            (data.TOW > new_data.TOW && data.WN < new_data.WN) ||
            (data.TOW < new_data.TOW && data.WN >= new_data.WN)
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

# Promote `raw_data` to validated `data`, re-anchoring the symbol counter to the
# current subframe (`preamble_length` symbols past its boundary). The TOW anchor is
# rebased into the same frame: promotion only happens with a TOW decoded in this
# subframe (`read_tlm_and_how_words` clears stale ones).
#
# `data` gets its own copy of every container (`publish_data`), so later subframes
# written into `raw_data`'s containers do not leak into `data` before being voted on.
function promote_data(state)
    promoted =
        GPSLNAVData(state.raw_data; tow_anchor_bit_count = state.constants.preamble_length)
    GNSSDecoderState(
        state;
        raw_data = promoted,
        data = publish_data(state.cache.storage, promoted),
        num_bits_since_time_reference = state.constants.preamble_length,
    )
end

# The voting tally `state.cache.old_data` is overwritten in place throughout.
function confirm_data(state, max_vote = 20)
    old_data = state.cache.old_data

    # Any entry with the same IODC?
    has_same_iodc = any(e.data.IODC == state.raw_data.IODC for e in old_data)

    # Entry with same IODC and same data
    matching_idx = findfirst(old_data) do entry
        entry.data.IODC == state.raw_data.IODC && compare_data(entry.data, state.raw_data)
    end

    if isnothing(matching_idx)
        if !has_same_iodc && state.data == GPSLNAVData()
            # New IODC and no data yet - add to cache and use data
            replace_vote_candidates!(old_data, VotedGPSLNAVData(0, state.raw_data))
            return promote_data(state)
        end
        # Same IODC with differing data, or a new IODC while data is in use: add as a
        # new entry, don't publish yet
        push_vote_candidate!(
            old_data,
            VotedGPSLNAVData(0, state.raw_data),
            GPSL1CA_MAX_VOTE_CANDIDATES,
        )
        return GNSSDecoderState(state; raw_data = GPSLNAVData())
    end

    # Upvote the match
    curr_score = old_data[matching_idx].vote
    new_vote = increment_voting(curr_score, max_vote)

    # Best score among entries with same IODC
    best_score = 0
    for entry in old_data
        if entry.data.IODC == state.raw_data.IODC
            best_score = max(best_score, entry.vote)
        end
    end

    if best_score > curr_score
        # Another entry leads: reject
        old_data[matching_idx] = VotedGPSLNAVData(new_vote, state.raw_data)
        return GNSSDecoderState(state; raw_data = GPSLNAVData())
    end

    # Best (or tied) score: publish
    if new_vote == max_vote && length(old_data) > 1
        # Max votes: keep only this entry
        replace_vote_candidates!(old_data, VotedGPSLNAVData(new_vote, state.raw_data))
    else
        old_data[matching_idx] = VotedGPSLNAVData(new_vote, state.raw_data)
    end

    promote_data(state)
end

function validate_data(state::GNSSDecoderState{<:GPSLNAVData})
    IODC = state.raw_data.IODC
    if has_positioning_messages(state.raw_data) &&
       IODC !== nothing &&
       IODC & 0xff == state.raw_data.IODE_subframe2 == state.raw_data.IODE_subframe3
        state = confirm_data(state)
    end
    return state
end

"""
$(TYPEDSIGNATURES)

Check if the GPS satellite is healthy and usable for positioning.

Healthy only if the 6-bit `sv_health` field of subframe 1 is zero (`"000000"`).

!!! warning

    Requires subframe 1 to be decoded; check that `state.data.sv_health` is not
    `nothing` before relying on this result.

# Arguments

  - `state::GNSSDecoderState{<:GPSLNAVData}`: GPS L1 decoder state with decoded data

# Returns

  - `Bool`: `true` if satellite health status indicates normal operation

# Example

```julia
state = GPSL1CADecoderState(1)
state = decode!(state, bits, num_bits)
if is_sat_healthy(state)
    # Safe to use for positioning
end
```

# See Also

  - [`GPSL1CADecoderState`](@ref): Create decoder state
  - [`decode!`](@ref): Decode navigation data
"""
function is_sat_healthy(state::GNSSDecoderState{<:GPSLNAVData})
    state.data.sv_health == 0
end
