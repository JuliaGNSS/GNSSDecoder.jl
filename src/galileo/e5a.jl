# UInt512 buffer for Galileo E5a F/NAV sync — holds the 512-symbol page window
# (500 syncro + 12 sync pattern) hard-sliced for the bit-pattern sync check.
BitIntegers.@define_integers 512

# Galileo E5a uses the shared K=7 NSC FEC (`GALILEO_VITERBI_POLY`, see
# `galileo.jl`) — the same code as E1B. After the 61×8 block deinterleave a F/NAV
# page carries 488 encoded symbols which the Viterbi decoder maps back to 244
# trellis steps: 238 information bits + 6 tail bits. AFF3CT's `ConvViterbiDecoder`
# is configured with K = 238, N = 488 and those polynomials; it applies the
# trellis termination internally and returns exactly the 238 information bits (the
# 6 tail bits are consumed by termination), which is the complete F/NAV page
# payload (page type + data + CRC) the parser expects. Cross-checked against
# GNSS-SDR (`Galileo_FNAV.h`) and PocketSDR (`decode_gal_FNAV` / `decode_gal_syms`,
# reshape 8×61, `syms[1::2] ^= 1`).
const GALILEO_E5A_VITERBI_K = 238
const GALILEO_E5A_VITERBI_N = 488

# Block deinterleaver dimensions for F/NAV: the ICD interleaver is 8 rows × 61
# columns. As with E1B (passed as `(30, 8)`), `deinterleave`'s first argument is
# the ICD column count and the second the ICD row count, so F/NAV uses `(61, 8)`.
const GALILEO_E5A_INTERLEAVER_ROWS = 61
const GALILEO_E5A_INTERLEAVER_COLS = 8

"""
    GalileoE5aConstants

GTRF constants and F/NAV message structure parameters for Galileo E5a signal decoding.

The physical constants are defined in the Galileo OS SIS ICD (Open Service Signal-In-Space
Interface Control Document) and are used for computing satellite positions and clock
corrections from broadcast ephemeris data.

# Fields

  - `syncro_sequence_length::Int`: Length of one F/NAV page in channel symbols (500 symbols = 10 s at 50 sps)
  - `preamble::UInt16`: F/NAV synchronisation pattern (101101110000 binary)
  - `preamble_length::Int`: Length of the sync pattern in symbols (12)
  - `PI::Float64`: Mathematical constant π = 3.1415926535898 (Galileo OS SIS ICD Table 68)
  - `Ω_dot_e::Float64`: Mean angular velocity of the Earth = 7.2921151467×10⁻⁵ rad/s
  - `c::Float64`: Speed of light = 2.99792458×10⁸ m/s
  - `μ::Float64`: Geocentric gravitational constant = 3.986004418×10¹⁴ m³/s²
  - `F::Float64`: Relativistic correction constant = -4.442807309×10⁻¹⁰ s/√m

# Reference

Galileo OS SIS ICD, Issue 2.2, §4.2 and Table 68
"""
Base.@kwdef struct GalileoE5aConstants <: AbstractGNSSConstants
    syncro_sequence_length::Int = 500
    preamble::UInt16 = 0b101101110000
    preamble_length::Int = 12
    PI::Float64 = GNSS_PI
    Ω_dot_e::Float64 = EARTH_ROTATION_RATE
    c::Float64 = SPEED_OF_LIGHT
    μ::Float64 = GALILEO_μ
    F::Float64 = GALILEO_F
end

"""
    GalileoE5aData

Decoded Galileo E5a F/NAV navigation message data.

Contains ephemeris, clock correction, signal health, group delay, ionospheric
correction, GST-UTC and GST-GPS conversion, and almanac parameters decoded from
the Galileo F/NAV message (the data component broadcast on E5a-I). All parameters
conform to the Galileo OS SIS ICD, Issue 2.2, §5.1.

Unlike I/NAV (E1B/E5b), F/NAV carries only the E5a signal-health (`E5a_HS`) and
data-validity (`E5a_DVS`) flags and a single broadcast group delay
(`BGD(E1, E5a)`); there is no Reduced CED and no E5b/E1-B field. Angular
quantities are stored in **radians** (the ICD broadcasts them in semi-circles;
the decoder multiplies by π), matching the convention used by [`GalileoE1BData`](@ref).

# Galileo System Time (GST) Fields

  - `WN::Int64`: Week Number (0-4095)
  - `TOW::Int64`: Time of Week at the start of the page (seconds, 0-604799)

# Satellite Identification (Word Type 1)

  - `SVID::Int`: Satellite Identifier (1-36 nominal range)

# Ephemeris Parameters (Word Types 2-4)

  - `t_0e::Float64`: Ephemeris reference time (seconds)
  - `M_0::Float64`: Mean anomaly at reference time (radians)
  - `e::Float64`: Eccentricity (dimensionless)
  - `sqrt_A::Float64`: Square root of semi-major axis (√m)
  - `Ω_0::Float64`: Longitude of ascending node at weekly epoch (radians)
  - `i_0::Float64`: Inclination angle at reference time (radians)
  - `ω::Float64`: Argument of perigee (radians)
  - `i_dot::Float64`: Rate of change of inclination angle (radians/s)
  - `Ω_dot::Float64`: Rate of change of right ascension (radians/s)
  - `Δn::Float64`: Mean motion difference from computed value (radians/s)
  - `C_uc::Float64`: Cosine harmonic correction to argument of latitude (rad)
  - `C_us::Float64`: Sine harmonic correction to argument of latitude (rad)
  - `C_rc::Float64`: Cosine harmonic correction to orbit radius (meters)
  - `C_rs::Float64`: Sine harmonic correction to orbit radius (meters)
  - `C_ic::Float64`: Cosine harmonic correction to inclination (rad)
  - `C_is::Float64`: Sine harmonic correction to inclination (rad)

# Signal-In-Space Accuracy (Word Type 1)

  - `SISA_e1_e5a::Int`: SISA index for dual frequency E1-E5a (Table 91/92; 255 = NAPA)

# Clock Correction Parameters (Word Type 1)

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
  - `num_bits_after_valid_syncro_sequence_after_last_TOW::Int`: Symbols since last TOW sync

# Signal Health and Data Validity (Word Type 1)

  - `signal_health_e5a::SignalHealth`: E5a signal health status (0=OK, 1=out of service, 2=Extended Operations Mode, 3=in test)
  - `data_validity_status_e5a::DataValidityStatus`: E5a data validity (0=valid, 1=working without guarantee)

# Broadcast Group Delay (Word Type 1)

  - `broadcast_group_delay_e1_e5a::Float64`: E1-E5a group delay correction (seconds)

# Ionospheric Correction (Word Type 1)

  - `a_i0::Float64`: Effective Ionisation Level 1st-order coefficient (sfu)
  - `a_i1::Float64`: Effective Ionisation Level 2nd-order coefficient (sfu/degree)
  - `a_i2::Float64`: Effective Ionisation Level 3rd-order coefficient (sfu/degree²)
  - `iono_storm_flag_region1..5::Bool`: Ionospheric Disturbance (storm) flags for regions 1-5

# GST-UTC Conversion (Word Type 4)

  - `A_0_utc::Float64`: Constant term of polynomial (s)
  - `A_1_utc::Float64`: 1st-order term of polynomial (s/s)
  - `Δt_LS::Int`: Leap Second count before leap second adjustment (s)
  - `t_0t::Int`: UTC data reference Time of Week (s)
  - `WN_0t::Int`: UTC data reference Week Number (8-bit, modulo 256)
  - `WN_LSF::Int`: Week Number of leap second adjustment (8-bit, modulo 256)
  - `DN::Int`: Day Number at end of which leap second becomes effective (1=Sunday … 7=Saturday)
  - `Δt_LSF::Int`: Leap Second count after leap second adjustment (s)

# GST-GPS Conversion / GGTO (Word Type 4)

  - `A_0G::Float64`: Constant term of GST-GPS offset polynomial (s)
  - `A_1G::Float64`: Rate of change of GST-GPS offset (s/s)
  - `t_0G::Int`: GGTO reference time (s)
  - `WN_0G::Int`: GGTO reference Week Number (6-bit)

# Almanac (Word Types 5-6)

  - `almanacs::Dictionary{Int,GalileoAlmanac}`: Decoded almanacs keyed by SVID.
    Galileo broadcasts three almanacs across the word-type-5/6 pair: SVID-1 (full
    in WT5), SVID-2 (split across WT5 and WT6), and SVID-3 (full in WT6). The
    in-flight SVID-2 partial lives in the decoder cache and is flushed here only
    once WT6 completes it with a consistent `IOD_a`. SVID-3, though fully carried
    in WT6, inherits its reference epoch (`WN_a`/`t_0a`) from the paired WT5; when
    that partial is missing (mid-stream acquisition or an IOD cutover), the SVID-3
    record is still stored with `WN_a`/`t_0a` left `nothing` — the decoder keeps
    whatever it can decode rather than discarding it. Because the epoch is shared
    by every almanac of a given `IOD_a`, a later WT5 back-fills it into any such
    partial record, so a one-shot WT6 orbit becomes usable even if that WT6 never
    reappears. **An almanac may therefore be incomplete**: any field can be
    `nothing`, and in particular a record's reference epoch (`WN_a`/`t_0a`) may be
    absent until a matching WT5 arrives — check the fields you need before using a
    record. F/NAV almanacs carry the E5a health (`signal_health_e5a`); the E5b/E1-B
    almanac-health fields are left `nothing`.

# Reference

Galileo OS SIS ICD, Issue 2.2, §5.1, Tables 75-80
"""
Base.@kwdef struct GalileoE5aData <: AbstractGalileoData
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

    SISA_e1_e5a::Union{Nothing,Int} = nothing

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

    signal_health_e5a::Union{Nothing,SignalHealth} = nothing
    data_validity_status_e5a::Union{Nothing,DataValidityStatus} = nothing

    broadcast_group_delay_e1_e5a::Union{Nothing,Float64} = nothing

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
end

function GalileoE5aData(
    data::GalileoE5aData;
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
    SISA_e1_e5a = data.SISA_e1_e5a,
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
    signal_health_e5a = data.signal_health_e5a,
    data_validity_status_e5a = data.data_validity_status_e5a,
    broadcast_group_delay_e1_e5a = data.broadcast_group_delay_e1_e5a,
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
)
    GalileoE5aData(
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
        SISA_e1_e5a,
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
        signal_health_e5a,
        data_validity_status_e5a,
        broadcast_group_delay_e1_e5a,
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
    )
end

# As with GalileoE1BData, the mutable `almanacs::Dictionary` field makes the
# default struct `==` (which falls back to `===`) too strict. Compare field-by-field.
function Base.:(==)(a::GalileoE5aData, b::GalileoE5aData)
    for f in fieldnames(GalileoE5aData)
        getfield(a, f) == getfield(b, f) || return false
    end
    return true
end

# `is_ephemeris_decoded` and `is_clock_correction_decoded` are per-constellation
# facts (identical fields for I/NAV and F/NAV), defined once on
# `AbstractGalileoData` in `galileo/galileo.jl`. Only the health-status check
# below is genuinely per-signal (E5a carries only E5a health).
function is_health_status_decoded(data::GalileoE5aData)
    !isnothing(data.signal_health_e5a) && !isnothing(data.data_validity_status_e5a)
end

function is_decoding_completed_for_positioning(data::GalileoE5aData)
    !isnothing(data.TOW) &&
        !isnothing(data.WN) &&
        !isnothing(data.broadcast_group_delay_e1_e5a) &&
        is_ephemeris_decoded(data) &&
        is_clock_correction_decoded(data) &&
        is_health_status_decoded(data)
end

"""
$(TYPEDEF)

Per-decoder cache for Galileo E5a F/NAV.

Holds the soft-symbol `CircularDeque{Float32}` (capacity = 500 + 12 = 512), the
in-flight almanac partial used to stitch SVID-2's halves across word types 5 and
6 (plus its 4-bit `Ω_0` MSB, which WT6 completes with a 12-bit LSB), and the
long-lived AFF3CT Viterbi decoder. Unlike I/NAV there is no even/odd page
stitching: each F/NAV page is a complete, independently CRC-protected word.

The decoder consumes *soft symbols* end-to-end: the sync hook hard-slices the
deque tail only for the 12-symbol sync-pattern match, while the K=7 NSC FEC is
undone on the raw `Float32` LLRs via AFF3CT.jl's `ConvViterbiDecoder`. The
deque-backed input boundary is identical to L1 C/A and E1B so the public API is
uniform.

# Fields

$(TYPEDFIELDS)
"""
struct GalileoE5aCache <: AbstractGNSSCache
    """
    Soft-symbol buffer (512 = 500 syncro + 12 sync pattern)
    """
    soft_buffer::CircularDeque{Float32}
    """
    SVID-2 almanac partial decoded from word type 5, completed by word type 6.
    """
    almanac_chain_partial::GalileoAlmanac
    """
    SVID-2 `Ω_0` most-significant 4 bits from word type 5 (LSB-aligned), combined
    with the 12 LSBs in word type 6. `nothing` when no WT5 partial is in flight.
    """
    almanac_chain_omega0_msb::Union{Nothing,Int}
    """
    AFF3CT K=7 NSC Viterbi decoder, built once and reused across pages.
    """
    viterbi_decoder::Aff3ct.ConvViterbiDecoder
    """
    488 encoded symbols copied out of `soft_buffer` (with the 180° polarity
    applied) for each synced page, reused rather than reallocated per page.
    """
    soft_page::Vector{Float32}
    """
    238 decoded bits unpacked for the CRC-24Q check, reused across pages.
    """
    crc_bits::Vector{Bool}
end

GalileoE5aCache() = GalileoE5aCache(
    CircularDeque{Float32}(512),
    GalileoAlmanac(),
    nothing,
    Aff3ct.ConvViterbiDecoder(
        GALILEO_E5A_VITERBI_K,
        GALILEO_E5A_VITERBI_N,
        GALILEO_VITERBI_POLY,
    ),
    Vector{Float32}(undef, GALILEO_E5A_VITERBI_N),
    Vector{Bool}(undef, GALILEO_E5A_VITERBI_K),
)

function GalileoE5aCache(
    cache::GalileoE5aCache;
    soft_buffer = cache.soft_buffer,
    almanac_chain_partial = cache.almanac_chain_partial,
    almanac_chain_omega0_msb = cache.almanac_chain_omega0_msb,
    viterbi_decoder = cache.viterbi_decoder,
    soft_page = cache.soft_page,
    crc_bits = cache.crc_bits,
)
    GalileoE5aCache(
        soft_buffer,
        almanac_chain_partial,
        almanac_chain_omega0_msb,
        viterbi_decoder,
        soft_page,
        crc_bits,
    )
end

function Base.:(==)(a::GalileoE5aCache, b::GalileoE5aCache)
    deques_equal(a.soft_buffer, b.soft_buffer) &&
        a.almanac_chain_partial == b.almanac_chain_partial &&
        a.almanac_chain_omega0_msb == b.almanac_chain_omega0_msb
end

"""
$(TYPEDSIGNATURES)

Create a decoder state for Galileo E5a F/NAV navigation messages.

Initializes a [`GNSSDecoderState`](@ref) configured for decoding Galileo E5a
(Open Service) F/NAV navigation messages. The decoder extracts ephemeris, clock
correction, ionospheric parameters, almanac, and health data from the 50 sps
F/NAV data stream broadcast on the E5a-I component using Viterbi decoding.

# Arguments

  - `prn::Int`: Pseudo-Random Noise code identifier (1-36 for Galileo satellites)

# Returns

  - `GNSSDecoderState{GalileoE5aData}`: Initialized decoder state for Galileo E5a

# Example

```julia
state = GalileoE5aDecoderState(21)  # Create decoder for PRN 21
state = decode(state, soft_symbols, num_symbols)
if is_sat_healthy(state)
    # Use state.data for positioning
end
```

# See Also

  - [`GNSSDecoderState`](@ref): The underlying state structure
  - [`decode`](@ref): Decode soft symbols using this state
  - [`reset_decoder_state`](@ref): Reset after signal loss
  - [`is_sat_healthy`](@ref): Check satellite health status
"""
function GalileoE5aDecoderState(prn)
    GNSSDecoderState(
        prn,
        GalileoE5aData(),
        GalileoE5aData(),
        GalileoE5aConstants(),
        GalileoE5aCache(),
        nothing,
        false,
    )
end

# Dispatch from a GNSSSignals system type, mirroring `GNSSDecoderState(::GalileoE1B, …)`.
# F/NAV is broadcast on the E5a-I (data) component — `GalileoE5aI` — while `GalileoE5aQ`
# is the dataless pilot, so only E5a-I maps to this decoder.
function GNSSDecoderState(system::GalileoE5aI, prn)
    GalileoE5aDecoderState(prn)
end

"""
$(TYPEDSIGNATURES)

Reset the Galileo E5a decoder state after a signal loss or reacquisition.

Clears the soft-symbol buffer and the time-of-week (TOW) field while preserving
other decoded ephemeris and clock data in `raw_data`, mirroring the E1B reset
semantics. The week number (`WN`) is intentionally not reset (it is broadcast
less frequently than TOW).
"""
function reset_decoder_state(state::GNSSDecoderState{<:GalileoE5aData})
    empty!(state.cache.soft_buffer)
    GNSSDecoderState(
        state;
        raw_data = GalileoE5aData(
            state.raw_data;
            TOW = nothing,
            num_bits_after_valid_syncro_sequence_after_last_TOW = nothing,
        ),
        data = GalileoE5aData(),
        num_bits_after_valid_syncro_sequence = nothing,
    )
end

packed_buffer_type(::GNSSDecoderState{<:GalileoE5aData}) = UInt512

"""
    galileo_e5a_viterbi(decoder, soft_page) -> UInt256

Recover one F/NAV page's 238-bit payload from its `soft_page` — the 488
polarity-corrected `Float32` LLR soft symbols of a Galileo E5a page (the
`syncro_sequence_length - preamble_length` encoded symbols between the leading
and trailing 12-symbol sync patterns). `decoder` is the cache's long-lived
`Aff3ct.ConvViterbiDecoder`, reused across pages.

The transmit FEC chain (Galileo OS SIS ICD, Issue 2.2, §4.1.4 / §4.2.5) is undone
in order on the soft symbols:

 1. **61×8 block deinterleave** of the 488 LLRs (`deinterleave` from `src/deinterleave.jl`).
 2. **Invert every second symbol** — the spec inverts the G2 output of the rate-1/2
    encoder. On soft symbols an inversion is a sign flip (negation), so confidence
    magnitudes are preserved.
 3. **K=7 NSC Viterbi** via AFF3CT.jl's `ConvViterbiDecoder`. AFF3CT's LLR sign
    convention matches ours (positive ⇒ bit 0), so the LLRs feed in directly. The
    decoder returns the 238 information bits (the 6 tail bits are consumed by
    trellis termination).

The 238 decoded bits are packed MSB-first into the low bits of a `UInt256`, ready
for `get_bits`/`get_twos_complement_num` field extraction.

Thin wrapper over the shared [`galileo_viterbi`](@ref) with E5a's 61×8 interleaver
shape and `UInt256` payload type.
"""
galileo_e5a_viterbi(
    decoder::Aff3ct.ConvViterbiDecoder,
    soft_page::AbstractVector{Float32},
) = galileo_viterbi(
    decoder,
    soft_page,
    GALILEO_E5A_INTERLEAVER_ROWS,
    GALILEO_E5A_INTERLEAVER_COLS,
    UInt256,
)

# Combine SVID-2's split right-ascension: WT5 carries the 4 MSBs, WT6 the 12 LSBs,
# of a 16-bit two's-complement value scaled by π·2⁻¹⁵ (semicircles → radians).
function combine_almanac_omega0(msb::Int, lsb::Int, PI::Float64)
    raw = (UInt16(msb) << 12) | UInt16(lsb)
    return get_twos_complement_num(raw, 16, 1, 16) * PI / (1 << 15)
end

# IOD-keyed epoch back-patch. WN_a/t_0a are shared by every almanac of a given
# IOD_a, so a freshly decoded WT5 (which carries them) can complete any earlier
# record left without its reference epoch — most notably an SVID-3 decoded from a
# WT6 whose paired WT5 was missed (mid-stream acquisition / IOD cutover). This
# lets a one-shot WT6 orbit become usable once a later WT5 arrives, even if that
# WT6 never reappears; without it the WT6's orbital block would be stranded.
# Returns `almanacs` unchanged when nothing matches, else a patched copy (the
# input dictionary, shared with `raw_data`, is never mutated in place).
function backpatch_almanac_epochs(
    almanacs::Union{Nothing,Dictionary{Int,GalileoAlmanac}},
    IOD_a::Int,
    WN_a::Int,
    t_0a::Int,
)
    isnothing(almanacs) && return almanacs
    patched = almanacs
    for SVID in keys(almanacs)
        alm = almanacs[SVID]
        if alm.IOD_a == IOD_a && (isnothing(alm.WN_a) || isnothing(alm.t_0a))
            patched === almanacs && (patched = copy(almanacs))
            set!(patched, SVID, GalileoAlmanac(alm; WN_a, t_0a))
        end
    end
    return patched
end

function decode_syncro_sequence(state::GNSSDecoderState{<:GalileoE5aData}, buffer)
    # The 488 encoded symbols sit between the leading 12-symbol sync pattern and
    # the trailing sync pattern of the next page (deque indices
    # preamble_length+1 .. syncro_sequence_length). Resolve the 180-degree
    # polarity ambiguity by negating the LLRs when the sync hook flagged the page
    # as inverted.
    deque = soft_buffer(state)
    sign = state.is_shifted_by_180_degrees ? -1.0f0 : 1.0f0
    soft_page = state.cache.soft_page
    @inbounds for i = 1:GALILEO_E5A_VITERBI_N
        soft_page[i] = sign * deque[state.constants.preamble_length+i]
    end
    bits = galileo_e5a_viterbi(state.cache.viterbi_decoder, soft_page)

    state = GNSSDecoderState(
        state;
        raw_data = GalileoE5aData(
            state.raw_data;
            num_pages_after_last_TOW = state.raw_data.num_pages_after_last_TOW + 1,
        ),
    )

    # F/NAV CRC-24Q is computed over the 214-bit (page type + data) prefix and
    # appended as bits 215-238; a clean page satisfies crc24q(all 238 bits) == 0.
    crc_bits = state.cache.crc_bits
    @inbounds for i = 1:GALILEO_E5A_VITERBI_K
        crc_bits[i] = get_bit(bits, GALILEO_E5A_VITERBI_K, i)
    end
    crc24q(crc_bits) == 0 || return state

    PI = state.constants.PI
    page_type = get_bits(bits, GALILEO_E5A_VITERBI_K, 1, 6)

    if page_type == 1
        SVID = Int(get_bits(bits, 238, 7, 6))
        IOD_nav1 = get_bits(bits, 238, 13, 10)
        t_0c = get_bits(bits, 238, 23, 14) * 60
        a_f0 = get_twos_complement_num(bits, 238, 37, 31) / (1 << 34)
        a_f1 = get_twos_complement_num(bits, 238, 68, 21) / (1 << 46)
        a_f2 = get_twos_complement_num(bits, 238, 89, 6) / Float64(1 << 59)
        SISA_e1_e5a = Int(get_bits(bits, 238, 95, 8))
        a_i0 = get_bits(bits, 238, 103, 11) / (1 << 2)
        a_i1 = get_twos_complement_num(bits, 238, 114, 11) / (1 << 8)
        a_i2 = get_twos_complement_num(bits, 238, 125, 14) / (1 << 15)
        iono_storm_flag_region1 = get_bit(bits, 238, 139)
        iono_storm_flag_region2 = get_bit(bits, 238, 140)
        iono_storm_flag_region3 = get_bit(bits, 238, 141)
        iono_storm_flag_region4 = get_bit(bits, 238, 142)
        iono_storm_flag_region5 = get_bit(bits, 238, 143)
        broadcast_group_delay_e1_e5a =
            get_twos_complement_num(bits, 238, 144, 10) / Float64(1 << 32)
        signal_health_e5a = SignalHealth(get_bits(bits, 238, 154, 2))
        WN = get_bits(bits, 238, 156, 12)
        TOW = get_bits(bits, 238, 168, 20)
        data_validity_status_e5a = DataValidityStatus(get_bit(bits, 238, 188))
        state = GNSSDecoderState(
            state;
            raw_data = GalileoE5aData(
                state.raw_data;
                SVID,
                IOD_nav1,
                t_0c,
                a_f0,
                a_f1,
                a_f2,
                SISA_e1_e5a,
                a_i0,
                a_i1,
                a_i2,
                iono_storm_flag_region1,
                iono_storm_flag_region2,
                iono_storm_flag_region3,
                iono_storm_flag_region4,
                iono_storm_flag_region5,
                broadcast_group_delay_e1_e5a,
                signal_health_e5a,
                data_validity_status_e5a,
                WN,
                TOW,
                num_pages_after_last_TOW = 1,
                num_bits_after_valid_syncro_sequence_after_last_TOW = state.num_bits_after_valid_syncro_sequence,
            ),
        )
    elseif page_type == 2
        IOD_nav2 = get_bits(bits, 238, 7, 10)
        M_0 = get_twos_complement_num(bits, 238, 17, 32) * PI / (1 << 31)
        Ω_dot = get_twos_complement_num(bits, 238, 49, 24) * PI / Float64(1 << 43)
        e = get_bits(bits, 238, 73, 32) / Float64(1 << 33)
        sqrt_A = get_bits(bits, 238, 105, 32) / (1 << 19)
        Ω_0 = get_twos_complement_num(bits, 238, 137, 32) * PI / (1 << 31)
        i_dot = get_twos_complement_num(bits, 238, 169, 14) * PI / Float64(1 << 43)
        WN = get_bits(bits, 238, 183, 12)
        TOW = get_bits(bits, 238, 195, 20)
        state = GNSSDecoderState(
            state;
            raw_data = GalileoE5aData(
                state.raw_data;
                IOD_nav2,
                M_0,
                Ω_dot,
                e,
                sqrt_A,
                Ω_0,
                i_dot,
                WN,
                TOW,
                num_pages_after_last_TOW = 1,
                num_bits_after_valid_syncro_sequence_after_last_TOW = state.num_bits_after_valid_syncro_sequence,
            ),
        )
    elseif page_type == 3
        IOD_nav3 = get_bits(bits, 238, 7, 10)
        i_0 = get_twos_complement_num(bits, 238, 17, 32) * PI / (1 << 31)
        ω = get_twos_complement_num(bits, 238, 49, 32) * PI / (1 << 31)
        Δn = get_twos_complement_num(bits, 238, 81, 16) * PI / Float64(1 << 43)
        C_uc = get_twos_complement_num(bits, 238, 97, 16) / Float64(1 << 29)
        C_us = get_twos_complement_num(bits, 238, 113, 16) / Float64(1 << 29)
        C_rc = get_twos_complement_num(bits, 238, 129, 16) / (1 << 5)
        C_rs = get_twos_complement_num(bits, 238, 145, 16) / (1 << 5)
        t_0e = get_bits(bits, 238, 161, 14) * 60
        WN = get_bits(bits, 238, 175, 12)
        TOW = get_bits(bits, 238, 187, 20)
        state = GNSSDecoderState(
            state;
            raw_data = GalileoE5aData(
                state.raw_data;
                IOD_nav3,
                i_0,
                ω,
                Δn,
                C_uc,
                C_us,
                C_rc,
                C_rs,
                t_0e,
                WN,
                TOW,
                num_pages_after_last_TOW = 1,
                num_bits_after_valid_syncro_sequence_after_last_TOW = state.num_bits_after_valid_syncro_sequence,
            ),
        )
    elseif page_type == 4
        IOD_nav4 = get_bits(bits, 238, 7, 10)
        C_ic = get_twos_complement_num(bits, 238, 17, 16) / Float64(1 << 29)
        C_is = get_twos_complement_num(bits, 238, 33, 16) / Float64(1 << 29)
        A_0_utc = get_twos_complement_num(bits, 238, 49, 32) / Float64(1 << 30)
        A_1_utc = get_twos_complement_num(bits, 238, 81, 24) / Float64(1 << 50)
        Δt_LS = Int(get_twos_complement_num(bits, 238, 105, 8))
        t_0t = Int(get_bits(bits, 238, 113, 8) * 3600)
        WN_0t = Int(get_bits(bits, 238, 121, 8))
        WN_LSF = Int(get_bits(bits, 238, 129, 8))
        DN = Int(get_bits(bits, 238, 137, 3))
        Δt_LSF = Int(get_twos_complement_num(bits, 238, 140, 8))
        t_0G = Int(get_bits(bits, 238, 148, 8) * 3600)
        A_0G = get_twos_complement_num(bits, 238, 156, 16) / Float64(1 << 35)
        A_1G = get_twos_complement_num(bits, 238, 172, 12) / Float64(1 << 51)
        WN_0G = Int(get_bits(bits, 238, 184, 6))
        TOW = get_bits(bits, 238, 190, 20)
        state = GNSSDecoderState(
            state;
            raw_data = GalileoE5aData(
                state.raw_data;
                IOD_nav4,
                C_ic,
                C_is,
                A_0_utc,
                A_1_utc,
                Δt_LS,
                t_0t,
                WN_0t,
                WN_LSF,
                DN,
                Δt_LSF,
                t_0G,
                A_0G,
                A_1G,
                WN_0G,
                TOW,
                num_pages_after_last_TOW = 1,
                num_bits_after_valid_syncro_sequence_after_last_TOW = state.num_bits_after_valid_syncro_sequence,
            ),
        )
    elseif page_type == 5
        IOD_a = Int(get_bits(bits, 238, 7, 4))
        WN_a = Int(get_bits(bits, 238, 11, 2))
        t_0a = Int(get_bits(bits, 238, 13, 10) * 600)
        # SVID-1: fully contained in word type 5 → flush immediately.
        SVID1 = Int(get_bits(bits, 238, 23, 6))
        almanac1 = GalileoAlmanac(;
            SVID = SVID1,
            Δsqrt_A = get_twos_complement_num(bits, 238, 29, 13) / (1 << 9),
            e = get_bits(bits, 238, 42, 11) / (1 << 16),
            ω = get_twos_complement_num(bits, 238, 53, 16) * PI / (1 << 15),
            δi = get_twos_complement_num(bits, 238, 69, 11) * PI / (1 << 14),
            Ω_0 = get_twos_complement_num(bits, 238, 80, 16) * PI / (1 << 15),
            Ω_dot = get_twos_complement_num(bits, 238, 96, 11) * PI / Float64(1 << 33),
            M_0 = get_twos_complement_num(bits, 238, 107, 16) * PI / (1 << 15),
            a_f0 = get_twos_complement_num(bits, 238, 123, 16) / Float64(1 << 19),
            a_f1 = get_twos_complement_num(bits, 238, 139, 13) / Float64(1 << 38),
            signal_health_e5a = SignalHealth(get_bits(bits, 238, 152, 2)),
            IOD_a,
            WN_a,
            t_0a,
        )
        # SVID-2: first half (orbital shape + Ω_0 MSB) in word type 5; the
        # remainder arrives in word type 6.
        SVID2 = Int(get_bits(bits, 238, 154, 6))
        almanac2_partial = GalileoAlmanac(;
            SVID = SVID2,
            Δsqrt_A = get_twos_complement_num(bits, 238, 160, 13) / (1 << 9),
            e = get_bits(bits, 238, 173, 11) / (1 << 16),
            ω = get_twos_complement_num(bits, 238, 184, 16) * PI / (1 << 15),
            δi = get_twos_complement_num(bits, 238, 200, 11) * PI / (1 << 14),
            IOD_a,
            WN_a,
            t_0a,
        )
        omega0_msb = Int(get_bits(bits, 238, 211, 4))

        # Complete any earlier record still missing its reference epoch (e.g. an
        # SVID-3 from a WT6 whose paired WT5 was missed) — WN_a/t_0a are shared by
        # all almanacs of this IOD_a.
        almanacs = backpatch_almanac_epochs(state.raw_data.almanacs, IOD_a, WN_a, t_0a)
        if SVID1 >= 1
            almanacs =
                isnothing(almanacs) ? Dictionary{Int,GalileoAlmanac}() : copy(almanacs)
            set!(almanacs, SVID1, almanac1)
        end
        valid_SVID2 = SVID2 >= 1
        state = GNSSDecoderState(
            state;
            raw_data = GalileoE5aData(state.raw_data; almanacs),
            cache = GalileoE5aCache(
                state.cache;
                almanac_chain_partial = valid_SVID2 ? almanac2_partial : GalileoAlmanac(),
                almanac_chain_omega0_msb = valid_SVID2 ? omega0_msb : nothing,
            ),
        )
    elseif page_type == 6
        IOD_a = Int(get_bits(bits, 238, 7, 4))
        # SVID-2 completion: combine the WT5 Ω_0 MSBs with the WT6 LSBs and add
        # the remaining orbital/clock/health terms. Only flush if the WT5 partial
        # is intact and its IOD_a matches.
        omega0_lsb = Int(get_bits(bits, 238, 11, 12))
        partial = state.cache.almanac_chain_partial
        msb = state.cache.almanac_chain_omega0_msb
        almanacs = state.raw_data.almanacs
        if !isnothing(msb) && !isnothing(partial.SVID) && partial.IOD_a == IOD_a
            completed2 = GalileoAlmanac(
                partial;
                Ω_0 = combine_almanac_omega0(msb, omega0_lsb, PI),
                Ω_dot = get_twos_complement_num(bits, 238, 23, 11) * PI / Float64(1 << 33),
                M_0 = get_twos_complement_num(bits, 238, 34, 16) * PI / (1 << 15),
                a_f0 = get_twos_complement_num(bits, 238, 50, 16) / Float64(1 << 19),
                a_f1 = get_twos_complement_num(bits, 238, 66, 13) / Float64(1 << 38),
                signal_health_e5a = SignalHealth(get_bits(bits, 238, 79, 2)),
            )
            almanacs =
                isnothing(almanacs) ? Dictionary{Int,GalileoAlmanac}() : copy(almanacs)
            set!(almanacs, completed2.SVID, completed2)
        end
        # SVID-3: orbital/clock/health are fully contained in word type 6, but its
        # almanac reference epoch (WN_a, t_0a) is broadcast only in the paired word
        # type 5. Inherit them from the cached WT5 partial when its IOD_a matches;
        # otherwise the SVID-3 record is still stored — with WN_a/t_0a left
        # `nothing` — so nothing decodable is discarded. This happens on mid-stream
        # acquisition or an IOD cutover that lands WT6 without its paired WT5; the
        # epoch is filled in later by `backpatch_almanac_epochs` when the matching
        # WT5 arrives (or wholesale by the next full WT5→WT6 cycle).
        SVID3 = Int(get_bits(bits, 238, 81, 6))
        if SVID3 >= 1
            shared_wn_a, shared_t_0a =
                (!isnothing(partial.SVID) && partial.IOD_a == IOD_a) ?
                (partial.WN_a, partial.t_0a) : (nothing, nothing)
            almanac3 = GalileoAlmanac(;
                SVID = SVID3,
                Δsqrt_A = get_twos_complement_num(bits, 238, 87, 13) / (1 << 9),
                e = get_bits(bits, 238, 100, 11) / (1 << 16),
                ω = get_twos_complement_num(bits, 238, 111, 16) * PI / (1 << 15),
                δi = get_twos_complement_num(bits, 238, 127, 11) * PI / (1 << 14),
                Ω_0 = get_twos_complement_num(bits, 238, 138, 16) * PI / (1 << 15),
                Ω_dot = get_twos_complement_num(bits, 238, 154, 11) * PI / Float64(1 << 33),
                M_0 = get_twos_complement_num(bits, 238, 165, 16) * PI / (1 << 15),
                a_f0 = get_twos_complement_num(bits, 238, 181, 16) / Float64(1 << 19),
                a_f1 = get_twos_complement_num(bits, 238, 197, 13) / Float64(1 << 38),
                signal_health_e5a = SignalHealth(get_bits(bits, 238, 210, 2)),
                IOD_a,
                WN_a = shared_wn_a,
                t_0a = shared_t_0a,
            )
            almanacs =
                isnothing(almanacs) ? Dictionary{Int,GalileoAlmanac}() : copy(almanacs)
            set!(almanacs, SVID3, almanac3)
        end
        state = GNSSDecoderState(
            state;
            raw_data = GalileoE5aData(state.raw_data; almanacs),
            cache = GalileoE5aCache(
                state.cache;
                almanac_chain_partial = GalileoAlmanac(),
                almanac_chain_omega0_msb = nothing,
            ),
        )
    end
    return state
end

function validate_data(state::GNSSDecoderState{<:GalileoE5aData})
    if is_decoding_completed_for_positioning(state.raw_data) &&
       state.raw_data.IOD_nav1 ==
       state.raw_data.IOD_nav2 ==
       state.raw_data.IOD_nav3 ==
       state.raw_data.IOD_nav4
        num_bits_after_valid_syncro_sequence = 0
        if state.data.TOW == state.raw_data.TOW
            num_bits_after_valid_syncro_sequence =
                state.num_bits_after_valid_syncro_sequence
        elseif !isnothing(
            state.raw_data.num_bits_after_valid_syncro_sequence_after_last_TOW,
        )
            # Re-reference the symbol counter to the page that carried the most
            # recent TOW. A F/NAV word is exactly one page (syncro_sequence_length
            # symbols + the leading sync pattern).
            num_bits_after_valid_syncro_sequence =
                state.num_bits_after_valid_syncro_sequence - (
                    state.raw_data.num_bits_after_valid_syncro_sequence_after_last_TOW -
                    state.constants.syncro_sequence_length -
                    state.constants.preamble_length
                )
        else # first successful decoding
            num_bits_after_valid_syncro_sequence =
                state.constants.preamble_length +
                state.raw_data.num_pages_after_last_TOW *
                state.constants.syncro_sequence_length
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

Check if the Galileo satellite is healthy and usable for positioning on E5a.

Examines both the E5a signal-health status (`signal_health_e5a`) and the E5a
data-validity status (`data_validity_status_e5a`) from word type 1. A satellite
is considered healthy only if the signal health is `signal_ok` and the data
validity is `navigation_data_valid`.

!!! warning

    This requires that word type 1 has been successfully decoded. Check that
    `state.data.signal_health_e5a` is not `nothing` before relying on the result.
"""
function is_sat_healthy(state::GNSSDecoderState{<:GalileoE5aData})
    state.data.signal_health_e5a == signal_ok &&
        state.data.data_validity_status_e5a == navigation_data_valid
end
