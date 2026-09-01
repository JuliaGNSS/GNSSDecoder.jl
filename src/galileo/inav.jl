# Shared Galileo I/NAV navigation message core — Galileo OS SIS ICD, Issue 2.2
# §4.3 (I/NAV Message Description).
#
# I/NAV is broadcast on *two* signal components, E1-B and E5b-I, and the ICD is
# explicit that they "use the same page layout since the service provided on
# these frequencies is a dual frequency service, using frequency diversity.
# Only page sequencing is different" (§4.3.1). Both run at 250 sps with
# 250-symbol (1 s) page parts, the same 10-bit page-sync pattern
# `0101100000`, the same K=7 rate-1/2 NSC FEC over the same 30×8 block
# interleaver, the same even/odd page pairing into 128-bit nominal words, and
# the same CRC-24Q over the 196 protected bits of the page pair.
#
# The two page-part layouts differ only in fields this decoder does not read
# (ICD Table 38): the E1-B odd part spends 64 bits on OSNMA(40) + SAR(22) +
# Spare(2) where E5b-I has a single 64-bit "Reserved 1" field, and its trailing
# 8 bits are the SSP where E5b-I has "Reserved 2". Crucially the *protected*
# prefix is 82 bits either way (Even/Odd + Page Type + Data(2/2) + those 64
# bits), the CRC sits at odd-part bits 83-106 either way, and neither trailing
# 8-bit field is covered by the CRC. So the CRC scope, the word extraction, and
# every word-type parser below are bit-for-bit identical on the two signals.
#
# Consequently this file holds the whole I/NAV decoder — framing, FEC, sync,
# CRC gate, and all word-type parsing, plus the shared `GalileoINAVData`
# container — and `e1b.jl` / `e5b.jl` are thin signal layers carrying only
# their own `*Constants` alias, decoder-state constructor, `get_signal_type`
# mapping, and health-facet selection in `is_sat_healthy` (E1-B/C health and
# data validity vs E5b's — both fields come from the same word type 5, so both
# are always decoded together). This mirrors how GPS L5I and L2C share
# `gps/cnav.jl` and BeiDou B1I and B3I share `beidou/dnav.jl`.
#
# The one substantive sequencing consequence: word types 16 (Reduced CED) and
# 17-20 (FEC2 RS CED) are broadcast on E1-B only (Table 40), so an E5b-I
# decoder simply never observes them. Nothing in the parser needs to know that.

# The packed word for a Galileo I/NAV page pair — the even part's 114 bits
# followed by the odd part's leading 106, the 220 bits the CRC-24Q gate runs
# over — is `UInt288`, defined in `gnss.jl` because two other decoders reach
# for the same width.

# Page-part geometry, as module-level constants rather than only
# `GalileoINAVConstants` fields: `try_sync` below hands them to
# `find_preamble_in_deque` once per symbol, where being compile-time constants is
# what lets the bit packing unroll (see that function's note).
const GALILEO_INAV_PAGE_PART_SYMBOLS = 250
const GALILEO_INAV_SYNC_SYMBOLS = 10

# Galileo I/NAV uses the shared K=7 NSC FEC (`GALILEO_VITERBI_POLY`, see
# `galileo.jl`). After the 30×8 block deinterleave a page carries 240 encoded
# symbols which the Viterbi decoder maps back to 120 trellis steps: 114
# information bits + 6 tail bits. AFF3CT's `ConvViterbiDecoder` is configured with
# K = 114, N = 240 and those polynomials; it applies the trellis termination
# internally and returns exactly the 114 information bits (the 6 tail bits are
# consumed by termination and never surface), which is precisely the per-page
# payload the I/NAV parser expects. The interleaver shape below is passed to the
# shared `galileo_viterbi` helper at decode time.
const GALILEO_INAV_VITERBI_K = 114
const GALILEO_INAV_VITERBI_N = 240

# Octets spanned by the 220 bits the CRC-24Q gate runs over — the 196 protected
# bits of a page pair followed by the 24-bit checksum itself, i.e. the even
# part's 114 bits plus the odd part's leading 106 — for the packed-word
# `crc24q`. 220 rounds up to 28 octets; the four leading bits are zero and so are
# neutral, CRC-24Q initialising its register to zero.
const GALILEO_INAV_CRC_OCTETS = cld(114 + 106, 8)

# Columns of the I/NAV block interleaver; the 8 rows every Galileo channel shares
# are `GALILEO_INTERLEAVER_ROWS` in `galileo.jl` (ICD 8×30).
const GALILEO_INAV_INTERLEAVER_COLUMNS = 30

"""
    GalileoINAVConstants{S}

GTRF constants and I/NAV message structure parameters shared by the two Galileo
signal components that broadcast I/NAV, E1-B and E5b-I. Aliased per signal as
[`GalileoE1BConstants`](@ref GNSSDecoder.GalileoE1BConstants) (`S = :GalileoE1B`)
and [`GalileoE5bConstants`](@ref GNSSDecoder.GalileoE5bConstants)
(`S = :GalileoE5bI`); the field values are identical and the tag `S` only selects
the signal reported by [`get_signal_type`](@ref) and the health facet checked by
[`is_sat_healthy`](@ref).

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
Base.@kwdef struct GalileoINAVConstants{S} <: AbstractGNSSConstants
    syncro_sequence_length::Int = GALILEO_INAV_PAGE_PART_SYMBOLS
    preamble::UInt16 = 0b0101100000
    preamble_length::Int = GALILEO_INAV_SYNC_SYMBOLS
    PI::Float64 = GNSS_PI
    Ω_dot_e::Float64 = EARTH_ROTATION_RATE
    c::Float64 = SPEED_OF_LIGHT
    μ::Float64 = GALILEO_μ
    F::Float64 = GALILEO_F
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
  - `Δi_0_red::Float64`: Inclination delta from nominal (rad)
  - `Ω_0_red::Float64`: Longitude of ascending node at weekly epoch (rad)
  - `λ_0_red::Float64`: Mean argument of latitude, M0+ω (rad)
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
"""
$(TYPEDEF)

Per-decoder cache for Galileo I/NAV (E1-B and E5b-I).

Holds the soft-symbol `CircularDeque{Float32}` (capacity = 250 + 10 = 260), the
scratch buffer the polarity-resolved page part is copied into, the cached
even-page bits used to stitch two consecutive pages into a word, and the
two-position almanac-chain state needed to merge word types 7-10. Soft-symbol
buffering is shared across all signals; the rest is Galileo-specific.

The Galileo decoder consumes *soft symbols* end-to-end: the page-sync hook
hard-slices only the two 10-bit preamble windows out of the deque, while the K=7
NSC FEC is undone on the raw `Float32` LLRs via AFF3CT.jl's `ConvViterbiDecoder`
(issue #37). The deque-backed input boundary is identical to L1 C/A so the public
API is uniform.

# Fields

$(TYPEDFIELDS)
"""
struct GalileoINAVCache <: AbstractGNSSCache
    """
    Soft-symbol buffer (260 = 250 syncro + 10 preamble)
    """
    soft_buffer::CircularDeque{Float32}
    """
    Polarity-resolved 240-symbol FEC window copied out per decoded page part
    """
    soft_page::Vector{Float32}
    """
    Bits of the even page part of a word, held until its odd partner arrives;
    `nothing` when no even part is in flight
    """
    even_page_part_bits::Union{Nothing,UInt128}
    """
    Almanac of chain position 1, part-filled by one word type and completed by
    the next (word types 7-10 spread three satellites over four words)
    """
    almanac_chain_pos1::GalileoAlmanac
    """
    Almanac of chain position 2, filled the same way
    """
    almanac_chain_pos2::GalileoAlmanac
    """
    Viterbi decoder and its scratch buffers, built once and reused across pages
    (cf. the GPS L1C-D LDPC decoders)
    """
    viterbi::GalileoViterbiScratch
end

GalileoINAVCache() = GalileoINAVCache(
    CircularDeque{Float32}(GALILEO_INAV_PAGE_PART_SYMBOLS + GALILEO_INAV_SYNC_SYMBOLS),
    Vector{Float32}(undef, GALILEO_INAV_VITERBI_N),
    nothing,
    GalileoAlmanac(),
    GalileoAlmanac(),
    GalileoViterbiScratch(GALILEO_INAV_VITERBI_K, GALILEO_INAV_VITERBI_N),
)

function GalileoINAVCache(
    cache::GalileoINAVCache;
    soft_buffer = cache.soft_buffer,
    soft_page = cache.soft_page,
    even_page_part_bits = cache.even_page_part_bits,
    almanac_chain_pos1 = cache.almanac_chain_pos1,
    almanac_chain_pos2 = cache.almanac_chain_pos2,
    viterbi = cache.viterbi,
)
    GalileoINAVCache(
        soft_buffer,
        soft_page,
        even_page_part_bits,
        almanac_chain_pos1,
        almanac_chain_pos2,
        viterbi,
    )
end

# `soft_page` and the Viterbi handle are scratch, not state — both are
# overwritten from the deque before they are read — so they are excluded, the way
# `GalileoE6BCache` excludes its `fec_window`.
function Base.:(==)(a::GalileoINAVCache, b::GalileoINAVCache)
    deques_equal(a.soft_buffer, b.soft_buffer) &&
        a.even_page_part_bits == b.even_page_part_bits &&
        a.almanac_chain_pos1 == b.almanac_chain_pos1 &&
        a.almanac_chain_pos2 == b.almanac_chain_pos2
end

"""
    GalileoINAVData

Decoded Galileo I/NAV navigation message data — the shared container for the
two signal components that broadcast I/NAV, E1-B and E5b-I.

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
  - `M_0::Float64`: Mean anomaly at reference time (rad)
  - `e::Float64`: Eccentricity (dimensionless)
  - `sqrt_A::Float64`: Square root of semi-major axis (√m)
  - `Ω_0::Float64`: Longitude of ascending node at weekly epoch (rad)
  - `i_0::Float64`: Inclination angle at reference time (rad)
  - `ω::Float64`: Argument of perigee (rad)
  - `i_dot::Float64`: Rate of change of inclination angle (rad/s)
  - `Ω_dot::Float64`: Rate of change of right ascension (rad/s)
  - `Δn::Float64`: Mean motion difference from computed value (rad/s)
  - `C_uc::Float64`: Cosine harmonic correction to argument of latitude (rad)
  - `C_us::Float64`: Sine harmonic correction to argument of latitude (rad)
  - `C_rc::Float64`: Cosine harmonic correction to orbit radius (meters)
  - `C_rs::Float64`: Sine harmonic correction to orbit radius (meters)
  - `C_ic::Float64`: Cosine harmonic correction to inclination (rad)
  - `C_is::Float64`: Sine harmonic correction to inclination (rad)

# Signal-In-Space Accuracy (Word Type 3)

  - `SISA_E1_E5b::Int`: SISA index for dual frequency E1-E5b (Table 91/92; 255 = NAPA)

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

  - `E1B_SHS::SignalHealth`: E1-B/C signal health status (0=OK, 1=out of service, 2=Extended Operations Mode, 3=in test)
  - `E5b_SHS::SignalHealth`: E5b signal health status
  - `E1B_DVS::DataValidityStatus`: E1-B data validity (0=valid, 1=working without guarantee)
  - `E5b_DVS::DataValidityStatus`: E5b data validity

# Broadcast Group Delay (Word Type 5)

  - `BGD_E1_E5a::Float64`: E1-E5a group delay correction (seconds)
  - `BGD_E1_E5b::Float64`: E1-E5b group delay correction (seconds)

# Ionospheric Correction (Word Type 5)

  - `a_i0::Float64`: Effective Ionisation Level 1st-order coefficient (sfu)
  - `a_i1::Float64`: Effective Ionisation Level 2nd-order coefficient (sfu/degree)
  - `a_i2::Float64`: Effective Ionisation Level 3rd-order coefficient (sfu/degree²)
  - `iono_storm_flag_region1..5::Bool`: Ionospheric Disturbance (storm) flags for regions 1-5

# GST-UTC Conversion (Word Type 6)

  - `A_0UTC::Float64`: Constant term of polynomial (s; the ICD writes `A0`)
  - `A_1UTC::Float64`: 1st-order term of polynomial (s/s; the ICD writes `A1`)
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

    All four stay `nothing` when the satellite broadcasts the ICD's "GGTO not
    valid" encoding — every one of the four fields all ones (5.1.8).

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
Base.@kwdef struct GalileoINAVData <: AbstractGalileoEphemerisData
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

    SISA_E1_E5b::Union{Nothing,Int} = nothing

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

    E1B_SHS::Union{Nothing,SignalHealth} = nothing
    E5b_SHS::Union{Nothing,SignalHealth} = nothing
    E1B_DVS::Union{Nothing,DataValidityStatus} = nothing
    E5b_DVS::Union{Nothing,DataValidityStatus} = nothing

    BGD_E1_E5a::Union{Nothing,Float64} = nothing
    BGD_E1_E5b::Union{Nothing,Float64} = nothing

    a_i0::Union{Nothing,Float64} = nothing
    a_i1::Union{Nothing,Float64} = nothing
    a_i2::Union{Nothing,Float64} = nothing
    iono_storm_flag_region1::Union{Nothing,Bool} = nothing
    iono_storm_flag_region2::Union{Nothing,Bool} = nothing
    iono_storm_flag_region3::Union{Nothing,Bool} = nothing
    iono_storm_flag_region4::Union{Nothing,Bool} = nothing
    iono_storm_flag_region5::Union{Nothing,Bool} = nothing

    A_0UTC::Union{Nothing,Float64} = nothing
    A_1UTC::Union{Nothing,Float64} = nothing
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

function GalileoINAVData(
    data::GalileoINAVData;
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
    SISA_E1_E5b = data.SISA_E1_E5b,
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
    E1B_SHS = data.E1B_SHS,
    E5b_SHS = data.E5b_SHS,
    E1B_DVS = data.E1B_DVS,
    E5b_DVS = data.E5b_DVS,
    BGD_E1_E5a = data.BGD_E1_E5a,
    BGD_E1_E5b = data.BGD_E1_E5b,
    a_i0 = data.a_i0,
    a_i1 = data.a_i1,
    a_i2 = data.a_i2,
    iono_storm_flag_region1 = data.iono_storm_flag_region1,
    iono_storm_flag_region2 = data.iono_storm_flag_region2,
    iono_storm_flag_region3 = data.iono_storm_flag_region3,
    iono_storm_flag_region4 = data.iono_storm_flag_region4,
    iono_storm_flag_region5 = data.iono_storm_flag_region5,
    A_0UTC = data.A_0UTC,
    A_1UTC = data.A_1UTC,
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
    GalileoINAVData(
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
        SISA_E1_E5b,
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
        E1B_SHS,
        E5b_SHS,
        E1B_DVS,
        E5b_DVS,
        BGD_E1_E5a,
        BGD_E1_E5b,
        a_i0,
        a_i1,
        a_i2,
        iono_storm_flag_region1,
        iono_storm_flag_region2,
        iono_storm_flag_region3,
        iono_storm_flag_region4,
        iono_storm_flag_region5,
        A_0UTC,
        A_1UTC,
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
Base.:(==)(a::GalileoINAVData, b::GalileoINAVData) = fields_equal(a, b)

# `is_ephemeris_decoded` and `is_clock_correction_decoded` are per-constellation
# facts (identical fields for I/NAV and F/NAV), defined once on
# `AbstractGalileoData` in `galileo/galileo.jl`. Only the health-status check
# below is genuinely per-signal only in *which* facet it reports: word type 5
# carries E1-B/C and E5b health and data validity together, so an I/NAV decoder
# on either component always has both. `is_sat_healthy` picks the facet matching
# the component being decoded (see `e1b.jl` / `e5b.jl`).
function is_health_status_decoded(data::GalileoINAVData)
    !isnothing(data.E1B_SHS) &&
        !isnothing(data.E5b_SHS) &&
        !isnothing(data.E1B_DVS) &&
        !isnothing(data.E5b_DVS)
end

# I/NAV broadcasts the time of week as seconds (word types 0, 5 and 6).
get_time_of_week(data::GalileoINAVData) = data.TOW

"""
$(TYPEDSIGNATURES)

The GST-to-GPS time offset from I/NAV word type 10, or `nothing` when the
satellite is not broadcasting a valid GGTO or `target` is not `GPST()`.
See `galileo_ggto_offset`.
"""
get_time_offset(state::GNSSDecoderState{<:GalileoINAVData}, target::TimeSystem) =
    galileo_ggto_offset(state, target)

function is_decoding_completed_for_positioning(data::GalileoINAVData)
    !isnothing(data.TOW) &&
        !isnothing(data.WN) &&
        !isnothing(data.BGD_E1_E5a) &&
        !isnothing(data.BGD_E1_E5b) &&
        is_ephemeris_decoded(data) &&
        is_clock_correction_decoded(data) &&
        is_health_status_decoded(data)
end

"""
$(TYPEDSIGNATURES)

Reset a Galileo I/NAV decoder state after a signal loss or reacquisition.

Clears the bit buffers and time-of-week (TOW) field while preserving other
decoded ephemeris and clock data in `raw_data`. This allows faster recovery
after brief signal outages without requiring a full re-decode of all pages.

!!! note

    The week number (`WN`) field is intentionally not reset as it is not
    broadcast as frequently as TOW. This may cause brief errors if a week
    rollover occurs during a signal outage.

# Arguments

  - `state::GNSSDecoderState{<:GalileoINAVData}`: Current Galileo I/NAV (E1-B or E5b-I) decoder state

# Returns

  - `GNSSDecoderState{<:GalileoINAVData}`: Reset decoder state with cleared buffers

# Example

```julia
# After detecting signal loss
state = reset_decoder_state(state)
# Continue decoding with preserved ephemeris
state = decode(state, new_bits, num_bits)
```

# See Also

  - [`GalileoE1BDecoderState`](@ref) / [`GalileoE5bDecoderState`](@ref): Create a fresh decoder state
  - [`decode`](@ref): Continue decoding after reset
"""
function reset_decoder_state(state::GNSSDecoderState{<:GalileoINAVData})
    # Reset bit buffers and TOW data field, while keeping the
    # remaining parameters in raw_data. This allows a GNSSReceiver
    # to use a satellite after a reacquisition without waiting for
    # the decoding of all data fields.
    # Note: WN is currently not reset as it is broadcast not as
    # frequently as the TOW and thus may increase the time until
    # the decoder is available again after an outage. This will
    # lead to erroneous decoder information for a few seconds after
    # reacquisition when a new week started during a signal outage.
    empty!(state.cache.soft_buffer)
    GNSSDecoderState(
        state;
        raw_data = GalileoINAVData(
            state.raw_data;
            TOW = nothing,
            num_bits_after_valid_syncro_sequence_after_last_TOW = nothing,
        ),
        data = GalileoINAVData(),
        num_bits_after_valid_syncro_sequence = nothing,
    )
end

# No `packed_buffer_type` method: I/NAV overrides `try_sync` and reads the two
# 10-symbol page-sync patterns straight from the soft buffer, so the 260-bit
# packed window the default would build once per symbol is never needed — the
# payload is consumed as soft symbols by the Viterbi decoder below.
"""
    try_sync(state::GNSSDecoderState{<:GalileoINAVData}) -> Union{Nothing,Bool}

I/NAV page-part sync: the 10-symbol pattern `0101100000` must appear at both ends
of the 260-symbol window (start of this page part and start of the next), both
upright or both inverted. Returns the resolved polarity, or `nothing` when there
is no sync.
"""
try_sync(state::GNSSDecoderState{<:GalileoINAVData}) = find_preamble_in_deque(
    soft_buffer(state),
    state.constants.preamble,
    GALILEO_INAV_SYNC_SYMBOLS,
    GALILEO_INAV_PAGE_PART_SYMBOLS,
)

# Record the polarity `try_sync` resolved; `decode_syncro_sequence` reads it back
# off the state when it copies the soft page part out of the deque.
function complement_buffer_if_necessary(
    state::GNSSDecoderState{<:GalileoINAVData},
    polarity_flipped::Bool,
)
    GNSSDecoderState(state; is_shifted_by_180_degrees = polarity_flipped), polarity_flipped
end

function decode_syncro_sequence(state::GNSSDecoderState{<:GalileoINAVData}, ::Bool)
    # The 240 encoded symbols are the soft-buffer entries between the leading
    # 10-bit page-sync preamble and the trailing preamble of the next page
    # (deque indices preamble_length+1 .. syncro_sequence_length). Resolve the
    # 180-degree polarity ambiguity by negating the LLRs when the sync hook
    # flagged the page as inverted (an inverted soft symbol = a negated one).
    soft_page = copy_soft_window!(
        state.cache.soft_page,
        soft_buffer(state),
        state.constants.preamble_length,
        GALILEO_INAV_VITERBI_N,
        state.is_shifted_by_180_degrees,
    )
    bits = galileo_viterbi(
        state.cache.viterbi,
        soft_page,
        GALILEO_INAV_INTERLEAVER_COLUMNS,
        UInt128,
    )
    is_even = !get_bit(bits, 114, 1)
    is_nominal_page = !get_bit(bits, 114, 2)
    state = GNSSDecoderState(
        state;
        raw_data = GalileoINAVData(
            state.raw_data;
            num_pages_after_last_TOW = state.raw_data.num_pages_after_last_TOW + 1,
        ),
    )
    if is_even
        state = GNSSDecoderState(
            state;
            cache = GalileoINAVCache(
                state.cache;
                even_page_part_bits = is_nominal_page ? bits : nothing,
            ),
        )
        return state
    end
    even_bits = state.cache.even_page_part_bits
    # `nothing`, not a zero sentinel: a nominal even page part whose 112 content
    # bits are all zero is legal (word type 0 with the Time field '00' and a zero
    # spare) and a `!= 0` test would silently discard its odd partner.
    if !isnothing(even_bits) && is_nominal_page
        data = get_bits(even_bits, 114, 3, 112) << 16 + get_bits(bits, 114, 3, 16)
        bits_to_check_CRC = UInt288(even_bits) << 106 + get_bits(bits, 114, 1, 106)
        if crc24q(bits_to_check_CRC, GALILEO_INAV_CRC_OCTETS) == 0
            data_type = get_bits(data, 128, 1, 6)
            if data_type == 0
                if get_bits(data, 128, 7, 2) == 2 # '10'
                    WN = get_bits(data, 128, 97, 12)
                    TOW = get_bits(data, 128, 109, 20)
                    state = GNSSDecoderState(
                        state;
                        raw_data = GalileoINAVData(
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
                    get_twos_complement_num(data, 128, 31, 32) *
                    state.constants.PI *
                    2.0^-31
                e = get_bits(data, 128, 63, 32) * 2.0^-33
                sqrt_A = get_bits(data, 128, 95, 32) / 1 << 19
                state = GNSSDecoderState(
                    state;
                    raw_data = GalileoINAVData(
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
                    get_twos_complement_num(data, 128, 17, 32) *
                    state.constants.PI *
                    2.0^-31
                i_0 =
                    get_twos_complement_num(data, 128, 49, 32) *
                    state.constants.PI *
                    2.0^-31
                ω =
                    get_twos_complement_num(data, 128, 81, 32) *
                    state.constants.PI *
                    2.0^-31
                i_dot =
                    get_twos_complement_num(data, 128, 113, 14) *
                    state.constants.PI *
                    2.0^-43
                state = GNSSDecoderState(
                    state;
                    raw_data = GalileoINAVData(
                        state.raw_data;
                        IOD_nav2,
                        Ω_0,
                        i_0,
                        ω,
                        i_dot,
                    ),
                )
            elseif data_type == 3
                IOD_nav3 = get_bits(data, 128, 7, 10)
                Ω_dot =
                    get_twos_complement_num(data, 128, 17, 24) *
                    state.constants.PI *
                    2.0^-43
                Δn =
                    get_twos_complement_num(data, 128, 41, 16) *
                    state.constants.PI *
                    2.0^-43
                C_uc = get_twos_complement_num(data, 128, 57, 16) / 1 << 29
                C_us = get_twos_complement_num(data, 128, 73, 16) / 1 << 29
                C_rc = get_twos_complement_num(data, 128, 89, 16) / 1 << 5
                C_rs = get_twos_complement_num(data, 128, 105, 16) / 1 << 5
                SISA_E1_E5b = Int(get_bits(data, 128, 121, 8))
                state = GNSSDecoderState(
                    state;
                    raw_data = GalileoINAVData(
                        state.raw_data;
                        IOD_nav3,
                        Ω_dot,
                        Δn,
                        C_uc,
                        C_us,
                        C_rc,
                        C_rs,
                        SISA_E1_E5b,
                    ),
                )
            elseif data_type == 4
                IOD_nav4 = get_bits(data, 128, 7, 10)
                SVID = Int(get_bits(data, 128, 17, 6))
                C_ic = get_twos_complement_num(data, 128, 23, 16) / 1 << 29
                C_is = get_twos_complement_num(data, 128, 39, 16) / 1 << 29
                t_0c = get_bits(data, 128, 55, 14) * 60
                a_f0 = get_twos_complement_num(data, 128, 69, 31) * 2.0^-34
                a_f1 = get_twos_complement_num(data, 128, 100, 21) * 2.0^-46
                a_f2 = get_twos_complement_num(data, 128, 121, 6) * 2.0^-59
                state = GNSSDecoderState(
                    state;
                    raw_data = GalileoINAVData(
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
                BGD_E1_E5a = get_twos_complement_num(data, 128, 48, 10) * 2.0^-32
                BGD_E1_E5b = get_twos_complement_num(data, 128, 58, 10) * 2.0^-32
                E5b_SHS = SignalHealth(get_bits(data, 128, 68, 2))
                E1B_SHS = SignalHealth(get_bits(data, 128, 70, 2))
                E5b_DVS = DataValidityStatus(get_bit(data, 128, 72))
                E1B_DVS = DataValidityStatus(get_bit(data, 128, 73))
                WN = get_bits(data, 128, 74, 12)
                TOW = get_bits(data, 128, 86, 20)
                state = GNSSDecoderState(
                    state;
                    raw_data = GalileoINAVData(
                        state.raw_data;
                        a_i0,
                        a_i1,
                        a_i2,
                        iono_storm_flag_region1,
                        iono_storm_flag_region2,
                        iono_storm_flag_region3,
                        iono_storm_flag_region4,
                        iono_storm_flag_region5,
                        BGD_E1_E5a,
                        BGD_E1_E5b,
                        E5b_SHS,
                        E1B_SHS,
                        E5b_DVS,
                        E1B_DVS,
                        WN,
                        TOW,
                        num_pages_after_last_TOW = 1,
                        num_bits_after_valid_syncro_sequence_after_last_TOW = state.num_bits_after_valid_syncro_sequence,
                    ),
                )
            elseif data_type == 6
                A_0UTC = get_twos_complement_num(data, 128, 7, 32) / 1 << 30
                A_1UTC = get_twos_complement_num(data, 128, 39, 24) * 2.0^-50
                Δt_LS = Int(get_twos_complement_num(data, 128, 63, 8))
                t_0t = Int(get_bits(data, 128, 71, 8) * 3600)
                WN_0t = Int(get_bits(data, 128, 79, 8))
                WN_LSF = Int(get_bits(data, 128, 87, 8))
                DN = Int(get_bits(data, 128, 95, 3))
                Δt_LSF = Int(get_twos_complement_num(data, 128, 98, 8))
                TOW = get_bits(data, 128, 106, 20)
                state = GNSSDecoderState(
                    state;
                    raw_data = GalileoINAVData(
                        state.raw_data;
                        A_0UTC,
                        A_1UTC,
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
                    get_twos_complement_num(data, 128, 96, 11) *
                    state.constants.PI *
                    2.0^-33
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
                    cache = GalileoINAVCache(
                        state.cache;
                        almanac_chain_pos1 = valid_SVID ? almanac_pos1 : GalileoAlmanac(),
                        almanac_chain_pos2 = GalileoAlmanac(),
                    ),
                )
            elseif data_type == 8
                IOD_a = Int(get_bits(data, 128, 7, 4))
                a_f0_pos1 = get_twos_complement_num(data, 128, 11, 16) / 1 << 19
                a_f1_pos1 = get_twos_complement_num(data, 128, 27, 13) * 2.0^-38
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
                    get_twos_complement_num(data, 128, 117, 11) *
                    state.constants.PI *
                    2.0^-33
                # Flush position-1 almanac if its WT7 partial is intact and IODa matches
                completed_pos1 = GalileoAlmanac(
                    state.cache.almanac_chain_pos1;
                    a_f0 = a_f0_pos1,
                    a_f1 = a_f1_pos1,
                    E5b_SHS = signal_health_e5b_pos1,
                    E1B_SHS = signal_health_e1b_pos1,
                )
                almanacs = state.raw_data.almanacs
                if completed_pos1.IOD_a == IOD_a && !isnothing(completed_pos1.SVID)
                    almanacs = _merge_keyed(almanacs, completed_pos1.SVID, completed_pos1)
                end
                almanac_pos2 = GalileoAlmanac(; SVID, Δsqrt_A, e, ω, δi, Ω_0, Ω_dot, IOD_a)
                valid_SVID = SVID >= 1
                state = GNSSDecoderState(
                    state;
                    raw_data = GalileoINAVData(state.raw_data; almanacs),
                    cache = GalileoINAVCache(
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
                a_f1_pos2 = get_twos_complement_num(data, 128, 55, 13) * 2.0^-38
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
                    E5b_SHS = signal_health_e5b_pos2,
                    E1B_SHS = signal_health_e1b_pos2,
                    WN_a,
                    t_0a,
                )
                almanacs = state.raw_data.almanacs
                if completed_pos2.IOD_a == IOD_a && !isnothing(completed_pos2.SVID)
                    almanacs = _merge_keyed(almanacs, completed_pos2.SVID, completed_pos2)
                end
                almanac_pos3 = GalileoAlmanac(; SVID, Δsqrt_A, e, ω, δi, IOD_a, WN_a, t_0a)
                valid_SVID = SVID >= 1
                state = GNSSDecoderState(
                    state;
                    raw_data = GalileoINAVData(state.raw_data; almanacs),
                    cache = GalileoINAVCache(
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
                    get_twos_complement_num(data, 128, 27, 11) *
                    state.constants.PI *
                    2.0^-33
                M_0_pos3 =
                    get_twos_complement_num(data, 128, 38, 16) * state.constants.PI /
                    1 << 15
                a_f0_pos3 = get_twos_complement_num(data, 128, 54, 16) / 1 << 19
                a_f1_pos3 = get_twos_complement_num(data, 128, 70, 13) * 2.0^-38
                signal_health_e5b_pos3 = SignalHealth(get_bits(data, 128, 83, 2))
                signal_health_e1b_pos3 = SignalHealth(get_bits(data, 128, 85, 2))
                # GGTO — all four fields all-ones means "not valid" (ICD
                # 5.1.8), so they are read raw and scaled by `galileo_ggto`.
                A_0G, A_1G, t_0G, WN_0G = galileo_ggto(
                    get_bits(data, 128, 87, 16),
                    get_bits(data, 128, 103, 12),
                    get_bits(data, 128, 115, 8),
                    get_bits(data, 128, 123, 6),
                )
                # Position 3 partial sits in cache pos1 (set by WT9)
                completed_pos3 = GalileoAlmanac(
                    state.cache.almanac_chain_pos1;
                    Ω_0 = Ω_0_pos3,
                    Ω_dot = Ω_dot_pos3,
                    M_0 = M_0_pos3,
                    a_f0 = a_f0_pos3,
                    a_f1 = a_f1_pos3,
                    E5b_SHS = signal_health_e5b_pos3,
                    E1B_SHS = signal_health_e1b_pos3,
                )
                almanacs = state.raw_data.almanacs
                if completed_pos3.IOD_a == IOD_a && !isnothing(completed_pos3.SVID)
                    almanacs = _merge_keyed(almanacs, completed_pos3.SVID, completed_pos3)
                end
                state = GNSSDecoderState(
                    state;
                    raw_data = GalileoINAVData(
                        state.raw_data;
                        almanacs,
                        A_0G,
                        A_1G,
                        t_0G,
                        WN_0G,
                    ),
                    cache = GalileoINAVCache(
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
                a_f1_red = get_twos_complement_num(data, 128, 123, 6) * 2.0^-35
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
                    raw_data = GalileoINAVData(state.raw_data; reduced_ced),
                )
            end
        end
    end
    return state
end

function validate_data(state::GNSSDecoderState{<:GalileoINAVData})
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
            num_bits_after_valid_syncro_sequence =
                state.num_bits_after_valid_syncro_sequence - (
                    state.raw_data.num_bits_after_valid_syncro_sequence_after_last_TOW -
                    2 * state.constants.syncro_sequence_length -
                    state.constants.preamble_length
                )
        else # first succesful decoding
            num_bits_after_valid_syncro_sequence =
                state.constants.preamble_length +
                (state.raw_data.num_pages_after_last_TOW + 1) *
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
