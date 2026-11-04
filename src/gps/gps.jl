# Definitions shared across GPS signals (counterpart of `galileo/galileo.jl` and
# `beidou/beidou.jl`). Included before the per-signal GPS decoders.

# WGS-84 constants (IS-GPS-200N §20.3.3.4.3; identical in IS-GPS-705J and IS-GPS-800J).
const GPS_μ = 3.986005e14        # Earth gravitational parameter (m³/s²)
const GPS_F = -4.442807633e-10   # relativistic correction constant (s/√m)

# GGTO Type ID (symbol `GGTO ID` since IRN-003; CNAV message type 35, IS-GPS-705J
# §20.3.3.8.1 / IS-GPS-200N §30.3.3.8.1; CNAV-2 subframe 3 page 2, IS-GPS-800J
# §3.5.4.2.1): 000 = no data, 001 = Galileo, 010 = GLONASS, 011-111 reserved
# (unusable). GLONASS is unreachable here: `GNSSSignals` defines no GLONASS `TimeSystem`.
const GPS_GGTO_ID_GALILEO = 1

# `WN_GGTO` is 13 bits, like the CNAV / CNAV-2 week number, so resolving it only
# matters across a rollover. See `resolve_reference_week`.
const GPS_GGTO_WN_MODULUS = 8192

"""
    gps_ggto_offset(state::GNSSDecoderState, target::TimeSystem) -> Union{Nothing,TimeOffsetParameters}

Normalise the GGTO parameter set of CNAV (message type 35) and CNAV-2
(subframe 3 page 2) into a `TimeOffsetParameters`, or `nothing` when the set is
absent, marked unavailable, or refers to a system other than `target`. GPS
counterpart of `beidou_bgto_offset`.

Only one set is carried at a time, tagged by `GGTO_ID`, so which target is
answerable is the satellite's choice.
"""
function gps_ggto_offset(state::GNSSDecoderState, target::TimeSystem)
    data = state.data
    # Only ID 001 is reachable; 000 is "no data", reserved codes are unusable.
    (data.GGTO_ID == GPS_GGTO_ID_GALILEO && target === GST()) || return nothing
    # Decoded together, so all present or all `nothing`; each is checked anyway to
    # narrow its type for the call below (see `galileo_ggto_offset`).
    (; A_0GGTO, A_1GGTO, A_2GGTO, t_GGTO, WN_GGTO) = data
    (
        isnothing(A_0GGTO) ||
        isnothing(A_1GGTO) ||
        isnothing(A_2GGTO) ||
        isnothing(t_GGTO) ||
        isnothing(WN_GGTO)
    ) && return nothing
    assemble_time_offset_parameters(
        target,
        A_0GGTO,
        A_1GGTO,
        A_2GGTO,
        t_GGTO,
        WN_GGTO,
        data.WN,
        GPS_GGTO_WN_MODULUS,
    )
end

# ---- Integrity Support Message (shared by CNAV and CNAV-2) ------------------
#
# IRN-IS-200N-004 defines one 236-bit ISM Packet (IS-GPS-200 §30.3.3.10); CNAV message
# type 40 (IRN-IS-705J-004 / IRN-IS-200N-004) and CNAV-2 subframe 3 page 8
# (IRN-IS-800J-004) carry it unchanged, so both decoders share the record and parser.

"""
    GPSIntegritySupportMessage

One GPS ISM Packet (GNSS ID `0100`; ARAIM integrity support data), complete in a
single CNAV message type 40 or CNAV-2 subframe 3 page 8, and the record of both
`GPSCNAVData.ism` and `GPSCNAV2Data.ism`. The packet is defined in IS-GPS-200
§30.3.3.10 (Figures 30-17 and 30-19, Table 30-XIa). IRN-IS-705J-004 and IRN-IS-200N-004
put it at message bits 39-274 of message type 40 in place of IS-GPS-705J Table 20-XIa
(packet bit `n` is message bit `38 + n`; see `parse_mt40`); IRN-IS-800J-004 puts it at
bits 15-250 of subframe 3 page 8 in place of IS-GPS-800J Table 3.5-9 (packet bit `n` is
page bit `14 + n`; see `parse_sf3_page8`). Only that layout is decoded.

Lookup-coded fields are converted to SI units (hours to seconds, per hour to per
second); a code the ICD marks Reserved decodes as `nothing`. Each packet applies to the
PRNs of its own inclusion mask at its own Service Level, so a satellite may broadcast
several. The decoder keeps one: the valid packet that arrived last. That is a known
limitation: the ICD has users take the parameters with the most recent `WN_ISM` /
`TOW_ISM` time stamp (IS-GPS-200 §30.3.3.10.1.2.2.2), and packets for different Service
Levels or inclusion masks overwrite each other here.

# Fields

  - `WN_ISM::Int`: ISM week number at the start of validity (13 bits, weeks).
  - `TOW_ISM::Int`: ISM time of week at the start of validity, in seconds
    (6 bits, LSB 4 hours, valid range 0-164 hours).
  - `t_correl::Union{Nothing,Float64}`: Correlation time constant of the errors the
    URA characterises (s; 0.25-12 hours).
  - `β_nom::Float64`: IAURA-independent additive term of the nominal pseudorange
    error bias (m; LSB 0.1 m, 0-1.5 m).
  - `γ_nom::Float64`: Scalar term of the nominal bias, scaling with IAURA
    (dimensionless; LSB 0.05, 0-0.75). The bias is
    `b_nom = β_nom + γ_nom · IAURA` (IS-GPS-200 §30.3.3.10.2).
  - `R_sat::Union{Nothing,Float64}`: Satellite fault rate, the onset rate of an
    instantaneous URE beyond ±4.42 IAURA on one satellite (1/s; 1e-8 to 1e-4 per
    hour).
  - `MFD_sat::Union{Nothing,Float64}`: Mean duration of a satellite fault without a
    timely notification (s; 0.25-8 hours).
  - `R_const::Union{Nothing,Float64}`: Constellation fault rate, the same onset rate
    for two or more satellites from a common cause (1/s; 3.16e-10 to 3.16e-6 per
    hour).
  - `MFD_const::Union{Nothing,Float64}`: Mean duration of a constellation fault
    without a timely notification issued to the user (s; 0.25-8 hours).
  - `Service_Level::Union{Nothing,Int}`: Service Level 1-5 of Table 30-XIb (code
    `000` is Level 1): 1 no data available (earlier parameters are invalid),
    2 non-safety-of-life, 3 H-ARAIM, 4 V-ARAIM, 5 both.
  - `GPS_PRN_Inclusion_Mask::UInt64`: 63-bit mask, PRN 1 in the most significant
    bit (packet bit 55), PRN 63 in the least; a set bit means the packet applies to
    that PRN (Table 30-XIc).
  - `ISM_CRC::UInt32`: The packet's 32-bit CRC, already verified against packet
    bits 1-204.

The 87 bits Reserved For ISM are not kept.

# Reference

IRN-IS-705J-004 (Figure 20-14a, §20.3.3.10); IRN-IS-800J-004 (Figure 3.5-8a,
§3.5.4.7, §3.5.4.7.2); IS-GPS-200 as amended by IRN-IS-200N-004, §30.3.3.10.1.2.2,
Figures 30-17 and 30-19, Tables 30-XIa to 30-XIc.
"""
Base.@kwdef struct GPSIntegritySupportMessage
    WN_ISM::Int
    TOW_ISM::Int
    t_correl::Union{Nothing,Float64}
    β_nom::Float64
    γ_nom::Float64
    R_sat::Union{Nothing,Float64}
    MFD_sat::Union{Nothing,Float64}
    R_const::Union{Nothing,Float64}
    MFD_const::Union{Nothing,Float64}
    Service_Level::Union{Nothing,Int}
    GPS_PRN_Inclusion_Mask::UInt64
    ISM_CRC::UInt32
end

# GPS ISM Packet lookup tables (IS-GPS-200 §30.3.3.10.1.2.2.3-.9 as amended by
# IRN-IS-200N-004), in the ICD's hours and per-hour values, indexed by code + 1; the
# codes past the end of each tuple are Reserved.
const ISM_T_CORREL_HOURS = (
    0.25, 0.33, 0.50, 0.67, 0.83, 1.00, 1.17, 1.33,
    1.50, 2.10, 3.00, 4.20, 6.00, 8.50, 12.00,
)
const ISM_R_SAT_PER_HOUR =
    (1.00e-8, 3.16e-8, 1.00e-7, 3.16e-7, 1.00e-6, 3.16e-6, 1.00e-5, 3.16e-5, 1.00e-4)
const ISM_R_CONST_PER_HOUR =
    (3.16e-10, 1.00e-9, 3.16e-9, 1.00e-8, 3.16e-8, 1.00e-7, 3.16e-7, 1.00e-6, 3.16e-6)
const ISM_MFD_HOURS = (0.25, 0.5, 1.0, 2.0, 4.0, 6.0, 8.0)  # MFD_sat and MFD_const
const ISM_SERVICE_LEVELS = 5  # codes 000-100 are Levels 1-5, 101-111 Reserved

const ISM_GNSS_ID_GPS = 0b0100
const ISM_PACKET_BITS = 236

# Table entry for `code` converted to SI, or `nothing` for a Reserved code.
function _ism_hours(table, code)
    hours = get(table, code + 1, nothing)
    isnothing(hours) ? nothing : 3600.0 * hours
end
function _ism_per_hour(table, code)
    rate = get(table, code + 1, nothing)
    isnothing(rate) ? nothing : rate / 3600.0
end

"""
    gps_ism_packet(word, word_length, start)
        -> Union{Nothing,GPSIntegritySupportMessage}

Decode the 236-bit ISM Packet whose bit 1 is bit `start` of `word` (IS-GPS-200
§30.3.3.10.1, Figures 30-17 and 30-19, Table 30-XIa, as amended by IRN-IS-200N-004).
Packet bit `n` sits at `start + n - 1`:

| Packet bits | Field |
| --- | --- |
| 1-4 | GNSS ID (`0100` GPS, `0000` test packet, others reserved) |
| 5-17 | `WN_ISM` |
| 18-23 | `TOW_ISM` |
| 24-27 | `t_correl` |
| 28-31 | `β_nom` |
| 32-35 | `γ_nom` |
| 36-39 | `R_sat` |
| 40-43 | `MFD_sat` |
| 44-47 | `R_const` |
| 48-51 | `MFD_const` |
| 52-54 | Service Level |
| 55-117 | GPS PRN Inclusion Mask |
| 118-204 | Reserved For ISM |
| 205-236 | ISM CRC over bits 1-204 |

Returns `nothing` unless the GNSS ID is GPS and the ISM CRC matches: the ICD defines no
other overlay, tells users to ignore a test packet, and calls the ISD valid only when
the CRC matches.
Shared by CNAV message type 40 (`start = 39`) and CNAV-2 subframe 3 page 8
(`start = 15`, IRN-IS-800J-004 Figure 3.5-8a).
"""
function gps_ism_packet(word::Unsigned, word_length::Int, start::Int)
    bits(n, len) = Int(get_bits(word, word_length, start + n - 1, len))
    bits(1, 4) == ISM_GNSS_ID_GPS || return nothing
    crc32q(word, word_length, start, ISM_PACKET_BITS) == 0 || return nothing
    service_level = bits(52, 3) + 1
    GPSIntegritySupportMessage(;
        WN_ISM = bits(5, 13),
        TOW_ISM = bits(18, 6) * 4 * 3600,
        t_correl = _ism_hours(ISM_T_CORREL_HOURS, bits(24, 4)),
        β_nom = bits(28, 4) * 0.1,
        γ_nom = bits(32, 4) * 0.05,
        R_sat = _ism_per_hour(ISM_R_SAT_PER_HOUR, bits(36, 4)),
        MFD_sat = _ism_hours(ISM_MFD_HOURS, bits(40, 4)),
        R_const = _ism_per_hour(ISM_R_CONST_PER_HOUR, bits(44, 4)),
        MFD_const = _ism_hours(ISM_MFD_HOURS, bits(48, 4)),
        Service_Level = service_level <= ISM_SERVICE_LEVELS ? service_level : nothing,
        GPS_PRN_Inclusion_Mask = UInt64(get_bits(word, word_length, start + 54, 63)),
        ISM_CRC = UInt32(get_bits(word, word_length, start + 204, 32)),
    )
end
