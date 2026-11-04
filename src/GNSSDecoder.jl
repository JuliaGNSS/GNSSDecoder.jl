module GNSSDecoder

using DocStringExtensions, GNSSSignals, BitIntegers, Dictionaries, DataStructures
using StaticStrings: AbstractStaticString, CStaticString
using GNSSSignals: SECONDS_PER_WEEK
using Unitful: Hz, ustrip
using Dates: Dates, Date, DateTime
import Aff3ct

export decode!,
    SlotDictionary,
    GPSL1CADecoderState,
    GPSLNAVData,
    GPSL1CAAlmanac,
    GPSL1C_DDecoderState,
    GPSCNAV2Data,
    GPSL1C_DReducedAlmanac,
    GPSL1C_DMidiAlmanac,
    GPSL1C_DDifferentialCorrection,
    GPSL5IDecoderState,
    GPSCNAVData,
    GPSCNAVReducedAlmanac,
    GPSCNAVMidiAlmanac,
    GPSCNAVClockDifferentialCorrection,
    GPSCNAVEphemerisDifferentialCorrection,
    GPSCNAVIntegritySupportMessage,
    GPSL2CMDecoderState,
    GalileoE1BDecoderState,
    GalileoE5aDecoderState,
    GalileoE5bDecoderState,
    GalileoE6BDecoderState,
    GalileoINAVData,
    GalileoReducedCED,
    GalileoFNAVData,
    GalileoE6BData,
    GalileoHASMessage,
    GalileoHASMask,
    GalileoHASSatelliteMask,
    GalileoHASCorrectionBlock,
    GalileoHASOrbitCorrection,
    GalileoHASClockCorrection,
    GalileoHASCodeBias,
    GalileoHASPhaseBias,
    GalileoAlmanac,
    SignalHealth,
    DataValidityStatus,
    HASStatus,
    BeiDouB1IDecoderState,
    BeiDouB3IDecoderState,
    BeiDouDNAVData,
    BeiDouDNAVAlmanac,
    BeiDouB1CDecoderState,
    BeiDouCNAV1Data,
    BeiDouB1CBGTO,
    BeiDouB2aDecoderState,
    BeiDouCNAV2Data,
    BeiDouB2bDecoderState,
    BeiDouCNAV3Data,
    BeiDouMidiAlmanac,
    BeiDouReducedAlmanac,
    is_sat_healthy,
    has_positioning_data,
    get_signal_type,
    OrbitClass,
    get_orbit_class,
    get_broadcast_time_of_week,
    get_decoded_time_of_week,
    get_broadcast_week,
    get_broadcast_inter_system_time_offset,
    get_broadcast_inter_system_time_offset_epoch,
    UTCParameters,
    LeapSecondParameters,
    InterSystemTimeOffsetParameters,
    get_utc_parameters,
    get_leap_second_parameters,
    get_inter_system_time_offset_parameters,
    GNSSDecoderState,
    reset_decoder_state!

# Shared decode framework, cross-signal accessors and bit-extraction helpers.
include("gnss.jl")
include("storage.jl")
include("bit_fiddling.jl")

# Signal-independent channel coding: CRC-24Q, BCH(51,8) TOI, block (de)interleaver, LDPC
# (consumes `crc24q`, so `crc.jl` comes first) and GF(2^8) Reed-Solomon.
include("coding/crc.jl")
include("coding/bch_toi.jl")
include("coding/deinterleave.jl")
include("coding/ldpc.jl")
include("coding/reed_solomon.jl")

# GPS. `gps.jl`: shared definitions; `cnav.jl`: CNAV core (identical on L5I and L2C),
# with `l5.jl` and `l2c.jl` as thin signal layers.
include("gps/gps.jl")
include("gps/l1ca.jl")
include("gps/l1c_d.jl")
include("gps/cnav.jl")
include("gps/l5.jl")
include("gps/l2c.jl")

# Galileo. `galileo.jl`: shared definitions; `inav.jl`: I/NAV core (identical on E1-B and
# E5b-I), with `e1b.jl` and `e5b.jl` as thin signal layers.
include("galileo/galileo.jl")
include("galileo/inav.jl")
include("galileo/e1b.jl")
include("galileo/e5b.jl")
include("galileo/e5a.jl")
include("galileo/e6b.jl")

# BeiDou. `beidou.jl`: shared definitions; `dnav.jl`: D1/D2 legacy NAV core (identical on
# B1I and B3I), with `b1i.jl` and `b3i.jl` as thin signal layers.
include("beidou/beidou.jl")
include("beidou/dnav.jl")
include("beidou/b1i.jl")
include("beidou/b3i.jl")
include("beidou/b1c.jl")
include("beidou/b2a.jl")
include("beidou/b2b.jl")
include("time_parameters.jl")
include("precompile.jl")
end
