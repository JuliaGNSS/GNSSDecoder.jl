module GNSSDecoder

using DocStringExtensions, GNSSSignals, BitIntegers, Dictionaries, DataStructures
using StaticStrings: AbstractStaticString, CStaticString
using GNSSSignals: SECONDS_PER_WEEK
using Unitful: Hz, ustrip
using Dates: Dates, Date
import Aff3ct

export decode!,
    GPSL1CADecoderState,
    GPSLNAVData,
    GPSL1C_DDecoderState,
    GPSCNAV2Data,
    GPSL5IDecoderState,
    GPSCNAVData,
    GPSL2CMDecoderState,
    GalileoE1BDecoderState,
    GalileoE5aDecoderState,
    GalileoE5bDecoderState,
    GalileoE6BDecoderState,
    GalileoINAVData,
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
    SignalHealth,
    DataValidityStatus,
    HASStatus,
    BeiDouB1IDecoderState,
    BeiDouB3IDecoderState,
    BeiDouDNAVData,
    BeiDouB1CDecoderState,
    BeiDouCNAV1Data,
    BeiDouB2aDecoderState,
    BeiDouCNAV2Data,
    BeiDouB2bDecoderState,
    BeiDouCNAV3Data,
    is_sat_healthy,
    has_positioning_data,
    get_signal_type,
    OrbitClass,
    get_orbit_class,
    get_broadcast_time_of_week,
    get_decoded_time_of_week,
    get_broadcast_week,
    get_broadcast_time_offset,
    get_broadcast_time_offset_epoch,
    get_group_delay,
    BroadcastEphemeris,
    get_broadcast_ephemeris,
    BroadcastEphemerisValidity,
    get_broadcast_ephemeris_validity,
    KlobucharParameters,
    BeiDouKlobucharParameters,
    NeQuickGParameters,
    BDGIMParameters,
    get_klobuchar_parameters,
    get_nequick_g_parameters,
    get_bdgim_parameters,
    UTCParameters,
    LeapSecondParameters,
    GGTOParameters,
    get_utc_parameters,
    get_leap_second_parameters,
    get_ggto_parameters,
    get_range_accuracy,
    GNSSDecoderState,
    reset_decoder_state!

# Shared decode framework and cross-signal accessors, then the generic
# bit-extraction helpers.
include("gnss.jl")
include("storage.jl")
include("bit_fiddling.jl")

# Signal-independent channel-coding primitives: CRC-24Q, the BCH(51,8) TOI
# codec, the block (de)interleaver, the shared LDPC pipeline (which consumes
# `crc24q`, so `crc.jl` comes first) and the GF(2^8) Reed-Solomon codec. Each
# carries its own unit tests under `test/`.
include("coding/crc.jl")
include("coding/bch_toi.jl")
include("coding/deinterleave.jl")
include("coding/ldpc.jl")
include("coding/reed_solomon.jl")

# GPS. `gps.jl` holds definitions shared across the GPS signals; `cnav.jl` is
# the shared CNAV core (broadcast identically on L5I and L2C), with `l5.jl` and
# `l2c.jl` as thin signal layers over it.
include("gps/gps.jl")
include("gps/l1ca.jl")
include("gps/l1c_d.jl")
include("gps/cnav.jl")
include("gps/l5.jl")
include("gps/l2c.jl")

# Galileo. `galileo.jl` holds definitions shared across the Galileo signals;
# `inav.jl` is the shared I/NAV core (broadcast identically on E1-B and E5b-I),
# with `e1b.jl` and `e5b.jl` as thin signal layers over it.
include("galileo/galileo.jl")
include("galileo/inav.jl")
include("galileo/e1b.jl")
include("galileo/e5b.jl")
include("galileo/e5a.jl")
include("galileo/e6b.jl")

# BeiDou. `beidou.jl` holds definitions shared across the BeiDou signals;
# `dnav.jl` is the shared D1/D2 legacy NAV core (broadcast identically on B1I
# and B3I), with `b1i.jl` and `b3i.jl` as thin signal layers over it.
include("beidou/beidou.jl")
include("beidou/dnav.jl")
include("beidou/b1i.jl")
include("beidou/b3i.jl")
include("beidou/b1c.jl")
include("beidou/b2a.jl")
include("beidou/b2b.jl")
include("decoded_time.jl")
include("ephemeris.jl")
include("ephemeris_validity.jl")
include("ionosphere.jl")
include("time_parameters.jl")
include("range_accuracy.jl")
include("precompile.jl")
end
