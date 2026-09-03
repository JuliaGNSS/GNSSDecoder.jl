module GNSSDecoder

using DocStringExtensions, GNSSSignals, BitIntegers, Dictionaries, DataStructures
using GNSSSignals: Hz, s, ustrip
import Aff3ct

export decode,
    GPSL1CADecoderState,
    GPSL1CAData,
    GPSL1CAAlmanac,
    GPSL1C_DDecoderState,
    GPSL1C_DData,
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
    GalileoE5aData,
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
    BeiDouB1CData,
    BeiDouB1CBGTO,
    BeiDouB2aDecoderState,
    BeiDouB2aData,
    BeiDouB2bDecoderState,
    BeiDouB2bData,
    BeiDouMidiAlmanac,
    BeiDouReducedAlmanac,
    is_sat_healthy,
    is_decoding_completed_for_positioning,
    get_signal_type,
    OrbitClass,
    get_orbit_class,
    get_time_of_week,
    GNSSTimeOffset,
    get_time_offset,
    GNSSDecoderState,
    reset_decoder_state

# Shared decode framework and cross-signal accessors, then the generic
# bit-extraction helpers.
include("gnss.jl")
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
end
