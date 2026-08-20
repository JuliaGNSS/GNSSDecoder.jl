module GNSSDecoder

using DocStringExtensions, GNSSSignals, BitIntegers, Dictionaries, DataStructures
# `Hz` is Unitful's, re-exported through GNSSSignals: the only decoder that
# states a rate itself is BeiDou D1/D2, whose symbol rate is a property of the
# satellite rather than of the signal type (see `beidou/dnav.jl`).
using GNSSSignals: Hz
import Aff3ct

export decode,
    GPSL1CADecoderState,
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
    GNSSDecoderState,
    reset_decoder_state
# v2 shared utilities — see issue #36. Used directly by issue #37 (Galileo
# soft-input migration) and #38 (GPS L1C-D).
export crc24q,
    BCH_TOI_CODEWORDS,
    BCHToiSync,
    sync_bch_toi,
    pack_hard_codeword,
    soft_to_hard_codeword,
    deinterleave!,
    interleave!

include("gnss.jl")
include("bit_fiddling.jl")

# Signal-independent deep-module utilities (issue #36). These carry their own
# unit tests under `test/`. They are included before the per-signal decoders
# because Galileo E1B (issue #37) consumes `crc24q` and `deinterleave!`
# directly.
include("crc.jl")
include("bch_toi.jl")
include("deinterleave.jl")

# Shared LDPC decode helper (Aff3ct BP decode → CRC-24Q gate → MSB-first
# packing), consumed by GPS L1C-D and the BeiDou B-CNAV decoders. Included
# after `crc.jl` (it consumes `crc24q`) and before every decoder that uses it.
include("ldpc.jl")

include("gps/l1ca.jl")

# Definitions shared across Galileo signals (the `SignalHealth` /
# `DataValidityStatus` enums, the `GalileoAlmanac` record, the common K=7 NSC
# `galileo_viterbi` FEC primitive). Included before the per-signal Galileo
# decoders, which all consume it — analogous to how `gnss.jl` precedes every
# signal. Itself depends only on the earlier shared utilities (`deinterleave`,
# `bit_fiddling`) and `Aff3ct`.
include("galileo/galileo.jl")
include("galileo/e1b.jl")

# Galileo E5a (F/NAV) decoder. Consumes the shared Galileo definitions above
# (`galileo/galileo.jl`) plus `crc24q`, `deinterleave`, and the generic `decode`
# framework hooks.
include("galileo/e5a.jl")

# GPS L1C-D (CNAV-2) decoder (issue #38). Included after the shared utilities
# above because it consumes `crc24q`, `sync_bch_toi`, the BCH TOI table, and
# the (de)interleaver, and after `gps/l1ca.jl` because it shares the generic
# `decode` framework hooks (`try_sync`, `decode_syncro_sequence`, …).
include("gps/l1c_d.jl")

# Shared GPS CNAV core (the 300-bit message broadcast identically on GPS L5I,
# IS-GPS-705J §20.3, and GPS L2C, IS-GPS-200N §30): FEC, window-Viterbi sync,
# CRC-24Q, and all per-message-type parsing, plus the shared `GPSCNAVData`
# container. Consumes only shared utilities included above — `crc24q` and the
# `UInt320` / `_merge_keyed` primitives in `gnss.jl` — so it is independent of
# the sibling signal files.
include("gps/cnav.jl")

# GPS L5I and GPS L2C signal layers — each a thin wrapper over `gps/cnav.jl`
# (their own `*Constants` type, decoder-state constructor, and `is_sat_healthy`
# health-bit selection). Included after the shared core they consume.
include("gps/l5.jl")
include("gps/l2c.jl")

# Definitions shared across BeiDou signals (BDCS constants, the GEO PRN
# partition selecting D1 vs D2 on B1I/B3I). Included before the per-signal
# BeiDou decoders, which all consume it — analogous to `galileo/galileo.jl`.
include("beidou/beidou.jl")

# Shared BeiDou legacy NAV core (the D1/D2 message broadcast identically on
# B1I, BDS-SIS-ICD-B1I-3.0 §5, and B3I, BDS-SIS-ICD-B3I-1.0 §5): BCH(15,11,1)
# word decoding, subframe/page parsing, SOW screening, data voting, and the
# shared `BeiDouDNAVData` container. Consumes `beidou/beidou.jl` and the
# shared GPS primitive `increment_voting` in `gps/l1ca.jl`, so it is included
# after both.
include("beidou/dnav.jl")

# BeiDou B1I and B3I signal layers — each a thin wrapper over `beidou/dnav.jl`
# (their own `*Constants` tag, decoder-state constructor, and `is_sat_healthy`
# selection). Included after the shared core they consume.
include("beidou/b1i.jl")
include("beidou/b3i.jl")

# BeiDou B1C (B-CNAV1) decoder. Included after `gps/l1c_d.jl` (it reuses the
# `UInt600` packed-word type and mirrors the BCH two-subframe sync pattern)
# and after the shared LDPC helper (`ldpc_decode_word`, `load_ldpc_decoder`).
include("beidou/b1c.jl")

# BeiDou B2a (B-CNAV2) decoder. Consumes the shared BeiDou definitions above
# plus the shared LDPC pipeline (`src/ldpc.jl`, `data/bcnv2.alist`) and the
# generic decode framework hooks.
include("beidou/b2a.jl")

# BeiDou B2b (B-CNAV3) decoder. Consumes the shared BeiDou definitions above
# plus the shared LDPC helper (`ldpc.jl`) and the `data/bcnv3.alist` binary
# image of the ICD's 64-ary LDPC(162, 81) code.
include("beidou/b2b.jl")
end
