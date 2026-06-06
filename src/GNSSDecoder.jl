module GNSSDecoder

using DocStringExtensions,
    GNSSSignals, BitIntegers, Dictionaries, DataStructures
import Aff3ct

export decode,
    GPSL1CADecoderState,
    GPSL1C_DDecoderState,
    GPSL1C_DData,
    GPSL1C_DReducedAlmanac,
    GPSL1C_DMidiAlmanac,
    GPSL1C_DDifferentialCorrection,
    GPSL5IDecoderState,
    GPSL5IData,
    GPSL5IReducedAlmanac,
    GPSL5IMidiAlmanac,
    GPSL5IClockDifferentialCorrection,
    GPSL5IEphemerisDifferentialCorrection,
    GPSL5IIntegritySupportMessage,
    GalileoE1BDecoderState,
    is_sat_healthy,
    GNSSDecoderState,
    reset_decoder_state
# v2 shared utilities — see issue #36. Used directly by issue #37 (Galileo
# soft-input migration) and #38 (GPS L1C-D).
export crc24q, BCH_TOI_CODEWORDS, BCHToiSync, sync_bch_toi, pack_hard_codeword,
       soft_to_hard_codeword, deinterleave!, interleave!

include("gnss.jl")
include("bit_fiddling.jl")

# Signal-independent deep-module utilities (issue #36). These carry their own
# unit tests under `test/`. They are included before the per-signal decoders
# because Galileo E1B (issue #37) consumes `crc24q` and `deinterleave!`
# directly.
include("crc.jl")
include("bch_toi.jl")
include("deinterleave.jl")

include("gps/l1ca.jl")
include("galileo/e1b.jl")

# GPS L1C-D (CNAV-2) decoder (issue #38). Included after the shared utilities
# above because it consumes `crc24q`, `sync_bch_toi`, the BCH TOI table, and
# the (de)interleaver, and after `gps/l1ca.jl` because it shares the generic
# `decode` framework hooks (`try_sync`, `decode_syncro_sequence`, …).
include("gps/l1c_d.jl")

# GPS L5I (CNAV) decoder. Included after `gps/l1ca.jl` (reuses its `UInt320`
# packed-word type) and after `gps/l1c_d.jl` (reuses the `_deque_slice` and
# `_merge_keyed` helpers), and consumes `crc24q` from the shared utilities
# above.
include("gps/l5.jl")
end
