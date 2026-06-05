module GNSSDecoder

using DocStringExtensions,
    GNSSSignals, BitIntegers, ViterbiDecoder, CRC, Dictionaries, DataStructures

export decode,
    GPSL1CADecoderState,
    GPSL1C_DDecoderState,
    GPSL1C_DData,
    GalileoE1BDecoderState,
    is_sat_healthy,
    GNSSDecoderState,
    reset_decoder_state
# v2 shared utilities — see issue #36. Used directly by issue #37 (Galileo
# soft-input migration) and #38 (GPS L1C-D).
export crc24q, BCH_TOI_CODEWORDS, BCHToiSync, sync_bch_toi, pack_hard_codeword,
       soft_to_hard_codeword, deinterleave!, interleave!

const galCRC24 = crc(spec(24, 0x864cfb, 0x000000, false, false, 0x000000, 0xcde703))

include("gnss.jl")
include("bit_fiddling.jl")
include("gps/l1ca.jl")
include("galileo/e1b.jl")

# New deep-module utilities (issue #36). These are signal-independent and
# carry their own unit tests under `test/`. The existing `galCRC24` above
# is *not* migrated here — issue #37 will switch Galileo to consume
# `crc24q` directly, after issue #35's framework refactor lands.
include("crc.jl")
include("bch_toi.jl")
include("deinterleave.jl")

# GPS L1C-D (CNAV-2) decoder (issue #38). Included after the shared utilities
# above because it consumes `crc24q`, `sync_bch_toi`, the BCH TOI table, and
# the (de)interleaver, and after `gps/l1ca.jl` because it shares the generic
# `decode` framework hooks (`try_sync`, `decode_syncro_sequence`, …).
include("gps/l1c_d.jl")
end
