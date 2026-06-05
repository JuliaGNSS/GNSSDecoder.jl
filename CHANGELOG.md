# Changelog

# [2.0.0](https://github.com/JuliaGNSS/GNSSDecoder.jl/compare/v1.3.0...v2.0.0) (2026-05-22)

### BREAKING CHANGES

* **api:** decoder now consumes `AbstractVector{<:Real}` soft symbols
  instead of a packed-bit `Unsigned`. `Float32` is canonical; sign carries
  the bit decision (positive ⇒ bit 0, negative ⇒ bit 1) and magnitude
  carries confidence (AFF3CT-LLR convention). Closes #35.
* **types:** `GPSL1Data` → `GPSL1CAData`, `GPSL1Constants` →
  `GPSL1CAConstants`, `GPSL1DecoderState` → `GPSL1CADecoderState`. The
  corresponding `GPSL1Almanac` / `VotedGPSL1Data` / `GPSL1Cache` renamed
  to the `*CA*` variants in step.
* **galileo:** Galileo E1B is now decoded from `Float32` soft symbols
  end-to-end. The K=7 NSC convolutional FEC is undone with AFF3CT.jl's
  `ConvViterbiDecoder` (`poly = [0o171, 0o133]`, K=114, N=240) instead of
  ViterbiDecoder.jl, and the `ViterbiDecoder` and `CRC` dependencies are
  dropped. Closes #37.

### Features

* **api:** new `decode(state, soft_symbols::AbstractVector{<:Real}, n)`
  entry point on `GNSSDecoderState`. Internal buffer is a
  `CircularDeque{Float32}` (DataStructures.jl) of length
  `syncro_sequence_length + preamble_length` per signal (308 for GPS
  L1 C/A; 260 for Galileo E1B).
* **layout:** source tree reorganised — `src/gpsl1.jl` → `src/gps/l1ca.jl`,
  `src/galileo_e1b.jl` → `src/galileo/e1b.jl`.
* **galileo:** Galileo E1B's `decode_syncro_sequence` reads the 240
  polarity-corrected `Float32` LLR soft symbols directly from the deque
  (sync hook hard-slices the tail only for the 10-bit page-sync preamble
  bit-pattern check). The 30×8 block deinterleave and "invert every second
  bit" steps now operate on soft symbols (inversion = negation), and the
  page parser, even/odd caching, and WT7→WT10 almanac chaining are
  unchanged. The inline `galCRC24` const is removed in favour of the shared
  `crc24q` (issue #36).
* **gpsl1c:** new GPS L1C-D (CNAV-2) decoder — TOI frame sync plus full
  subframe-2 parsing. `GPSL1C_DDecoderState(prn)` wires up a 1852-symbol
  soft buffer, the BCH(51,8) TOI codeword table for frame sync
  (`sync_bch_toi`), and two lazily-built AFF3CT LDPC belief-propagation
  decoders (SF2 K=600/N=1200, SF3 K=274/N=548) from the committed
  `data/cnv2_sf*.alist` matrices. Decodes subframe 2 (clock, ephemeris and
  accuracy parameters per IS-GPS-800G Figure 3.5-1 / Table 3.5-1) into the
  new `GPSL1C_DData`, with 38×46 block deinterleaving, CRC-24Q validation,
  and 180° polarity resolution. Subframe 3 is LDPC-decoded and CRC-checked
  but only counted as a received page (`num_sf3_pages_received`); per-page
  field parsing is deferred to #39. Closes #38.
* **gpsl1c:** parse GPS L1C-D (CNAV-2) subframe-3 pages. After the SF3 CRC
  passes, `decode_subframe3` dispatches on the 6-bit page number (bits 9-14;
  bits 1-8 are the transmitting PRN) and merges the parsed fields into
  `GPSL1C_DData` immutably, mirroring the subframe-2 style. Implements the
  IRN-IS-800J layouts: page 1 — UTC parameters + Klobuchar ionospheric
  coefficients + ISC (L1C/A, L2C, L5I5, L5Q5); page 2 — GGTO + Earth
  orientation parameters; page 3 — reduced almanac (six 33-bit per-SV packets);
  page 4 — Midi almanac; page 5 — clock + ephemeris differential correction;
  page 6 — 29-character ASCII text. Almanacs and differential corrections are
  stored as `Dictionaries.Dictionary` keyed by SV PRN (`reduced_almanacs`,
  `midi_almanacs`, `differential_corrections`) via the new
  `GPSL1C_DReducedAlmanac`, `GPSL1C_DMidiAlmanac`, and
  `GPSL1C_DDifferentialCorrection` structs; both almanac types are
  single-page-per-SV (no Galileo-style IOD chaining). Unknown/reserved pages
  (7 SV-config, 8 ISM, …) are silently ignored while still counted in
  `num_sf3_pages_received`. Semi-circle quantities are converted to radians on
  decode. Pre-IRN-J page-1 recordings (no ISC fields) are out of scope.
  Closes #39.
* **gpsl1c:** validate the L1C-D pipeline end-to-end against a Spirent
  GSS7000 recording, fixing two bit-exact bugs that synthetic round-trip
  tests could not catch: (1) the BCH(51,8) TOI LFSR emitted the complement
  of the register LSB — the corrected table matches the recording's
  subframe-1 symbols bit-for-bit (PocketSDR's `LFSR` emits `CHIP[R & 1]`
  with `CHIP = (-1, 1)`, i.e. the LSB itself); (2) the CNAV-2 LDPC
  codewords are systematic `[info | parity]` per IS-GPS-800G §3.2.3.3, but
  Aff3ct's alist loader auto-derives info-bit positions that select the
  parity columns — `GPSL1C_DCache` now forces `info_bits_pos = 0:K-1`
  before building the BP decoders. The gated Spirent fixture test now
  parses the GSS-CNAVDATA container (16-byte header, 225-byte post-FEC
  symbol blocks, per-epoch satellite-channel round-robin) and asserts the
  decoded SF2/SF3 fields against Spirent's own pre-FEC field dump in both
  polarities; the synthetic CI test encodes systematically over GF(2) from
  the alist instead of via Aff3ct's LU-based encoder.
* **gpsl1c:** commit two Spirent-derived recording fixtures so the
  end-to-end L1C-D test runs automatically in CI. The post-FEC channel
  symbols of PRN 1 and PRN 2 (69 CNAV-2 frames each) are stored as
  packed-bit files under `test/data/` (15,525 bytes per satellite — the
  ±1 hard symbols carry one bit each, so packed bits beat `Float32`
  literals 32:1) and unpacked to ±1 `Float32` soft symbols at the test
  boundary. PRN 1 is golden-checked field-by-field against Spirent's
  pre-FEC dump (including exact reduced-/midi-almanac entries) in both
  polarities; PRN 2 cross-checks satellite-specific ephemeris against
  the midi-almanac entry decoded from PRN 1's stream and requires the
  constellation-wide subframe-3 broadcast to decode identically from
  both satellites. The previous env-var-gated full-recording test
  (`GPS_L1C_D_FIXTURE_DIR`) is removed: it never ran in CI, the fixture
  provenance byte-compare against the source recording was verified
  once at extraction time, and all 31 channels share the decode path
  the two committed PRNs already cover.
* **cleanup:** remove the dead `String`-based `deinterleave` and
  `invert_every_second_bit` helpers left over from the hard-bit era —
  Galileo E1B now runs the shared `Float32` deinterleaver and negation;
  the Galileo-ICD example test drives the same soft path the decoder
  uses. Documentation (README, manual index) updated to list GPS L1C-D
  and the soft-symbol input convention.

### Internal

* `GNSSDecoderState` drops `raw_buffer`, `buffer`, and `num_bits_buffered`
  from its top-level fields. The packed-bit `UInt320` / `UInt288` buffers
  live transiently inside the per-signal cache and are populated at sync
  time by `pack_buffer_into_cache!`. The packed buffer's layout
  (oldest bit at MSB, newest at LSB) matches the v1 `raw_buffer`, so the
  existing word-extraction / parity-check / subframe-parsing helpers are
  reused unchanged.
* Galileo E1B's `decode_syncro_sequence` reads the hard-sliced packed
  buffer from `state.cache.complemented_buffer[]` instead of
  `state.buffer`; everything downstream of that read (deinterleave,
  Viterbi, page parser) is byte-for-byte unchanged.

# [1.3.0](https://github.com/JuliaGNSS/GNSSDecoder.jl/compare/v1.2.0...v1.3.0) (2026-06-03)


### Features

* **deps:** bump GNSSSignals to v2.2 ([78e8597](https://github.com/JuliaGNSS/GNSSDecoder.jl/commit/78e8597d0f258bd8f4af94ad21be812d8fbf6d61))

# [1.2.0](https://github.com/JuliaGNSS/GNSSDecoder.jl/compare/v1.1.0...v1.2.0) (2026-05-17)


### Bug Fixes

* discard Galileo E1B almanac data with SVID = 0 ([b2f5554](https://github.com/JuliaGNSS/GNSSDecoder.jl/commit/b2f55540ec1cdead5be9e61e2eab4cf73d8ac0be))
* store Galileo E1B IOD_a, WN_a, and t_0a for each satellite almanac separately ([1b006e8](https://github.com/JuliaGNSS/GNSSDecoder.jl/commit/1b006e891ea284b6c6a91005078da1adf76949b1))


### Features

* **galileo:** decode remaining I/NAV bits from word types 3-10 and 16 ([922be0e](https://github.com/JuliaGNSS/GNSSDecoder.jl/commit/922be0e06ee7e4ff44a5fe46fa55aa3078c9ad54))

# [1.1.0](https://github.com/JuliaGNSS/GNSSDecoder.jl/compare/v1.0.0...v1.1.0) (2026-05-08)


### Features

* implement GPS L1 subframe 4 and 5 decoding ([3d77cc2](https://github.com/JuliaGNSS/GNSSDecoder.jl/commit/3d77cc2b0d3b9dd05014b5d5018063f65a08eeef))

# [0.2.0](https://github.com/JuliaGNSS/GNSSDecoder.jl/compare/v0.1.5...v0.2.0) (2026-01-04)


### Bug Fixes

* handling of new IODCs ([bd26ac8](https://github.com/JuliaGNSS/GNSSDecoder.jl/commit/bd26ac8f808eb6b6852028c9474c5ff72889af0c))


### Features

* voting based gps ca data validation ([47e75dd](https://github.com/JuliaGNSS/GNSSDecoder.jl/commit/47e75dd61e452904e25a2006bce54f4f83aba687))

## [0.1.5](https://github.com/JuliaGNSS/GNSSDecoder.jl/compare/v0.1.4...v0.1.5) (2025-12-31)


### Bug Fixes

* num_bits_after_valid_syncro_sequence computation after missed pages ([b5dd500](https://github.com/JuliaGNSS/GNSSDecoder.jl/commit/b5dd500b4309f10740c2e0ca96a181481851700d))
