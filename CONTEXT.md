# GNSSDecoder Context

Glossary of domain terms used in this package. Implementation details live in
code; this file is for naming and meaning only.

## Signals decoded

- **GPS L1 C/A** (`GPSL1CA`) — 50 bps LNAV broadcast on L1 C/A. No FEC; parity
  protection per 30-bit word. Subframe = 300 bits = 10 words = 6 seconds.
- **GPS L1C-D** (`GPSL1C_D`) — 100 sps CNAV-2 broadcast on the L1C data
  channel. Subframe == frame == 1800 symbols == 18 seconds. BCH(51,8)-coded
  TOI in the first 52 symbols; rate-1/2 LDPC over an interleaved block of 1748
  symbols for subframes 2 + 3; 24-bit CRC inside each.
- **GPS L5I** (`GPSL5I`) — 50 bps CNAV broadcast on L5 in-phase. K=7 rate-1/2
  non-systematic convolutional code (G1 = 171₈, G2 = 133₈), convolved
  *continuously* across message boundaries (no tail bits). Message = 300 bits
  = 600 channel symbols at 100 sps = 6 seconds (IS-GPS-705).
- **GPS L2C** (`GPSL2CM`) — 25 bps CNAV broadcast on the L2 CM (civil-moderate)
  code; the time-multiplexed L2 CL code is a dataless pilot. The CNAV message
  is *bit-for-bit identical* to GPS L5I's (IS-GPS-200N §30 ≡ IS-GPS-705J
  §20.3.3): same preamble, FEC, CRC-24Q, message-type layouts, π and `TOW × 6`.
  Only the signal layer differs — 25 bps → 50 sps, 300-bit message = 600
  symbols = 12 seconds — which is purely time-domain, so the symbol-domain
  decoder is the same 600-symbol message / 616-symbol sync window. Shares the
  GPS CNAV core (`src/gps/cnav.jl`) and the `GPSCNAVData` container with L5I;
  the only decode difference is that `is_sat_healthy` reports the L2 health bit
  (MT10 bit 53) instead of the L5 bit (54). `GNSSDecoderState(::GPSL2CM, prn)`
  maps here (`GPSL2CL` is the pilot); `GPSL2CMDecoderState(prn)` is the
  equivalent direct constructor.
- **Galileo E1B** (`GalileoE1B`) — 250 sps I/NAV nominal pages over a K=7
  rate-1/2 convolutional code plus 30×8 block interleaver. Page = 250 channel
  symbols = 1 second. Two consecutive pages (even+odd) carry one word.
- **Galileo E5a** (`GalileoE5aDecoderState`) — 50 sps F/NAV broadcast on the
  E5a-I component. Page = 500 channel symbols = 10 seconds = a 12-symbol sync
  pattern + 488 encoded symbols. Same K=7 rate-1/2 NSC convolutional code as
  E1B (G1 = 0o171, G2 = 0o133, G2 inverted) but a 61×8 block interleaver; one
  page decodes to 238 information bits (page type + data + CRC). Unlike I/NAV,
  each page is a complete, independently CRC-protected word (no even/odd
  stitching). Word types 1-4 carry clock/iono/health (WT1), ephemeris (WT2-3),
  and GST-UTC/GGTO + Cic/Cis (WT4); word types 5-6 carry the almanac chain.
  F/NAV rides on the E5a-I (data) component, so `GNSSDecoderState(::GalileoE5aI,
  prn)` maps here (E5a-Q is the dataless pilot); `GalileoE5aDecoderState(prn)` is
  the equivalent direct constructor.

## Frame structure terms

The same word means different things in different ICDs. Within this codebase:

- **subframe** — for GPS L1 C/A: one of five 300-bit blocks (LNAV). For L1C-D:
  a block inside a CNAV-2 frame. For L5I and Galileo: not used.
- **message** — L5I: one 300-bit CNAV unit (6 seconds), self-delimiting via
  preamble + CRC; there is no frame/subframe hierarchy.
- **frame** — L1C-D: one 18-second cycle = 1800 channel symbols. Contains
  *subframe 1* (TOI, 52 sym), *subframe 2* (CED+iono, 1200 sym), *subframe 3*
  (variable, 548 sym).
- **page** — Galileo I/NAV unit (1 sec, 250 channel symbols). Each page has
  even and odd halves; a *word* spans two consecutive pages.
- **syncro sequence** — package-internal: the smallest navigable unit the
  decoder synchronises on (subframe for GPS L1 C/A, page for Galileo, frame
  for L1C-D, message for L5I).

## TOI (Time of Interval)

L1C-D only. 9-bit counter modulo 400, broadcast in subframe 1 of every L1C-D
frame. Increments by 1 every 18-second frame. Resets at the 2-hour epoch
boundary. Carried as a BCH(51,8)-encoded 52-symbol codeword. The TOI value at
position `k` corresponds to the SV time at the *next* frame's leading edge.

## Soft symbol

A real-valued sample (`Float32`) representing one channel-symbol's worth of
post-correlation, post-Costas-demodulation data. Convention used everywhere in
this package:

- positive value ⇒ transmitted bit 0
- negative value ⇒ transmitted bit 1
- magnitude ⇒ confidence (proportional to SNR × coherent integration)

This matches AFF3CT's LLR convention. Glue code from Tracking.jl typically
takes `Float32(real(prompt))` after polarity has been resolved by the bit-sync
detector.

## IOD (Issue of Data)

A version stamp on broadcast ephemeris. The decoder uses it to confirm that
multiple subframes/pages carry a consistent ephemeris set before publishing.

- **IODC** — L1 C/A clock: 10 bits. Subframe 1 holds the full IODC; subframes
  2 and 3 carry an 8-bit IODE that must match IODC[3:10].
- **IODE_Sub_2 / IODE_Sub_3** — L1 C/A 8-bit ephemeris IODs; must equal each
  other and the lower 8 bits of IODC.
- **IODnav** — Galileo I/NAV: 10 bits, present in word types 1–4. All four
  must match for a publishable ephemeris.
- **IOD_a** — Galileo I/NAV almanac IOD, per chained almanac word (WT7–WT10).

## CRC-24Q

The CRC polynomial used by Galileo I/NAV pages, GPS L1C-D subframes 2 and 3,
and several other CNAV variants. Polynomial 0x864cfb, init 0, no reflection,
xor-out 0. Shared across signals — implemented once in this package.

## Decoder state vs cache

`GNSSDecoderState` is an immutable struct rebuilt on each `decode` step. It
references a mutable `cache` field for state that genuinely changes in place
(soft-symbol buffers, FEC decoder handles, voting tallies). New decoded
fields live in the immutable `data` / `raw_data`; "I'm still partial" state
lives in the `cache`.

## Sync mechanisms (per signal)

- **L1 C/A**: fixed 8-bit preamble `10001011` at the start of every subframe;
  TOW-continuity check across two subframes to confirm.
- **Galileo E1B**: fixed 10-bit page-sync pattern `0101100000` at the start of
  every page (in the encoded symbol stream).
- **Galileo E5a**: fixed 12-symbol F/NAV sync pattern `101101110000` at the
  start of every page; matched at both ends of the 500-symbol page window, in
  either polarity (180-degree ambiguity), exactly like E1B.
- **L1C-D**: no fixed preamble. Sync via BCH match on the 52-symbol TOI fields
  of *two consecutive subframes*: pick TOI such that subframe N's BCH matches
  TOI=t and subframe N+1's matches TOI=t+1 (mod 400). Handles polarity
  ambiguity by accepting full-inverted matches too.
- **L5I**: the 8-bit preamble `10001011` only exists in the *decoded* bit
  domain (the FEC runs continuously across messages), so each sync attempt
  Viterbi-decodes the buffered 616-symbol window and requires the preamble at
  both ends of the decoded 308-bit window plus a clean CRC-24Q.
- **L2C**: identical to L5I — the same shared GPS CNAV `try_sync`
  Viterbi-decodes the 616-symbol window and gates on preamble + CRC-24Q. (The
  symbol rate, 50 vs 100 sps, does not enter the symbol-domain sync.)
