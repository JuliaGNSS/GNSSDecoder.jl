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
- **Galileo E1B** (`GalileoE1B`) — 250 sps I/NAV nominal pages over a K=7
  rate-1/2 convolutional code plus 30×8 block interleaver. Page = 250 channel
  symbols = 1 second. Two consecutive pages (even+odd) carry one word.

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
- **L1C-D**: no fixed preamble. Sync via BCH match on the 52-symbol TOI fields
  of *two consecutive subframes*: pick TOI such that subframe N's BCH matches
  TOI=t and subframe N+1's matches TOI=t+1 (mod 400). Handles polarity
  ambiguity by accepting full-inverted matches too.
- **L5I**: the 8-bit preamble `10001011` only exists in the *decoded* bit
  domain (the FEC runs continuously across messages), so each sync attempt
  Viterbi-decodes the buffered 616-symbol window and requires the preamble at
  both ends of the decoded 308-bit window plus a clean CRC-24Q.
