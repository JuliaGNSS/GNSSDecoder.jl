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
- **Galileo E1B / E5b** (`GalileoE1B`, `GalileoE5bI`) — 250 sps I/NAV nominal
  pages over a K=7 rate-1/2 convolutional code plus 30×8 block interleaver.
  Page part = 250 channel symbols = 1 second; two consecutive parts (even+odd)
  carry one 128-bit word, CRC-24Q over the 196 protected bits of the pair. The
  OS SIS ICD (Issue 2.2 §4.3.1) states E1-B and E5b-I use "the same page layout
  … only page sequencing is different", so they share one decoder core
  (`src/galileo/inav.jl`) and the `GalileoINAVData` container, the way L5I/L2C
  share the GPS CNAV core. The page parts differ only in fields this decoder
  never reads: E1-B's odd part spends 64 bits on OSNMA(40) + SAR(22) + Spare(2)
  where E5b-I has one "Reserved 1" field, and ends in the 8-bit SSP where E5b-I
  has "Reserved 2" — both 64 + 8, both leaving the CRC-protected prefix at 82
  bits and neither trailing field CRC-protected. Word types 16 (Reduced CED) and
  17-20 (FEC2 RS CED) are E1-B only (Table 40), so an E5b decoder simply never
  sees them. The only decode difference is that `is_sat_healthy` reports the E5b
  facet of word type 5 instead of the E1-B/C facet; `GalileoE5bDecoderState(prn)`
  is the E5b constructor (`GalileoE5bQ` is the dataless pilot).
- **Galileo E5a** (`GalileoE5aDecoderState`) — 50 sps F/NAV broadcast on the
  E5a-I component. Page = 500 channel symbols = 10 seconds = a 12-symbol sync
  pattern + 488 encoded symbols. Same K=7 rate-1/2 NSC convolutional code as
  E1B (G1 = 0o171, G2 = 0o133, G2 inverted) but a 61×8 block interleaver; one
  page decodes to 238 information bits (page type + data + CRC). Unlike I/NAV,
  each page is a complete, independently CRC-protected word (no even/odd
  stitching). F/NAV numbers its units *page types*, not word types (the I/NAV
  term): page types 1-4 carry clock/iono/health (PT1), ephemeris (PT2-3), and
  GST-UTC/GGTO + Cic/Cis (PT4); page types 5-6 carry the almanac chain.
  F/NAV rides on the E5a-I (data) component, so `GNSSDecoderState(::GalileoE5aI,
  prn)` maps here (E5a-Q is the dataless pilot); `GalileoE5aDecoderState(prn)` is
  the equivalent direct constructor.
- **Galileo E6B** (`GalileoE6B`) — 1000 sps C/NAV, the Signal-in-Space channel
  of the Galileo High Accuracy Service (HAS SIS ICD Issue 1.0). Page = 1000
  symbols = 1 second = a 16-symbol sync pattern `1011011101110000` + 984 encoded
  symbols; same K=7 rate-1/2 NSC code as I/NAV and F/NAV (the ICD says so
  explicitly) but a 123×8 block interleaver, giving 486 information bits =
  Reserved(14) + HAS Page(448) + CRC-24Q(24), the CRC over the leading 462.
  Uniquely among the signals here it carries *no ephemeris at all*: the HAS
  message is a set of corrections to another signal's broadcast navigation data,
  so `is_decoding_completed_for_positioning` is always `false` and
  `is_sat_healthy` reports the HAS *service* status instead of a satellite
  health bit. Also uniquely, one message spans *satellites*: it is cut into
  k ≤ 32 pages of 424 bits, encoded with a systematic RS(255, 32) code over
  GF(256) ("HPVRS"), and different encoded pages are handed to different
  satellites — a receiver reassembles a message from any k pages with the same
  Message ID by inverting the k×k generator submatrix their Page IDs select
  (`src/reed_solomon.jl`). E6-C is the dataless pilot, so only
  `GNSSDecoderState(::GalileoE6B, prn)` exists;
  `GalileoE6BDecoderState(prn)` is the equivalent direct constructor.
- **BeiDou B1I / B3I** (`BeiDouB1I` / `BeiDouB3I`) — the legacy BDS message,
  identical on both signals (BDS-SIS-ICD-B1I-3.0 / BDS-SIS-ICD-B3I-1.0 §5), so
  they share one decoder core (`src/beidou/dnav.jl`) and the `BeiDouDNAVData`
  container, the way L5I/L2C share the GPS CNAV core. Two formats selected by
  PRN (`is_beidou_geo`): **D1** (MEO/IGSO, PRN 6-58) at 50 bps under the NH20
  secondary code — subframes 1-3 carry the fundamental nav set, subframes 4-5
  are 24-page almanac/health/time-offset cycles (incl. the AmEpID/AmID-gated
  expanded almanac for SV 31-63); **D2** (GEO, PRN 1-5 and 59-63) at 500 bps —
  the same fundamental set spread over pages 1-10 of subframe 1 (D2 subframes
  2-5, the legacy wide-area differential/integrity service, are not decoded).
  Subframe = 300 bits = 10 × 30-bit words; word 1 = 15 raw bits + one
  BCH(15,11,1) block, words 2-10 = two bit-interleaved BCH(15,11,1) blocks.
  No IODs: `validate_data` promotes by broadcast-repetition voting keyed on
  (t_0c, t_0e) plus a SOW-vs-elapsed-symbols screen.
- **BeiDou B1C** (`BeiDouB1C_D`) — 100 sps B-CNAV1 on the B1C data component
  (BDS-SIS-ICD-B1C-1.0). Frame = 1800 symbols = 18 s: subframe 1 = 72 symbols
  (PRN as BCH(21,6) ++ SOH as BCH(51,8), concatenated, no CRC); subframes 2+3
  are 64-ary-LDPC-coded (1200 + 528 symbols) and block-interleaved together
  (36 × 48, staggered write). SF2 (600 bits, LDPC(200,100), CRC-24Q) carries
  WN/HOW + the complete ephemeris + clock + group delays behind IODC/IODE in
  one CRC-protected block; SF3 (264 bits, LDPC(88,44), CRC-24Q) is a paged
  channel (page types 1-4: iono/BDT-UTC, reduced almanacs, EOP/BGTO, midi
  almanac). The pilot `BeiDouB1C_P` is dataless, so `GNSSDecoderState(::BeiDouB1C_D,
  prn)` maps here.
- **BeiDou B2a** (`BeiDouB2aI`) — 200 sps B-CNAV2 on the B2a data component
  (BDS-SIS-ICD-B2a-1.0). Frame = 600 symbols = 3 s = 24-symbol preamble
  0xE24DE8 + one 64-ary-LDPC(96,48)-coded 288-bit message (PRN, MesType,
  SOW, 234 data bits, CRC-24Q). Message types 10/11 carry ephemeris I/II
  (paired via adjacent SOWs, |ΔSOW| = 3 s, since MT11 has no IOD), 30-34 the
  clock set (IODC) + iono/UTC/reduced-almanac/EOP/BGTO variants, 40 the midi
  almanac. B2aQ is the dataless pilot.
- **BeiDou B2b** (`BeiDouB2bI`) — 1000 sps B-CNAV3 on B2b_I
  (BDS-SIS-ICD-B2b-1.0). Frame = 1000 symbols = 1 s = 16-symbol preamble
  0xEB90 + 6 unencoded PRN symbols + 6 reserved + one
  64-ary-LDPC(162,81)-coded 486-bit message (MesType, SOW, 436 data bits,
  CRC-24Q). Message types 10 (complete ephemeris), 30 (complete clock set +
  iono/UTC), 40 (BGTO + almanacs). B-CNAV3 broadcasts no IODs — every message
  is atomic behind its CRC. B-CNAV3 is the MEO/IGSO message; the BDS-3 GEO
  satellites put the PPP-B2b service (message types 1-7, its own ICD) in the
  same framing on the same component, which this decoder frames and CRC-checks
  but does not parse — the same treatment as the legacy D2 differential
  subframes.

### The B-CNAV LDPC codes are decoded through their binary image

B1C, B2a and B2b are coded with *non-binary* LDPC codes over GF(2⁶) (the ICDs'
LDPC(200,100), (88,44), (96,48) and (162,81)). This package decodes them with
the same binary belief-propagation machinery every other LDPC signal here uses,
applied to the **binary image** of the ICD's own parity-check matrix: each
non-zero GF(2⁶) entry expands into the 6×6 GF(2) matrix of "multiply by that
element", so a bit sequence is a codeword of the image iff its symbol sequence
is a codeword of the non-binary code (`scripts/generate_beidou_alist.jl` builds
`data/bcnv*.alist` this way and cross-checks each against an independent
non-binary encoder). The code definition is therefore exact — no approximation
enters what counts as a valid frame.

The **decoder** is weaker, though, and by an amount this package has not
measured. Binary BP over the image is strictly worse than non-binary BP over
GF(2⁶): the 6×6 expansion plants length-4 cycles densely through the Tanner
graph, and the constraint binding a symbol's six bits together is discarded.
Both effects bite hardest for short, high-order codes, which is exactly what
these are, so expect a loss on the order of a dB rather than a fraction of one
— showing up as a raised frame-erasure rate at low C/N₀, since the CRC gate
turns a failed decode into a dropped frame rather than into bad data. Quantifying
it needs a reference non-binary (FFT-QSPA) decoder to compare against, which
this repository does not have.

## Field naming

Two rules, in this order:

 1. **A quantity that is the same across signals gets one package-wide name**,
    and the field documentation names the ICD symbol wherever that name departs
    from it. One spelling per quantity is what lets a consumer read a clock, an
    ephemeris or a UTC offset out of any decoder here without dispatching on the
    constellation; a per-ICD spelling would push a naming difference into every
    downstream reader as a real branch. `PositionVelocityTime.jl` evaluates the
    clock polynomial as `data.a_f0 + data.a_f1 * (t - data.t_0c) + …` for every
    decoder it supports, and that one expression only stays one expression while
    the containers agree.
 2. **A quantity that genuinely differs keeps its ICD name**, even where two
    ICDs' names then look alike or unalike. `SOW` (BDT seconds of week) is not
    `TOW`; `ΔUT1` (BDS, UT1 − UTC) is not `ΔUT_GPS` (GPS, UT1 − GPS time);
    `T_GD1`, `T_GD_B2ap` and `ISC_L5I5` are three different group delays.

The quantities where rule 1 bites — the ICDs disagree, so some ICD is departed
from and the field lists say which:

| Quantity | ICD spellings | Package |
| --- | --- | --- |
| SV clock polynomial | GPS/Galileo `a_f0`…, BDS `a_0`… (but BDS midi almanac `af0`) | `a_f0`, `a_f1`, `a_f2` |
| Clock reference epoch | GPS/BDS `t_oc`, Galileo `t_0c` | `t_0c` |
| Ephemeris reference epoch | GPS/BDS `t_oe`, Galileo `t_0e` | `t_0e` |
| Almanac reference epoch | GPS/BDS `t_oa`, Galileo `t_0a` | `t_0a` |
| Data prediction epoch | GPS `t_op`, BDS `t_op` (B2b figure prints it inside SISAIoc) | `t_op` |
| UTC polynomial | GPS `A0`/`A0-n`, Galileo `A0`, BDS `A_0UTC` | `A_0UTC`, `A_1UTC`, `A_2UTC` |
| UTC reference epoch | GPS/BDS `t_ot`, Galileo `t_0t` | `t_0t` |
| UTC reference week | GPS `WN_t`/`WN_ot`, BDS `WN_ot`, Galileo `WN_0t` | `WN_0t` |
| GNSS-GNSS offset polynomial | GPS `A0GGTO`…, BDS `A_0BGTO`… | `A_0GGTO`…, `A_0BGTO`… |
| GNSS-offset target id | GPS CNAV `GNSS ID`, GPS L1C `GGTO ID`, BDS `GNSS ID` | `GNSS_ID`; `GGTO_ID` on L1C-D only, which IRN-IS-800J-003 renamed |
| Week number | GPS LNAV `WN` (was `trans_week` here) | `WN` |
| Satellite identifier | Galileo `SVID`, GPS/BDS `PRN` | `SVID` in Galileo containers, `PRN_a` in almanac records |
| Message type of the last frame | BDS `MesType`, GPS `Message type ID` | `last_message_type` |
| Per-satellite almanac store | — (package-internal) | `almanacs` |
| BDGIM ionosphere | BDS `α_1`…`α_9` | `α_bdgim_1`…`α_bdgim_9` |

The last row is rule 2 wearing rule 1's clothes: BDS writes `α` for BDGIM and
GPS writes `α` for Klobuchar, so following both ICDs literally would give one
name to two different models. The qualifier keeps them apart.

Where the whole set is the ICD's own spelling and no note is needed: the
Keplerian elements and harmonic corrections (`M_0`, `sqrt_A`, `e`, `Ω_0`, `i_0`,
`ω`, `Ω_dot`, `i_dot`, `Δn`, `C_uc` … `C_is`), which every ICD spells alike; and
the per-signal integrity and health words, which no two constellations share
(`SISAI_oe`, `SISMAI`, `HS`, `AIF`, `DIF`, `SIF`, `URAI`, `SatH1` for BeiDou;
`E1B_SHS`, `E1B_DVS`, `BGD_E1_E5a` for Galileo; `ura_ed_index`, `ISC_L1CA` for
GPS). Two cautions inside that second group, both learned the hard way against
the documents:

  - Galileo signal health is `SHS`, not `HS`. The OS SIS ICD writes `E1BSHS` /
    `E1-BSHS`, `E5bSHS`, `E5aSHS` (Tables 30, 46, 48-51, 82, 83, 86) — *Signal
    Health Status*, S included. Only the underscore is this package's. The
    `*_DVS` half needs no such warning; the ICD spells that one as we do.
  - BeiDou writes its accuracy indices without underscores (`SISAIoe`,
    `SISAIocb`, `SISMAI`) and this package inserts one for legibility, the same
    departure as `E1B_DVS`. But `t_op` is *not* one of them: B2b §7.15 item (5)
    calls it "`top`: the time of week for data prediction", so it keeps the
    package-wide `t_op` and never a `SISAI_` prefix.

## Group delays and positioning readiness

`is_decoding_completed_for_positioning` requires a signal's own single-band
group delay **when that field rides in the block the gate already waits for**,
and excludes it otherwise. The question is never "is it important" — it always
is, at metre scale — but "does requiring it cost time to first fix".

- **Required.** BeiDou B1C's `T_GD_B1Cp` and `ISC_B1Cd` are in the same
  CRC-protected subframe 2 as the ephemeris and clock, broadcast in every
  18-second frame; BeiDou B2b's `T_GD_B2bI` is in the same MT30 as the clock;
  BeiDou B1I/B3I's `T_GD1` is in the same D1 subframe 1 / D2 page 1 as the
  clock. All cost nothing extra, and all are needed because the BDS clock term
  `a_0` is referenced to **B3I**, not to the signal being tracked (B1I §5.2.4.10,
  B1C §7.6.1, B2b §7.6) — a receiver that skips them eats a bias rather than
  losing a refinement.
- **Excluded.** BeiDou B2a's `T_GD_B2ap` / `ISC_B2ad` are in MT30 alone while
  the clock and IODC are in all of MT30-34, and the B2a ICD declines to schedule
  MT30: "the broadcast order of the B-CNAV2 message types may be dynamically
  adjusted, however Message Types 10 and 11 shall be broadcast continuously
  together" (§6.2). Gating on MT30 would gate on an unbounded interval. GPS
  CNAV's `T_GD` / `ISC_*` are excluded for the same reason.
- **Never required.** A correction for a band this decoder is not on — B1C's
  `T_GD_B2ap`, B1I/B3I's `T_GD2` (there is no B2I signal in GNSSSignals and none
  planned — JuliaGNSS/GNSSSignals.jl#156), `ISC_L5Q5` on an L1 fix. Such a field
  is still decoded and published; only the gate ignores it.

One wrinkle where the two rules meet: B1I and B3I share a `BeiDouDNAVData`, and
this gate dispatches on the data, so one check serves both signals. A B3I user
needs no group delay at all (its clock is already B3I-referenced), yet `T_GD1` is
required there too — free, since it is in the block the gate waits for anyway,
and the conservative side of the trade, because the alternative lets a B1I
consumer see a "ready" decoder with no group delay.

One further readiness condition, unrelated to group delays and on the B-CNAV
signals only: the 2-bit `sat_type` encodes 1 = GEO, 2 = IGSO, 3 = MEO and
reserves 0 (B1C/B2a/B2b Table 7-6), and it is what selects the reference the
broadcast ephemeris is expressed against — `A = A_ref + ΔA` with
`A_ref = 27 906 100 m` for MEO against `42 162 200 m` for IGSO/GEO. A satellite
broadcasting the reserved code leaves its own semi-major axis unknowable, with
14 256 100 m between the two candidates, so `is_known_sat_type`
(`beidou/beidou.jl`) keeps it out of the gate rather than let a consumer position
against a guessed orbit class. B1I and B3I are unaffected: D1/D2 NAV broadcasts
`sqrt_A` outright and carries no orbit-type field. The almanac records keep the
raw field either way — a reserved orbit type in somebody else's almanac is no
reason to discard the rest of the page.

## Frame structure terms

The same word means different things in different ICDs. Within this codebase:

- **subframe** — for GPS L1 C/A: one of five 300-bit blocks (LNAV). For L1C-D:
  a block inside a CNAV-2 frame. For BeiDou B1I/B3I: one of five 300-bit
  blocks (6 s in D1, 0.6 s in D2). For B1C: one of the three FEC blocks of a
  B-CNAV1 frame (72 + 1200 + 528 symbols). For L5I and Galileo: not used.
- **message** — L5I: one 300-bit CNAV unit (6 seconds), self-delimiting via
  preamble + CRC; there is no frame/subframe hierarchy. BeiDou B2a/B2b: the
  288-/486-bit unit inside one frame, likewise typed (MesType) and
  CRC-delimited.
- **frame** — L1C-D: one 18-second cycle = 1800 channel symbols. Contains
  *subframe 1* (TOI, 52 sym), *subframe 2* (CED+iono, 1200 sym), *subframe 3*
  (variable, 548 sym). BeiDou B1C: the analogous 18-second, 1800-symbol
  B-CNAV1 cycle. BeiDou B2a/B2b: one 600-/1000-symbol preamble-delimited
  broadcast unit (3 s / 1 s). BeiDou B1I/B3I: five subframes = 1500 bits (the
  ICD term; the decoder unit is the subframe).
- **page** — Galileo I/NAV unit (1 sec, 250 channel symbols). Each page has
  even and odd halves; a *word* spans two consecutive pages, and I/NAV types its
  content by *word type*. Galileo F/NAV inverts the relationship: a page is the
  complete CRC-protected unit and the ICD types it by *page type* (Tables 30-36),
  never "word type". For Galileo E6-B
  C/NAV: one 1000-symbol (1 s) broadcast unit — and, separately, one 424-bit
  slice of a HAS *message*, either non-encoded (`M_i`) or RS-encoded (`C_i`);
  the HAS Page ID names the latter. For BeiDou
  B1I/B3I: the Pnum-stamped cycle position of a D1 subframe 4/5 (1-24) or a
  D2 subframe (Pnum1 1-10). For B1C: a PageID-typed subframe 3 (1-4).
- **syncro sequence** — package-internal: the smallest navigable unit the
  decoder synchronises on (subframe for GPS L1 C/A and BeiDou B1I/B3I, page
  for Galileo — including one C/NAV page on E6-B — frame for L1C-D, B1C, B2a,
  and B2b, message for L5I).

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
  must match for a publishable ephemeris. Galileo HAS corrections name it as
  their `IOD_ref` (8 bits for the GPS IODE/IODC they correct instead).
- **IOD Set ID** — Galileo HAS: a 5-bit counter over the *set* of reference IODs
  of every satellite in a mask. It changes when any one satellite's reference IOD
  does; it is not itself an issue of data, but the key that ties a correction
  block to the ephemeris issue it corrects (HAS SIS ICD §7.6).
- **Mask ID** — Galileo HAS: a 5-bit counter over the set of corrected
  satellites, signals and reference navigation messages. Correction blocks
  reference a mask by this ID, and cannot be parsed at all without it — every
  block's *length* is derived from the mask.
- **IOD_a** — Galileo I/NAV almanac IOD, per chained almanac word (WT7–WT10).
- **BeiDou**: B-CNAV1 carries IODC/IODE inside the single CRC-protected
  subframe 2, so no cross-block stitching is needed; B-CNAV2 gates the MT10/11
  ephemeris pair by adjacent SOWs and the clock via IODC/IODE (`IODE == IODC &
  0xff`, §7.4.3); B-CNAV3 broadcasts no IODs (atomic messages); the legacy
  D1/D2 message has only AODC/AODE *ages*, so consistency comes from
  broadcast-repetition voting keyed on (t_0c, t_0e).

## CRC-24Q

The CRC polynomial used by Galileo I/NAV page pairs and E6-B C/NAV pages, GPS L1C-D subframes 2 and 3,
the BeiDou B-CNAV blocks (B1C subframes 2 and 3, every B2a and B2b message),
and several other CNAV variants. Polynomial 0x864cfb, init 0, no reflection,
xor-out 0. Shared across signals — implemented once in this package.

## Signal metadata

Facts about the *signal* a decoder demodulates — its name and ids, band,
constellation, navigation-message symbol rate, time scale — as opposed to the
navigation *data* it decodes. All of it lives in GNSSSignals and none of it is
restated here: each signal file states the constants → signal mapping
(`get_signal_type`) and `src/gnss.jl` forwards GNSSSignals' accessors for a
`GNSSDecoderState` through it.

The mapping names the data-bearing component of a signal pair, since that is
what carries the navigation message (`GPSL2CM`, not the `GPSL2CL` pilot;
`GalileoE5aI`, not `GalileoE5aQ`), and reports the signal an approximation
approximates (the E1B BOC(1,1) decoder is a `GalileoE1B`). It is keyed on the
constants type because that, not the data type, is what tells apart decoders
sharing a data container: L5-I and L2C-M both decode into a `GPSCNAVData`.

## Time system

The atomic (leap-second-free) scale a constellation counts its broadcast time
in: GPS Time (`GPST`) for GPS, Galileo System Time (`GST`) for Galileo. Every
decoded week number and time of week is in this scale, so turning a decoded
WN/TOW pair into an absolute instant needs the scale's epoch
(`get_system_start_time`) and its offset from TAI (`get_tai_offset`) — both
available from the decoder state, both sourced from GNSSSignals. The Galileo
epoch is not a UTC midnight: it is 13 s before it (`1999-08-21T23:59:47`).

## Decoder state vs cache

`GNSSDecoderState` is an immutable struct rebuilt on each `decode` step. It
references a mutable `cache` field for state that genuinely changes in place
(soft-symbol buffers, FEC decoder handles, voting tallies). New decoded
fields live in the immutable `data` / `raw_data`; "I'm still partial" state
lives in the `cache`.

## Sync mechanisms (per signal)

- **L1 C/A**: fixed 8-bit preamble `10001011` at the start of every subframe;
  TOW-continuity check across two subframes to confirm.
- **Galileo E1B / E5b**: fixed 10-bit page-sync pattern `0101100000` at the
  start of every page part (in the encoded symbol stream).
- **Galileo E6B**: fixed 16-symbol C/NAV sync pattern `1011011101110000` at the
  start of every page, matched at both ends of the 1016-symbol window in either
  polarity, then hardened by Viterbi-decoding the page and requiring CRC-24Q —
  the RS erasure decoder downstream trusts every collected page absolutely, so
  a page must be error-free before it enters the store.
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
- **BeiDou B1I/B3I**: fixed 11-bit preamble `11100010010` at the start of
  every 300-bit subframe (unencoded bits 1-11 of word 1), matched at both
  ends of the window in either polarity, plus a SOW-vs-elapsed-symbols screen.
- **BeiDou B1C**: no fixed preamble; like L1C-D, sync is a BCH match on the
  72-symbol subframe 1 of *two consecutive frames* — the own-PRN BCH(21,6)
  codeword plus a BCH(51,8) SOH codeword at frame N and SOH+1 (mod 200) at
  frame N+1, in either polarity.
- **BeiDou B2a**: fixed 24-symbol preamble `0xE24DE8` at both ends of the
  600-symbol frame window, either polarity; the LDPC + CRC-24Q + own-PRN
  gates then confirm inside `decode_syncro_sequence`.
- **BeiDou B2b**: fixed 16-symbol preamble `0xEB90` at both ends of the
  1000-symbol frame window, either polarity, strengthened by the 6 unencoded
  own-PRN symbols and the LDPC + CRC-24Q gate before any state is updated.
