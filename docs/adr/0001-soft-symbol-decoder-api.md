# ADR 0001 — Soft-symbol decoder API for v2

Status: **Accepted** (2026-05-22)

## Context

Until v1, `GNSSDecoder` exposed `decode(state, bits::Unsigned, num_bits::Int)`
— each call passed hard-decision bits packed into an unsigned integer. This
was natural when only GPS L1 C/A (no FEC) and Galileo E1B (decoded inside via
hard-decision Viterbi) needed support.

Adding GPS L1C-D forces the question. L1C-D uses rate-1/2 LDPC for subframes 2
and 3. LDPC decoding gives a meaningful coding gain (~2 dB) only with
soft-decision inputs; hard-decision LDPC throws the whole reason for picking
LDPC over a simpler code in the air. The same argument applies to Galileo E1B
(soft Viterbi gain ≈ 2 dB over hard) and GPS L5I (same convolutional code as
Galileo E1B's). Once *one* signal needs soft inputs and Tracking.jl v2 begins
exposing soft prompts via `get_filtered_prompts`, leaving the other decoders
on a hard-bit API costs us the soft information at the boundary on every
single signal — the worst place to throw it away.

## Decision

Bump to v2.0 and replace the hard-bit API with a soft-input API across every
decoder:

```julia
decode(state, soft_symbols::AbstractVector{<:Real}, num_symbols::Int)
```

Conventions:

- `Float32` everywhere — matches AFF3CT's LLR type with no per-call
  conversion.
- Positive value ⇒ transmitted bit 0, negative ⇒ bit 1, magnitude ⇒
  confidence. See [[soft-symbol convention in CONTEXT.md]].
- Internal storage: `CircularDeque{Float32}` (DataStructures.jl).
  Per-signal buffer length: 308 (L1 C/A: 300+8), 260 (E1B: 250+10), 1852
  (L1C-D: 1800+52), TBD for L5I.

All FEC moves to AFF3CT.jl:

- LDPC BP decoder for L1C-D subframes 2 and 3 (separate `.alist` parity-check
  matrices shipped under `data/`).
- Soft Viterbi for Galileo E1B (and later L5I), replacing the
  `ViterbiDecoder.jl` dependency.

LDPC `.alist` files load lazily inside each `GPSL1C_DDecoderState(prn)`
constructor — no `__init__` hook (juliac-friendly).

## Alternatives considered

1. **Two parallel APIs** (hard for L1 C/A + Galileo, soft only for L1C-D).
   Rejected: doubles the API surface, leaves Galileo on hard-decision Viterbi
   for no good reason, asks callers to remember which signal uses which API.

2. **Soft for L1C-D only; migrate the rest in a follow-up release.** Rejected
   in favour of a single v2 bump that coincides with the GNSSSignals.jl v2 and
   Tracking.jl v2 bumps the rest of the JuliaGNSS stack is taking. Migrating
   each decoder in its own minor bump churns downstream users twice.

3. **Stay on hard bits, do hard LDPC for L1C-D.** Rejected: the whole point of
   LDPC over a simpler code is the soft-decision gain; on Spirent recordings
   the BER is low enough that hard-decision LDPC would work, but real-receiver
   performance would be substantially worse.

## Consequences

Positive

- Single API across all signals; LLR convention shared across the FEC stack.
- L1C-D, Galileo E1B and L5I all get the soft-decision coding gain.
- One FEC library (AFF3CT) replaces ViterbiDecoder.jl; LDPC is added for free.
- Matches how Tracking.jl v2 already exposes prompts.

Negative

- Breaking API change. Every downstream caller (GNSSReceiver.jl, pvt.jl,
  bench scripts) updates once at the v2 bump.
- Tracking.jl-to-decoder glue must now project complex prompts onto the data
  axis: `Float32(real(prompt))` plus polarity from bit-sync. Documented in
  [[CONTEXT.md soft symbol section]].
- AFF3CT is a heavier dependency than ViterbiDecoder.jl (binary blob, C++
  bindings). Mitigated by the fact that L1C-D LDPC requires it anyway.

## References

- IS-GPS-800G §3.2.3.3–3.2.3.5 (LDPC, interleaver)
- Galileo OS SIS ICD Issue 2.2 §4.3.3 (E1B FEC)
- IS-GPS-705 §3.5 (L5 CNAV)
- AFF3CT.jl: `LDPCBPDecoder`, `ViterbiDecoder`
- PocketSDR `python/sdr_nav.py` — reference implementation of L1C-D
  multi-subframe BCH sync
