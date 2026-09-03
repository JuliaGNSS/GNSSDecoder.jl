# Changelog

# [4.1.0](https://github.com/JuliaGNSS/GNSSDecoder.jl/compare/v4.0.0...v4.1.0) (2026-09-03)


### Features

* precompile every decoder ([fe5a514](https://github.com/JuliaGNSS/GNSSDecoder.jl/commit/fe5a5144acb4d926dc494bb75143e963e7832c7f)), closes [GNSSReceiver.jl#107](https://github.com/GNSSReceiver.jl/issues/107)

# [4.0.0](https://github.com/JuliaGNSS/GNSSDecoder.jl/compare/v3.14.0...v4.0.0) (2026-09-03)


* feat(beidou)!: gate positioning on exactly what it needs ([6d9ce8b](https://github.com/JuliaGNSS/GNSSDecoder.jl/commit/6d9ce8bbbea18ddabab0fa86a2486d020b59428d))
* refactor!: one name per quantity, raw broadcast values, corrected docs ([cd2b9c5](https://github.com/JuliaGNSS/GNSSDecoder.jl/commit/cd2b9c53999cad9b9ba594707eaef1cc7df49350)), closes [#37](https://github.com/JuliaGNSS/GNSSDecoder.jl/issues/37)


### Bug Fixes

* **beidou:** a t_op is seconds on every signal that carries one ([bdd32ae](https://github.com/JuliaGNSS/GNSSDecoder.jl/commit/bdd32aea2bb8d45c0cacb2af3ffb1bce9243549a))
* **beidou:** anchor the B2a/B2b symbol counter to the promoted frame ([57cf047](https://github.com/JuliaGNSS/GNSSDecoder.jl/commit/57cf047a888bad918d8cc39ea7e80c486afcb992))
* **beidou:** carry B1C's hour across the SOH wrap, and forget it on reset ([fa32580](https://github.com/JuliaGNSS/GNSSDecoder.jl/commit/fa325807188a4de2d999b71b03122ad030d1b21a))
* **beidou:** keep the D1/D2 SOW screen armed across vote-round clears ([d78cd1b](https://github.com/JuliaGNSS/GNSSDecoder.jl/commit/d78cd1b5c31e82b2a8d1a22835537240344507ae))
* correct four pre-existing decode defects ([985c57a](https://github.com/JuliaGNSS/GNSSDecoder.jl/commit/985c57a54f183f7804ee0672d530ec9f4db08028))
* **gpsl1c:** stamp TOI 0 on the next interval, and break the BCH tie in-band ([07461ac](https://github.com/JuliaGNSS/GNSSDecoder.jl/commit/07461ac1c5e5f7e89a8f3032319446326648db98))
* **gps:** no GPS signal broadcasts a BeiDou GGTO ([946fc4d](https://github.com/JuliaGNSS/GNSSDecoder.jl/commit/946fc4d344388b1f06b2b0a17616087c69e44d2a))
* include the defined scale offset in a broadcast time offset ([7104470](https://github.com/JuliaGNSS/GNSSDecoder.jl/commit/7104470f300af3d0bae434efef6fcee5b94bbde5))
* resolve the truncated reference week of a broadcast time offset ([78dcd78](https://github.com/JuliaGNSS/GNSSDecoder.jl/commit/78dcd78e3bf8afd73936293868e99ab7f8ad8aba))


### Features

* answer orbit class, time of week and inter-system offset per signal ([46e471d](https://github.com/JuliaGNSS/GNSSDecoder.jl/commit/46e471d840d7b69e7ae8a518c0b7def5caf3938e))
* **beidou:** decode B1C's subframe-1 BCH fields instead of matching them ([9ea8df4](https://github.com/JuliaGNSS/GNSSDecoder.jl/commit/9ea8df49c9dd9f910ae6741b9bee722800fcd3d8))
* forward the TAI epoch onto decoder states ([d24b432](https://github.com/JuliaGNSS/GNSSDecoder.jl/commit/d24b432c6ccf4397c91e8b0e6ea23aa13e4b13ee)), closes [JuliaGNSS/GNSSSignals.jl#157](https://github.com/JuliaGNSS/GNSSSignals.jl/issues/157)
* **gpsl1c:** parse the subframe-2 integrity status flag and WN_op ([78112c9](https://github.com/JuliaGNSS/GNSSDecoder.jl/commit/78112c961dc063c596686178b40dd5a67f7e3d02))


### Performance Improvements

* decode FEC blocks without allocating ([949b1c7](https://github.com/JuliaGNSS/GNSSDecoder.jl/commit/949b1c754ae7ee0cfbbb6e2875f81ed67a1ac960))


### BREAKING CHANGES

* unexports `crc24q`, `BCH_TOI_CODEWORDS`, `BCHToiSync`,
`sync_bch_toi`, `pack_hard_codeword`, `soft_to_hard_codeword`, `deinterleave!`
and `interleave!` — still reachable and documented as `GNSSDecoder.crc24q` etc.,
or via `using GNSSDecoder: crc24q`. Adds `GalileoReducedCED` to the export list,
and renames public fields on GPSL1CAData, GPSL1CAAlmanac,
GPSCNAVData, GPSL1C_DData, GalileoINAVData, GalileoE5aData, GalileoAlmanac,
GalileoE6BData and its HAS correction blocks, BeiDouDNAVData, BeiDouDNAVAlmanac,
BeiDouB1CData, BeiDouB1CBGTO, BeiDouB2aData and BeiDouB2bData; the mapping is
listed above. `GPSL1CAData.ura` and `BeiDouDNAVData.ura` are removed in favour
of the raw `URA_index` / `URAI`, so a consumer wanting metres must apply
IS-GPS-200N Table 20-I or BDS-SIS-ICD-B1I-3.0 Table 5-4 itself.
`GPSL1CAData.IODC`, `.IODE_Sub_2`, `.IODE_Sub_3`, `.sv_health`,
`GPSL1CAAlmanac.sv_health` and the two `sv_health_*_25` lists change from
`String`/`Vector{String}` holding binary literals to `Int64`/`Vector{Int64}`
holding the raw broadcast words — code comparing them against a string literal
compiles and silently evaluates false.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
* `is_decoding_completed_for_positioning` now additionally
requires `T_GD_B1Cp` and `ISC_B1Cd` on BeiDou B1C and `T_GD_B2bI` on BeiDou
B2b, and no longer requires `T_GD2` on BeiDou B1I/B3I (a relaxation: the flag
can only turn `true` earlier than before, and in practice at the same subframe). No decode timing changes — on both signals those fields are decoded by the
same constructor call as the clock the gate already required — but code that
constructs a `BeiDouB1CData` or `BeiDouB2bData` by hand must now populate them
for the flag to read `true`.
* `is_decoding_completed_for_positioning` now reads `false` on
BeiDou B1C, B2a and B2b while the broadcast `sat_type` is the reserved value 0,
where B1C and B2a previously ignored the field and B2b only required it to be
present. A satellite transmitting a reserved orbit type is no longer offered for
positioning; live satellites transmit 1, 2 or 3, so this changes nothing against
real signals, but code that constructs one of the three data containers by hand
must set `sat_type` to an orbit class for the flag to read `true`.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>

# [3.14.0](https://github.com/JuliaGNSS/GNSSDecoder.jl/compare/v3.13.0...v3.14.0) (2026-08-21)


### Bug Fixes

* 32-bit-safe wide-field extraction across the CNAV-family parsers ([9b89f75](https://github.com/JuliaGNSS/GNSSDecoder.jl/commit/9b89f759027f45d2594039a8815d95cb82750132))


### Features

* decode the five BeiDou open-service navigation messages ([f7b34ef](https://github.com/JuliaGNSS/GNSSDecoder.jl/commit/f7b34ef3c25ecd270171f3f568d516920fd91ae1))

# [3.13.0](https://github.com/JuliaGNSS/GNSSDecoder.jl/compare/v3.12.0...v3.13.0) (2026-08-20)


### Features

* support Tracking 8 ([c660c40](https://github.com/JuliaGNSS/GNSSDecoder.jl/commit/c660c40b278eb02ac95be3707f6cbc13c5b0e73f)), closes [JuliaGNSS/Tracking.jl#229](https://github.com/JuliaGNSS/Tracking.jl/issues/229)

# [3.12.0](https://github.com/JuliaGNSS/GNSSDecoder.jl/compare/v3.11.0...v3.12.0) (2026-08-18)


### Features

* support Tracking 7 ([4147c59](https://github.com/JuliaGNSS/GNSSDecoder.jl/commit/4147c59a96e02ec161ca5b3a0472f54994e0d1fc)), closes [JuliaGNSS/Tracking.jl#223](https://github.com/JuliaGNSS/Tracking.jl/issues/223)

# [3.11.0](https://github.com/JuliaGNSS/GNSSDecoder.jl/compare/v3.10.0...v3.11.0) (2026-08-13)


### Bug Fixes

* **gpsl1:** screen HOW time of week for plausibility ([19dea2d](https://github.com/JuliaGNSS/GNSSDecoder.jl/commit/19dea2d49939cede6ea382dfdcf965640f1efa94)), closes [#82](https://github.com/JuliaGNSS/GNSSDecoder.jl/issues/82)


### Features

* **gpsl1:** predict the HOW time of week from elapsed symbols ([a7b31c7](https://github.com/JuliaGNSS/GNSSDecoder.jl/commit/a7b31c7fc638199a4ebf4e594f9eb48be8303f30)), closes [#82](https://github.com/JuliaGNSS/GNSSDecoder.jl/issues/82) [#82](https://github.com/JuliaGNSS/GNSSDecoder.jl/issues/82)

# [3.10.0](https://github.com/JuliaGNSS/GNSSDecoder.jl/compare/v3.9.0...v3.10.0) (2026-08-12)


### Features

* report signal metadata for a decoder state ([7456410](https://github.com/JuliaGNSS/GNSSDecoder.jl/commit/745641095528bdc20f0d25c0f5828a7103b749d9))

# [3.9.0](https://github.com/JuliaGNSS/GNSSDecoder.jl/compare/v3.8.0...v3.9.0) (2026-08-05)


### Features

* support Tracking.jl v6 ([185004e](https://github.com/JuliaGNSS/GNSSDecoder.jl/commit/185004e87eb4677eb60ce45c2851cfb3782f6f3f))

# [3.8.0](https://github.com/JuliaGNSS/GNSSDecoder.jl/compare/v3.7.0...v3.8.0) (2026-08-03)


### Features

* support Tracking.jl v5 ([91e6936](https://github.com/JuliaGNSS/GNSSDecoder.jl/commit/91e6936ac0cf92f59759ca3eeca3e3d9ed6a0edf))

# [3.7.0](https://github.com/JuliaGNSS/GNSSDecoder.jl/compare/v3.6.0...v3.7.0) (2026-07-07)


### Features

* export is_decoding_completed_for_positioning ([e4c9d87](https://github.com/JuliaGNSS/GNSSDecoder.jl/commit/e4c9d876a6493eb673fd64f1e6a9121760b6f6af))

# [3.6.0](https://github.com/JuliaGNSS/GNSSDecoder.jl/compare/v3.5.0...v3.6.0) (2026-07-05)


### Features

* report nav-message data frequency for a decoder state ([3112804](https://github.com/JuliaGNSS/GNSSDecoder.jl/commit/3112804454d4d91f1e8fa363f566a8a1a2dc65fb))

# [3.5.0](https://github.com/JuliaGNSS/GNSSDecoder.jl/compare/v3.4.0...v3.5.0) (2026-07-05)


### Features

* introduce per-constellation data supertypes ([9e5cba6](https://github.com/JuliaGNSS/GNSSDecoder.jl/commit/9e5cba66a2ff2dfc9d5a846ab000cf4c7f1f58e9))

# [3.4.0](https://github.com/JuliaGNSS/GNSSDecoder.jl/compare/v3.3.1...v3.4.0) (2026-07-03)


### Features

* **galileo_e1b:** support the E1B BOC(1,1) signal ([ec2c597](https://github.com/JuliaGNSS/GNSSDecoder.jl/commit/ec2c597d97eedd7fb3e2f2950d7ae8c85a2e74a6))

## [3.3.1](https://github.com/JuliaGNSS/GNSSDecoder.jl/compare/v3.3.0...v3.3.1) (2026-07-03)


### Bug Fixes

* **gpscnav:** don't require T_GD for positioning readiness ([b4a2bd5](https://github.com/JuliaGNSS/GNSSDecoder.jl/commit/b4a2bd5330795ee9e23f792e9cae938679f4662f))

# [3.3.0](https://github.com/JuliaGNSS/GNSSDecoder.jl/compare/v3.2.0...v3.3.0) (2026-07-03)


### Features

* **galileo_e5a:** back-patch almanac reference epochs from a later WT5 ([1728ed4](https://github.com/JuliaGNSS/GNSSDecoder.jl/commit/1728ed45b2d75cc76d21b119ca46d4c94009b092))

# [3.2.0](https://github.com/JuliaGNSS/GNSSDecoder.jl/compare/v3.1.0...v3.2.0) (2026-07-01)


### Features

* **gpsl2c:** add GPS L2C (CNAV) decoder ([bc9a607](https://github.com/JuliaGNSS/GNSSDecoder.jl/commit/bc9a60726c83384de00f89efa3a4d731104a4626))

# [3.1.0](https://github.com/JuliaGNSS/GNSSDecoder.jl/compare/v3.0.0...v3.1.0) (2026-06-30)


### Features

* **galileo_e5a:** add Galileo E5a F/NAV decoder ([86e6c2b](https://github.com/JuliaGNSS/GNSSDecoder.jl/commit/86e6c2b70ea720efe8b8dca658bbe2fd8f509236)), closes [#83](https://github.com/JuliaGNSS/GNSSDecoder.jl/issues/83)

# [3.0.0](https://github.com/JuliaGNSS/GNSSDecoder.jl/compare/v2.1.0...v3.0.0) (2026-06-24)


* fix(l1c_d)!: correct CNAV-2 EOP ΔUT_GPS per IS-GPS-800J Table 3.5-5 ([7fb25a7](https://github.com/JuliaGNSS/GNSSDecoder.jl/commit/7fb25a7e6ce2eb972a976f15d0a3e450b0d15190))
* refactor!: unify nav-data field names across signals ([f0bd94f](https://github.com/JuliaGNSS/GNSSDecoder.jl/commit/f0bd94fc37236c02c476dd92039c99a8fcce8236)), closes [#60](https://github.com/JuliaGNSS/GNSSDecoder.jl/issues/60) [#60](https://github.com/JuliaGNSS/GNSSDecoder.jl/issues/60) [#60](https://github.com/JuliaGNSS/GNSSDecoder.jl/issues/60)


### BREAKING CHANGES

* the GPSL1C_DData fields ΔUT1/ΔUT1_dot are renamed to
ΔUT_GPS/ΔUT_GPS_dot, and the decoded ΔUT_GPS value is now twice the
previous (incorrect) value due to the 2⁻²⁴ → 2⁻²³ scale-factor fix.
Downstream code reading data.ΔUT1 must be updated.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
* renamed public struct fields on GPSL1CAData, GPSL1CAAlmanac,
GPSL5IData and GPSL1C_DData. Downstream code reading the old names must update.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>

# [2.1.0](https://github.com/JuliaGNSS/GNSSDecoder.jl/compare/v2.0.0...v2.1.0) (2026-06-23)


### Bug Fixes

* **gpsl5:** use standard GPS L5 FEC polynomials 0o171/0o133 ([6bdaaa2](https://github.com/JuliaGNSS/GNSSDecoder.jl/commit/6bdaaa2f3d9671c84a5d2eaf1a56b8cb68170e68))


### Features

* **gpsl5:** add GPS L5I (CNAV) decoder on the v2 soft-symbol API ([#29](https://github.com/JuliaGNSS/GNSSDecoder.jl/issues/29)) ([57f17ff](https://github.com/JuliaGNSS/GNSSDecoder.jl/commit/57f17ff6685032c13e0ea102d8a80b04558e6745))


### Performance Improvements

* **gpsl5:** reuse preallocated Viterbi scratch; clarify decode comments ([b717275](https://github.com/JuliaGNSS/GNSSDecoder.jl/commit/b717275b69f1ec2c247058af74bd6e4b92d54a5d))

# [2.0.0](https://github.com/JuliaGNSS/GNSSDecoder.jl/compare/v1.3.0...v2.0.0) (2026-06-23)


* feat(api)!: v2 soft-symbol decoder API + GPS L1 C/A migration ([#35](https://github.com/JuliaGNSS/GNSSDecoder.jl/issues/35)) ([1de348c](https://github.com/JuliaGNSS/GNSSDecoder.jl/commit/1de348cc1dcf6e554186b9d21c481569b2e79397)), closes [#37](https://github.com/JuliaGNSS/GNSSDecoder.jl/issues/37)


### Bug Fixes

* clamp drain_after_sync! to buffer length to survive mid-frame reset ([3c85f9c](https://github.com/JuliaGNSS/GNSSDecoder.jl/commit/3c85f9c153d25a36ba0ec05a186aee5cf782ede4))
* **test:** restore corrupted GPSL1DATA chunk 2 (issue [#35](https://github.com/JuliaGNSS/GNSSDecoder.jl/issues/35)) ([58dbcc9](https://github.com/JuliaGNSS/GNSSDecoder.jl/commit/58dbcc96a36f3cf4c0d8be7509380c53d6500ff9))


### Features

* **gpsl1c:** add GPS L1C-D decoder — TOI sync + subframe 2 ([#38](https://github.com/JuliaGNSS/GNSSDecoder.jl/issues/38)) ([3a9be1a](https://github.com/JuliaGNSS/GNSSDecoder.jl/commit/3a9be1a17c6782049fa6c6e72a292f1a43d7bf02)), closes [#39](https://github.com/JuliaGNSS/GNSSDecoder.jl/issues/39)
* **gpsl1c:** parse L1C-D subframe 3 pages ([#39](https://github.com/JuliaGNSS/GNSSDecoder.jl/issues/39)) ([d9c7233](https://github.com/JuliaGNSS/GNSSDecoder.jl/commit/d9c72337ad8ee9354ca7f58ece8438817fd9a288))
* **v2:** add LDPC alist files + shared decoder utilities (issue [#36](https://github.com/JuliaGNSS/GNSSDecoder.jl/issues/36)) ([bf27fd0](https://github.com/JuliaGNSS/GNSSDecoder.jl/commit/bf27fd0e91f232e78912e158c02ecf2925f410ad))


### BREAKING CHANGES

* decode now accepts AbstractVector{<:Real} soft symbols
instead of packed unsigned bits. Sign convention positive => bit 0,
negative => bit 1, magnitude => confidence (matches AFF3CT LLR).

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
