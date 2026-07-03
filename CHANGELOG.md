# Changelog

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
