# Changelog

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
