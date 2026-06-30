# GNSSDecoder.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaGNSS.github.io/GNSSDecoder.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaGNSS.github.io/GNSSDecoder.jl/dev/)
[![Build Status](https://github.com/JuliaGNSS/GNSSDecoder.jl/actions/workflows/ci.yml/badge.svg?branch=master)](https://github.com/JuliaGNSS/GNSSDecoder.jl/actions/workflows/ci.yml?query=branch%3Amaster)
[![codecov](https://codecov.io/gh/JuliaGNSS/GNSSDecoder.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaGNSS/GNSSDecoder.jl)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Decodes various GNSS satellite signals from soft symbols (`Float32`: sign =
bit decision, magnitude = confidence) as produced by a tracking loop such as
`Tracking.jl` (feed `get_soft_bits`). By convention positive ⇒ bit 0 and
negative ⇒ bit 1, but this is only a convention: the decoder resolves the 180°
polarity ambiguity itself, so the opposite sign decodes the same data, and the
magnitudes need not be normalized.

Currently implemented:
 * GPS L1 C/A (LNAV)
 * GPS L1C-D (CNAV-2)
 * GPS L2C (CNAV)
 * GPS L5I (CNAV)
 * Galileo E1B (I/NAV)
 * Galileo E5a (F/NAV)

## Installation

```julia
julia> ]
pkg> add GNSSDecoder
```

## Documentation

For usage examples and API reference, see the [documentation](https://JuliaGNSS.github.io/GNSSDecoder.jl/stable/).

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE.md) file for details.
