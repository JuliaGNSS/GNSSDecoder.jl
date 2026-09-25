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
 * Galileo E1B (I/NAV), including the BOC(1,1) approximation
 * Galileo E5b (I/NAV)
 * Galileo E5a (F/NAV)
 * Galileo E6B (C/NAV — the High Accuracy Service corrections)
 * BeiDou B1I (D1/D2 NAV)
 * BeiDou B3I (D1/D2 NAV)
 * BeiDou B1C (B-CNAV1)
 * BeiDou B2a (B-CNAV2)
 * BeiDou B2b (B-CNAV3)

## Installation

```julia
julia> ]
pkg> add GNSSDecoder
```

## Usage

```julia
using GNSSDecoder

state = GPSL1CADecoderState(25)   # every buffer is allocated here, once
for chunk in soft_symbol_chunks   # e.g. `get_soft_bits` from Tracking.jl
    state = decode!(state, chunk, length(chunk))
end
is_decoding_completed_for_positioning(state) && is_sat_healthy(state)
```

`decode!` allocates nothing: it **overwrites** the buffers the decoder state
was constructed with (soft-symbol buffer, FEC scratch, the almanac stores and
other containers behind `state.raw_data` and `state.data`). Always continue
with the state it returns and do not keep using the old one — take a
`copy(state)` if you need a snapshot.

## Standalone executables (`juliac --trim`)

On Julia 1.12 and later every decoder compiles under
[JuliaC](https://github.com/JuliaLang/JuliaC.jl)'s `--trim=safe`, so you can
build a small standalone binary around it. `test/trim/main.jl` is a working
example, and CI builds it on every change:

```sh
julia -e 'using Pkg; Pkg.Apps.add("JuliaC")'
julia --project=test/trim -e 'using Pkg; Pkg.instantiate()'
juliac --output-exe gnssdecoder_trim --bundle build --trim=safe --experimental \
    --project test/trim test/trim/main.jl
build/bin/gnssdecoder_trim
```

The LDPC parity-check matrices (GPS L1C-D, BeiDou B1C/B2a/B2b) are embedded
in the package image, so the bundled binary does not need the source tree.

## Documentation

For usage examples and API reference, see the [documentation](https://JuliaGNSS.github.io/GNSSDecoder.jl/stable/).

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE.md) file for details.
