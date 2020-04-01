<<<<<<< README.md
# GNSSDecoder.jl

Decodes GPSL1 satellite signals and computes the Satellite position in WGS84.
## Features
* Decoding of Satellite Signals
* Satellite Position
* Pseudorange
## Usage

#### Install: 
```julia
julia> ]
(@v1.4) pkg> add GNSSDecoder


```
#### Decoding:
```julia
using GNSSDecoder, Tracking

decoder = GNSSDecoderState()
carrier_doppler = 1000Hz
code_phase = 50
sample_frequency = 2.5e6Hz
prn = 1
state = TrackingState(GPSL1, carrier_doppler, code_phase)
for x in 1:1000
    # get signal
    track_res = track(signal, state, prn , sample_freq)
    state = get_state(track_res)
    decode(decoder, get_bits(track_res), get_num_bits(track_res))
end
```
Output: 
```t
Decoding subframe 1...
Decoding subframe 2...
Decoding subframe 3...
Decoding subframe 4...
Decoding subframe 5...
DECODING COMPLETED!
```
#### Satellite Position:
 Before using this function, you should check if computing is possible by using `can_get_sat_position`.
```julia
julia> can_get_sat_position(decoder)
true
```
You can compute the raw satellite position in two ways using `sat_position(dc::GNSSDecoderState, code_phase = 0)`. 

You can use only the decoder for simple and on 20ms accurate computing by using `sat_position(dc::GNSSDecoderState)`
```julia
julia> sat_position(decoder)
3-element StaticArrays.SArray{Tuple{3},Float64,1,3} with indices SOneTo(3):
 1.7583851256102066e7
 6.145733031320688e6
 1.8857912273122527e7
```
The precise position needs the codephase of the satellite signal. This module can´t get the codephase, but 
[Tracking](https://github.com/JuliaGNSS/Tracking.jl "JuliaGNSS:Tracking.jl") does.
```julia
julia> using Tracking, GNSSSignals, Unitful: s, Hz
julia> code_phase = Tracking.get_code_phase(track_res)
0.012411024547025855
julia> sat_postition = sat_position(decoder, code_phase)
3-element StaticArrays.SArray{Tuple{3},Float64,1,3} with indices SOneTo(3):
 1.758386392765011e7
 6.145759965494706e6
 1.885789227720346e7
```

#### Pseudorange:
You can compute the Pseudorange of the Satellite by using `pseudo_range(dc::GNSSDecoderState, code_phase::Float64)`. 
```julia
julia> using GNSSSignals
julia> pseudo_range(decoder, code_phase)
1.3825384018646529e9
```