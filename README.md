# GNSSDecoder.jl

Decodes GPSL1 satellite signals and computes the Satellite position in WGS84.
## Features
* Decoding of Satellite Signals (Ephemeris)
## Usage

#### Install: 
```julia
julia> ]
(@v1.4) pkg> add GNSSDecoder
```


### Prerequiries:

To set up the Workspace, the module
* [Tracking.jl](https://github.com/JuliaGNSS/Tracking.jl "JuliaGNSS:Tracking.jl")

is required. The definition of sample frequency and the size of the Array of Signals finalize the Preparation. 
```julia
using GNSSDecoder, Tracking
using Tracking: Hz, s
```

### Initialization
Before Decoding, the signals have to be acquired for each PRN. The decoder is initialized in this state.
```julia
decoder = GNSSDecoderState(PRN = 1) #Initialization of Decoder
```



### Decoding
Before Decoding, the signals have to be acquired for each used PRN to initialize primary `TrackingStates`. This function is defined in the Module [Tracking.jl](https://github.com/JuliaGNSS/Tracking.jl "JuliaGNSS:Tracking.jl"). 
```julia
# ´carrier_doppler´: Doppler of carrier signal
# ´code_phase´: code phase
tracking_states = TrackingState(GPSL1, carrier_doppler, code_phase)
```


The filtered Bits for decoding are supplied by the Tracking module.
```julia
# `iterations`: number of iterations
for i in 1:iterations
    # Read in new signal
    
    #sample_freq: sample frequency in Hz
    track_res = track(signal, tracking_states, decoder.PRN , sample_freq)
    tracking_states = get_state(track_res)
    decode(decoder, get_bits(track_res), get_num_bits(track_res))
end
```


The following Output for each PRN will show: 
```t
DECODING...
Decoding subframe 1...
Decoding subframe 2...
Decoding subframe 3...
Decoding subframe 4...
Decoding subframe 5...
DECODING COMPLETED!
```
