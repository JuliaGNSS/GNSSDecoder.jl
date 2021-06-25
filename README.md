# GNSSDecoder.jl

Decodes GPSL1 satellite signals.

## Usage

#### Install: 
```julia
julia> ]
pkg> add GNSSDecoder
```

### Initialization
The decoder must be initialized beforehand.
```julia
decoder = GNSSDecoderState(PRN = 1) #Initialization of Decoder
```

### Decoding
Pass bits to decoder as an integer value and let the decoder decode the message.
```julia
for i in 1:iterations
    # Track signal for example with Tracking.jl
    track_res = track(signal, track_state, decoder.PRN , sampling_freq)
    track_state = get_state(track_res)
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
