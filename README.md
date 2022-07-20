# GNSSDecoder.jl

Decodes various GNSS satellite signals.
Currently implemented:
 * GPS L1

## Usage

#### Install: 
```julia
julia> ]
pkg> add git@github.com:JuliaGNSS/GNSSDecoder.jl.git
```

### Initialization
The decoder must be initialized beforehand.
```julia
decoder = GPSL1DecoderState(1) #Initialization of decoder with PRN = 1
```

### Decoding
Pass bits to decoder as an unsigned integer value and let the decoder decode the message.
```julia
for i in 1:iterations
    # Track signal for example with Tracking.jl
    track_res = track(signal, track_state, decoder.PRN , sampling_freq)
    track_state = get_state(track_res)
    decoder = decode(decoder, get_bits(track_res), get_num_bits(track_res))
end
```

The data can be retrieved by
```julia
decoder.data
```

Note that GNSSDecoder decodes each time a complete subframe has been retrieved.
`decoder.raw_data` holds the raw data. `decoder.data` hold data that has been checked for consistency.
`decoder.num_bits_after_valid_subframe` counts the number of bits after a valid subframe has been retrieved.
