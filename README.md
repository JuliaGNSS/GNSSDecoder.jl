[![pipeline status](https://git.rwth-aachen.de/nav/GNSSDecoder.jl/badges/master/pipeline.svg)](https://git.rwth-aachen.de/nav/GNSSDecoder.jl/commits/master)
[![coverage report](https://git.rwth-aachen.de/nav/GNSSDecoder.jl/badges/master/coverage.svg)](https://git.rwth-aachen.de/nav/GNSSDecoder.jl/commits/master)

# Decode GNSS signals.

# Usage

decode = init_decode()

i = 1;

while true

       read!(stream,signal)
       global track, track_res = track(signal)
       decode(track_res.data_bits, track_res.num_bits)

       if i == 1500 + decode.found_preambles.preamble_pos
          break
       end

       if track_res.num_bits != 0
          global i +=1
       end
end

This loop will stop the moment a complete navigation message (i.e. 1500 bits) has been received (i.e. 1500 bits after the first preamble was reveceived)

decode.found_preambles.preamble_pos stands for the number of bits that were received before the first complete subframe
