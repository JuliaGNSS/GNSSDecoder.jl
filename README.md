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
