module GNSSDecoder
    
    using DocStringExtensions, GNSSSignals, Parameters, FixedPointNumbers, StaticArrays, LinearAlgebra
    using Unitful: Hz
    
    const SUBFRAME_LENGTH = 300 # Size of Subframes
    const BUFFER_LENGTH = SUBFRAME_LENGTH + 2 + 8 # 2 Bits previous subfram + 300 Bits current subframe + 8 Bits following TOW
    const WORD_LENGTH = 30

    export decode,
        GNSSDecoderState,
        GPSData,
        GPSL1Constants
    
    include("gpsl1_structs.jl")
    
    """
        stores element in a bitbuffer

        $(SIGNATURES)
        
        # Buffer must be of Type BitArray{1}
    """
    function push_buffer(buffer::BitArray{1}, val::Bool)
        buffer = buffer >> 1
        buffer.chunks[1] = buffer.chunks[1] + val
        return buffer
    end


    """
        Decodes the bits inserted using a bit buffer
        $(SIGNATURES)

        # Details
        checks, if the buffer has the expected size of 310 bits (length of subframe + preamble + 2 prevbits),
        then tries to find a preamble.
        
    """
    function decode(dc::GNSSDecoderState, data_bits::UInt64, num_bits, debug = false)     
        
        for i = num_bits - 1:-1:0
            
            current_bit = data_bits & (1 << i)
            
            dc.buffer = push_buffer(dc.buffer, current_bit > 0)
            
            # Count if Buffer is new to full stored buffer
            if dc.num_bits_buffered != BUFFER_LENGTH 
                dc.num_bits_buffered += 1
            end

            if !isnothing(dc.data.TOW)
                dc.num_bits_after_valid_subframe += 1
            end

            # Every time the preamble is found at the correct positions
            # decode the current subframe
            if dc.num_bits_buffered == BUFFER_LENGTH && find_preamble(dc.buffer)
                rev_buf = reverse(dc.buffer)
                words_in_subframe  = map(word_idx -> rev_buf[word_idx * WORD_LENGTH + 3:(word_idx + 1) * WORD_LENGTH + 2], 0:9)
                dc.prev_29 = rev_buf[1]
                dc.prev_30 = rev_buf[2]
                decode_words(dc, read_id(dc.buffer), words_in_subframe, debug)
                if dc.data.integrity_status_flag == false
                    dc.num_bits_after_valid_subframe = 8
                end
            end
        end # end of for-loop 
    end # end of decode()

    
    include("subframe_decoding.jl") 
    include("find_preamble.jl") 
    include("decode_words.jl") 
    include("parity_check.jl")
    include("bin_2_dec.jl")
end# end of module