module GNSSDecoder
    
    using DocStringExtensions, GNSSSignals, Parameters, FixedPointNumbers, StaticArrays, LinearAlgebra
    using Unitful: Hz

    BUFFER_LENGTH = 1502 #Size of Buffer - All 5 Subframes + the last 2 Bits are stored
    SUBFRAME_LENGTH = 300 #Size of Subframes

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
        checks, if the buffer has the expected size of 1508 bits (length of subframe + preamble),
        then tries to find a preamble.
        
    """
    function decode(dc::GNSSDecoderState, data_bits::UInt64, num_bits)
    
        
        a = num_bits - dc.nb_prev
        dc.nb_prev = num_bits
        
        if a !== 0
            num_bits = a
        end
        
        
        for i = num_bits - 1:-1:0
            
            current_bit = data_bits & (1 << i)
            
            dc.buffer = push_buffer(dc.buffer, current_bit > 0)
            
            # * Count if Buffer is new to full stored buffer
            if dc.num_bits_buffered != BUFFER_LENGTH
                dc.num_bits_buffered += 1
            end

    
            # * Find first preamble
            if (!dc.preamble_found) && (dc.num_bits_buffered == BUFFER_LENGTH)
                dc.preamble_found = find_preamble(dc.buffer)
            end
        
            # * Begin decoding after preamble is confirmed
            if dc.preamble_found
                rev_buf = reverse(dc.buffer)
                words_in_subframes  = map( ##Splits Buffer in Subframes and Words
                        subfr_it->map(
                        word_it->rev_buf[(word_it * 30) + 3 + (SUBFRAME_LENGTH * subfr_it):((word_it + 1) * 30 ) + 2 + (SUBFRAME_LENGTH * subfr_it)],
                        0:9),
                        0:4) # gets the words of the subframe (word length = 30)
                 dc.prev_29 = rev_buf[1]
                 dc.prev_30 = rev_buf[2]
                 decode_words(dc, words_in_subframes)
                dc.preamble_found = false
            end
        end # end of for-loop 
    end # end of decode()

    
    include("subframe_decoding.jl")
    include("find_preamble.jl")
    include("decode_words.jl")
    include("parity_check.jl")
    include("bin_2_dec.jl")
end# end of module