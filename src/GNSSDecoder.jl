module GNSSDecoder
    
    using DocStringExtensions, GNSSSignals, Parameters, FixedPointNumbers, StaticArrays, LinearAlgebra
    using Unitful: Hz
    
    BUFFER_LENGTH = 310 ##NEEDS Modification e.g 308 Bit if there are two präamble, then the word wise decoding starts.
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
        checks, if the buffer has the expected size of 310 bits (length of subframe + preamble + 2 prevbits),
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

            current_id = 0
            # * Find first preamble
            if (!dc.preamble_found) && (dc.num_bits_buffered == BUFFER_LENGTH)
                dc.preamble_found = find_preamble(dc.buffer)
            end

            #Check Id of current subframe
            if dc.preamble_found
                current_id = read_id(dc.buffer)
                if current_id != sum(dc.subframes_decoded_new)+1 #checks if the current_id is correct
                    dc.preamble_found = false
                end
            end

            # * Begin decoding after preamble is confirmed
            if dc.preamble_found
                rev_buf = reverse(dc.buffer)
                words_in_subframe  = map(wrd_it -> rev_buf[ wrd_it*30+3 : (wrd_it+1)*30+2 ],0:9) # gets the words of the subframe (word length = 30)
                dc.prev_29 = rev_buf[1]
                dc.prev_30 = rev_buf[2]
                decode_words(dc, words_in_subframe)
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