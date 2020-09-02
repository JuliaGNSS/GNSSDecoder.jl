
PREAMBLE = BitArray([1,1,0,1,0,0,0,1]) # Reversed preamble through reversed buffer


FRAME_POSITIONS = [300, 600, 900, 1200, 1500]

"""
    Accesses Buffer at a position with an offset and length of searched word
    $(SIGNATURES)

    #Details
    Helped troubleshooting though you don´t get confused about indizes
"""
function access_buffer(buffer, pos, offset, len)
    return buffer[pos - offset - len + 1:pos - offset]
end
"""
    Searches for preamble of Subframe
    $(SIGNATURES)

    #Details
    
    Checks whether there is an preamble at each FRAME_POSITION of the buffer or not. Checks for right order of subframes too by checking the ID. 
    To prevent combining data of different datasets the Order must always be subframe 1,2,3,4,5. Because the Buffer is reversed, it is checked for 
    the sequence 5,4,3,2,1. 
"""
function find_preamble(buffer)
    expected_seq = 5
    for pos in FRAME_POSITIONS
        # Verify preamble
        data = access_buffer(buffer, pos, 0, 8)
        if data != PREAMBLE && data != .!PREAMBLE
            return false
        end
        # -> Preambles at expected Position

        # Verify order of sequence numbers
        inverting_bit = access_buffer(buffer, pos, 30, 1)[1]
        id = access_buffer(buffer, pos, 49, 3)
        id = reverse(id) #Databit are reversed in Buffer
        if inverting_bit #If last Bit of previous Word is set the next Word gets inverted -> ID gets inverted
            id = .!id
        end

        id = bin2dec(id) 
        if id != expected_seq
            return false
        end
        expected_seq -= 1
    end
    println("DECODING...")
    return true
end

