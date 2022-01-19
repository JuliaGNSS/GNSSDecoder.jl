
PREAMBLE = BitArray([1,1,0,1,0,0,0,1]) # Reversed preamble through reversed buffer

FRAME_POSITIONS = [8,308] #NEEDED for subframe wise decoding


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
    
    Checks whether there is an preamble at each FRAME_POSITION of the buffer or not. 
"""
function find_preamble(buffer)
    for pos in FRAME_POSITIONS 
        # Verify preamble
        data = access_buffer(buffer, pos, 0, 8)
        if data != PREAMBLE && data != .!PREAMBLE
            return false
        end
        # -> Preambles at expected Position
    end
    #println("DECODING...")
    return true
end

"""
    Reads ID of Subframe in sorted buffer
    $(SIGNATURES)

    #Details
    
    Reads the ID of the, in buffer stored, subframe. Buffer must be in the order, that at the begining and the end are preambles.
"""
function read_id(buffer)
    pos = FRAME_POSITIONS[2]
    inverting_bit = access_buffer(buffer, pos, 30, 1)[1]
    id = access_buffer(buffer, pos, 49, 3)
    id = reverse(id) #Databit are reversed in Buffer
    if inverting_bit #If last Bit of previous Word is set the next Word gets inverted -> ID gets inverted
        id = .!id
    end

    id = bin2dec(id)
end


