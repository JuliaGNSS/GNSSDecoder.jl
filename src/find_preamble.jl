
_preamble = BitArray([true, true, false, true, false, false, false, true])
_inverted_preamble = BitArray([false, false, true, false, true, true, true, false])

"""
    Searches for preamble of Subframe
    $(SIGNATURES)

    #Details
    
    Checks whether there is an preamble at beginning and end of the buffer or not. 
"""
function find_preamble(buffer)

    buf_begin = buffer[1:8]
    buf_end = buffer[buffer.len - 9:buffer.len - 2]

    if buf_end == _preamble || buf_end == _inverted_preamble
        if buf_begin == _preamble || buf_begin == _inverted_preamble
            return true

        else
            return false
        end
    else
        return false
    end
end

    
