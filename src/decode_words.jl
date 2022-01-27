struct ParityError <: Exception
    var
    text
end

"""
    Checks the Buffer for Parity Errors
    $(SIGNATURES)

    ´dc´: decoder, struct for all ephemeris data
    ´words_in_subframe´: buffer of 310 bits, which are again sliced in 10 Words á 30 Bits and 8 Bits prev preamble
     
    # Details
    # Checks every Word for Parity Errors
"""
function check_buffer(dc::GNSSDecoderState, words_in_subframe)
    
    parity_buffer = words_in_subframe
    
    dc.data_integrity = true
    for i in 1:length(parity_buffer)
        if i != 1
            dc.prev_29 = parity_buffer[i - 1][29]
            dc.prev_30 = parity_buffer[i - 1][30]
        end
        if dc.prev_30 
            parity_buffer[i][1:24] = .!parity_buffer[i][1:24]
        end
        if !parity_check(parity_buffer[i], dc.prev_29, dc.prev_30)
            @warn("Parity Error: Word , ", i, dc.PRN)
            dc.data_integrity = false
        end
    end

    return dc
end

"""
    Checks Data for Errors (IODC and matching IODEs)
    $(SIGNATURES)

    ´data´: Data

    # Details
    # Checks if IODE of subframe 1 matches the 8 LSB of th in Subframe 1 transmitted IODC
"""
function check_data_sub1_and_sub2(data::GPSData,
    subfr_1_data::Subframe_1_Data)

    has_error = !isnothing(data.IODE_Sub_2) && subfr_1_data.IODC[3:10] != data.IODE_Sub_2 ||
        !isnothing(data.IODE_Sub_3) && subfr_1_data.IODC[3:10] != data.IODE_Sub_3 # IODE and the 8 LSB of IODC must match
    if has_error
        @info "New data required, IODC and IODE of Subframe 2 and 3 do not match (CEI data set cutover)"
    end
    return !has_error
end

"""
    Checks Data for Errors (IODC and matching IODEs)
    $(SIGNATURES)

    ´data´: Data

    # Details
    # Checks if IODE of subframe 2 matches the 8 LSB of th in Subframe 1 transmitted IODC
"""
function check_data_sub2(data::GPSData,
    subfr_2_data::Subframe_2_Data)

    has_error = !isnothing(data.IODC) && subfr_2_data.IODE != data.IODC[3:10] # IODE and the 8 LSB of IODC must match
    if has_error
        @info "New data required, IODC and IODE of Subframe 2 do not match (CEI data set cutover)"
    end
    return !has_error
end

"""
    Checks Data for Errors (IODC and matching IODEs)
    $(SIGNATURES)

    ´data´: Data

    # Details
    # Checks if IODE of subframe 3 matches the 8 LSB of th in Subframe 1 transmitted IODC
"""
function check_data_sub3(data::GPSData,
    subfr_3_data::Subframe_3_Data)

    has_error = !isnothing(data.IODC) && subfr_3_data.IODE != data.IODC[3:10] # IODE and the 8 LSB of IODC must match
    if has_error
        @info "New data required, IODC and IODE of Subframe 3 do not match (CEI data set cutover)"
    end
    return !has_error
end



"""
    Decodes words of subframe 1-3, TLM and HOW of subframe 1-5
    $(SIGNATURES)

    ´dc´: decoder, struct for all Satellite Data
    ´words_in_subframe´: buffer of 310 bits, which are again sliced in 10 Words á 30 Bits and 8 Bits prev preamble


    # Details
    # Decodes complete subframe from the buffer, saves data for position computing in ´dc.data´. It returns the decoder, ´dc´. 
    # The number of saved Bits is resetted to 10 due to the Buffer length of 308. Those 10 Bits will be the previous 2 Bits + 8 Bits TOW.  
"""
function decode_words(dc::GNSSDecoderState, subframe_id, words_in_subframe, debug)

    dc = check_buffer(dc, words_in_subframe)
    if dc.data_integrity == true
        prev_subframes_decoded = copy(dc.subframes_decoded)
        if subframe_id == 1
            TLM_HOW, subframe_1_data = decode_subframe_1(words_in_subframe, debug)
            dc.data = GPSData(dc.data, subframe_1_data)
            dc.data = GPSData(dc.data, TLM_HOW)
            are_correct = check_data_sub1_and_sub2(dc.data, subframe_1_data)
            dc.subframes_decoded[1] = true
            dc.subframes_decoded[2] = dc.subframes_decoded[2] && are_correct
            dc.subframes_decoded[3] = dc.subframes_decoded[3] && are_correct
            if debug
                println("DECODED SUB1!")
            end
        elseif subframe_id == 2
            TLM_HOW, subframe_2_data = decode_subframe_2(words_in_subframe, debug)
            is_correct = check_data_sub2(dc.data, subframe_2_data)
            dc.subframes_decoded[2] = is_correct
            if is_correct
                dc.data = GPSData(dc.data, subframe_2_data)
                dc.data = GPSData(dc.data, TLM_HOW)
                if debug
                    println("DECODED SUB2!")
                end
            end
        elseif subframe_id == 3
            TLM_HOW, subframe_3_data = decode_subframe_3(words_in_subframe, debug)
            is_correct = check_data_sub3(dc.data, subframe_3_data)
            dc.subframes_decoded[3] = is_correct
            if is_correct
                dc.data = GPSData(dc.data, subframe_3_data)
                dc.data = GPSData(dc.data, TLM_HOW)
                if debug
                    println("DECODED SUB3!")
                end
            end
        elseif subframe_id == 4
            TLM_HOW = decode_subframe_4(words_in_subframe, debug)
            dc.subframes_decoded[4] = true
            dc.data = GPSData(dc.data, TLM_HOW)
            if debug
                println("DECODED SUB4!")
            end
        elseif subframe_id == 5
            TLM_HOW = decode_subframe_5(words_in_subframe, debug)
            dc.data = GPSData(dc.data, TLM_HOW)
            dc.subframes_decoded[5] = true
            if debug
                println("DECODED SUB5!")
            end
        end

        if prev_subframes_decoded != dc.subframes_decoded == [true, true, true, true, true] && debug
            println("PRN",dc.PRN," DATA COMPLETED OR RENEWED!")
        end
    end

    dc.num_bits_buffered = 10 # Due to Buffer length of 1502 and the need to save the last ten bits, only 1500 Bits have to been read in the next loop.
    return dc
end

