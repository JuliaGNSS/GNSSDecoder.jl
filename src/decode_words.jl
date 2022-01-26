struct ParityError <: Exception
    var
    text
end

PRINT_ON = false
new_DATA_PRINT = true

"""
    Controls the Buffer for Parity Errors
    $(SIGNATURES)

    ´dc´: decoder, struct for all ephemeris data
    ´words_in_subframe´: buffer of 310 bits, which are again sliced in 10 Words á 30 Bits and 8 Bits prev preamble
     
    # Details
    # Controls every Word for Parity Errors
"""
function buffer_control(dc::GNSSDecoderState, words_in_subframe)
    

    parity_buffer = words_in_subframe
    
    dc.data_integrity = true
    for i in 1:length(parity_buffer)                    
        if i == 1
            if dc.prev_30 
                parity_buffer[i][1:24] = .!parity_buffer[i][1:24]
            end
            if parity_check(parity_buffer[i], dc.prev_29, dc.prev_30)
                
            else
                @warn("Parity Error: Word , ", i, dc.PRN)
                dc.data_integrity = false
            end
        else
            dc.prev_29 = parity_buffer[i - 1][29]
            dc.prev_30 = parity_buffer[i - 1][30]

            if dc.prev_30 
                parity_buffer[i][1:24] = .!parity_buffer[i][1:24]
            end
            if parity_check(parity_buffer[i], dc.prev_29, dc.prev_30)
            else
                @warn("Parity Error: Word , ", i, dc.PRN)
                dc.data_integrity = false
            end
        end
    end

    return dc
end



"""
    Controls Data for Errors (IODC and matching IODEs)
    $(SIGNATURES)

    ´TLM_HOW_Data´: Data of TLM and HOW Words of last read in subframe
    ´subframe_1_data´: Data of Subframe 1
    ´subframe_2_data´: Data of Subframe 2
    ´subframe_3_data´: Data of Subframe 3

    # Details
    #checks if IODE of subframe 2 and 3 matches the 8 LSB of th in Subframe 1 transmitted IODC
"""
function control_data(TLM_HOW_Data::TLM_HOW_Data_Struct,
    subfr_1_data::Subframe_1_Data,
    subfr_2_data::Subframe_2_Data,
    subfr_3_data::Subframe_3_Data)

    status = true
    if subfr_2_data.IODE != subfr_1_data.IODC[3:10] # IODE and the 8 LSB of IODC must match
        @info "new data required, IODC and IODE of Subframe 2 do not match (CEI data set cutover)"
        status = false
    end
    

    if subfr_3_data.IODE != subfr_1_data.IODC[3:10] # IODE and the 8 LSB of IODC must match
        @info "new data required, IODC and IODE of Subframe 3 do not match (CEI data set cutover)"
        status = false
    end


    return status
end

"""
    Controls Data for Errors (IODC and matching IODEs)
    $(SIGNATURES)

    ´next_data´: Data of next potential data

    # Details
    #checks if IODE of subframe 2 matches the 8 LSB of th in Subframe 1 transmitted IODC
"""
function control_data_sub2(next_data::GPSData,
    subfr_2_data::Subframe_2_Data)

    status = true
    if subfr_2_data.IODE != next_data.IODC[3:10] # IODE and the 8 LSB of IODC must match
        @info "new data required, IODC and IODE of Subframe 2 do not match (CEI data set cutover)"
        status = false
    end
    return status
end

"""
    Controls Data for Errors (IODC and matching IODEs)
    $(SIGNATURES)

    ´next_data´: Data of next potential data

    # Details
    #checks if IODE of subframe 3 matches the 8 LSB of th in Subframe 1 transmitted IODC
"""
function control_data_sub3(next_data::GPSData,
    subfr_3_data::Subframe_3_Data)

    status = true
    if subfr_3_data.IODE != next_data.IODC[3:10] # IODE and the 8 LSB of IODC must match
        @info "new data required, IODC and IODE of Subframe 3 do not match (CEI data set cutover)"
        status = false
    end
    return status
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
function decode_words(dc::GNSSDecoderState, words_in_subframe)

    control_state = false

    dc = buffer_control(dc, words_in_subframe) 
    if dc.data_integrity == true
        if dc.subframes_decoded_new == [0,0,0,0,0]
            TLM_HOW, subframe_1_data = decode_subframe_1(words_in_subframe)
            dc.data_next = GPSData(TLM_HOW, subframe_1_data)
            control_state = true
            dc.subframes_decoded_new[1] = control_state
            dc.data = GPSData(dc.data, TLM_HOW)
            if PRINT_ON
                println("DECODED SUB1!")
            end
        elseif dc.subframes_decoded_new == [1,0,0,0,0]
            TLM_HOW, subframe_2_data = decode_subframe_2(words_in_subframe)
            control_state = control_data_sub2(dc.data_next, subframe_2_data)
            dc.subframes_decoded_new[2] = control_state
            if control_state
                dc.data_next = GPSData(dc.data_next, TLM_HOW, subframe_2_data)
                dc.data = GPSData(dc.data, TLM_HOW)
                if PRINT_ON
                    println("DECODED SUB2!")
                end
            end
        elseif dc.subframes_decoded_new == [1,1,0,0,0]
            TLM_HOW, subframe_3_data = decode_subframe_3(words_in_subframe)
            control_state = control_data_sub3(dc.data_next, subframe_3_data)
            dc.subframes_decoded_new[3] = control_state
            if control_state
                dc.data_next = GPSData(dc.data_next, TLM_HOW, subframe_3_data)
                dc.data = GPSData(dc.data, TLM_HOW)
                if PRINT_ON
                    println("DECODED SUB3!")
                end
            end
        elseif dc.subframes_decoded_new == [1,1,1,0,0]
            TLM_HOW = decode_subframe_4(words_in_subframe)
            control_state = true
            dc.subframes_decoded_new[4] = control_state
            dc.data = GPSData(dc.data, TLM_HOW)
            if PRINT_ON
                println("DECODED SUB4!")
            end
        elseif dc.subframes_decoded_new == [1,1,1,1,0]
            TLM_HOW = decode_subframe_5(words_in_subframe)
            dc.data_next = GPSData(dc.data_next, TLM_HOW)
            control_state = true
            dc.subframes_decoded_new[5] = control_state
            if PRINT_ON
                println("DECODING OF FRAME FINISHED!")
            end
        end

        if dc.subframes_decoded_new == [1,1,1,1,1]
            dc.subframes_decoded = dc.subframes_decoded_new
            dc.subframes_decoded_new = [0,0,0,0,0]
            dc.data = dc.data_next
            if PRINT_ON || new_DATA_PRINT
                println("PRN",dc.PRN," NEW DECODED DATA!")
            end
        elseif !control_state
            dc.subframes_decoded_new = [0,0,0,0,0]
        end
    end

    dc.num_bits_buffered = 10 # Due to Buffer length of 1502 and the need to save the last ten bits, only 1500 Bits have to been read in the next loop.
    return dc
end

