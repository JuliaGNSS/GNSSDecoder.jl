struct ParityError <: Exception
    var
    text
end


"""
    Controls the Buffer for Parity Errors
    $(SIGNATURES)

    ´dc´: decoder, struct for all ephemeris data
    ´words_in_subframes´: buffer of 1500 bits, sliced in 5 subframes á 300 bits, which are again sliced in 10 Words á 30 Bits
     
    # Details
    # Controls every Word for Parity Errors
"""
function buffer_control(dc::GNSSDecoderState, words_in_subframes)
    

    parity_buffer = vcat(words_in_subframes[1], 
                        words_in_subframes[2],
                        words_in_subframes[3],
                        words_in_subframes[4],
                        words_in_subframes[5])
    
    dc.data_integrity = true
    for i in 1:length(parity_buffer)                    
        if i == 1
            if dc.prev_30 
                parity_buffer[i][1:24] = .!parity_buffer[i][1:24]
            end
            if parity_check(parity_buffer[i], dc.prev_29, dc.prev_30)
                
            else
                @warn("Parity Error: Word , ", i)
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
                @warn("Parity Error: Word , ", i)
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
    Decodes words of subframe 1-3, TLM and HOW of subframe 1-5
    $(SIGNATURES)

    ´dc´: decoder, struct for all Satellite Data
    ´words_in_subframes´: buffer of 1500 bits, sliced in 5 subframes á 300 bits, which are again sliced in 10 Words á 30 Bits


    # Details
    # Decodes complete subframe from the buffer, saves data for position computing in ´dc.data´. It returns the decoder, ´dc´. 
    #The number of saved Bits is resetted to 2 due to the Buffer length of 1502. Those 2 Bits will be the previous 2 Bits to the next 5 Subframes.  
"""
function decode_words(dc::GNSSDecoderState, words_in_subframes)

    dc = buffer_control(dc, words_in_subframes)
    if dc.data_integrity == true
        TLM_HOW, subframe_1_data = decode_subframe_1(words_in_subframes[1])
        dc.subframes_decoded[1] = true
        TLM_HOW, subframe_2_data = decode_subframe_2(words_in_subframes[2])
        dc.subframes_decoded[2] = true
        TLM_HOW, subframe_3_data = decode_subframe_3(words_in_subframes[3])
        dc.subframes_decoded[1] = true
        TLM_HOW = decode_subframe_4(words_in_subframes[4])
        dc.subframes_decoded[4] = true
        TLM_HOW = decode_subframe_5(words_in_subframes[5])
        dc.subframes_decoded[5] = true


        if control_data(TLM_HOW, subframe_1_data, subframe_2_data, subframe_3_data)  # Checks for Data Errors
            data = create_data(TLM_HOW, subframe_1_data, subframe_2_data, subframe_3_data) # Creates GPSData from single, subframe specific data and the TLM Word
            dc.data = data
            println(" DECODING COMPLETED!")
        end
    end

    dc.num_bits_buffered = 2 # Due to Buffer length of 1502 and the need to save the last two bits, only 1500 Bits have to been read in the next loop.
    return dc
end

