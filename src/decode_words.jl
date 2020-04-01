struct ParityError <: Exception
    var
    text
end

"""
    Decodes the buffer, decodes words for subframe 1-3
    $(SIGNATURES)

    ´dc´: decoder, struct for all ephemeris data
    ´words´: buffer of 308 bits, containing one complete subframe and the last preamble for parity checking
     
    # Details
    # Decodes complete subframe from the buffer, saves data for position computing in ´dc.data´ 
"""
function decode_words(dc::GNSSDecoderState, words)

    # _invert_preamble = BitArray([false, false, true, false, true, true, true, false])

    # if words[1][1:8] == _invert_preamble
    #     words[1][1:24] = .!words[1][1:24]
    # end

    # * Invert word if D30 of previous word is set
    for i in 1:length(words)
        if i == 1
            if dc.prev_30 
                words[i][1:24] = .!words[i][1:24]
            end
            if parity_check(words[i], dc.prev_29, dc.prev_30)
                dc.data_integrity = true
            else
                @warn("Parity Error: Word , ", i)
                dc.data_integrity = false
            end
        else
            dc.prev_29 = words[i - 1][29]
            dc.prev_30 = words[i - 1][30]

            if dc.prev_30 
                words[i][1:24] = .!words[i][1:24]
            end
            if parity_check(words[i], dc.prev_29, dc.prev_30)
                dc.data_integrity = true
            else
                @warn("Parity Error: Word , ", i)
                dc.data_integrity = false
            end
        end
    end
    


    # * Decode HOW - Word
    tow = words[2][1:17]
    dc.data.tow = convert(Float64, bin2dec(tow))
    tow_computed = dc.data.tow * 6
    # if tow_computed > 302400
    #     tow_computed = tow_computed - 604800
    # elseif tow_computed < -302400
    #     tow_computed = tow_computed + 604800
    # end
    dc.data.tow_computed = tow_computed
    dc.data.alert = words[2][18]
    if dc.data.alert
        @warn "Signal URA may be worse than indicated in subframe 1 - Use satellite at own risk!"
    end
    dc.data.antiSpoof = words[2][19]
    




    subfr_bits = words[2][20:22] 

    if subfr_bits == [0,0,1]
        dc.subframe_count = 1
        println("Decoding subframe ", dc.subframe_count, "...")

        # * Decoding Word 3
        # Transmission Week
        trans_week = bin2dec(words[3][1:10])
        dc.data.transmissionWeekNumber = convert(Float64, trans_week)

        # Codes on L2 Channel
        codeonl2 = bin2dec(words[3][11:12])
        dc.data.l2codes = convert(Float64, codeonl2)

        if codeonl2 == 3 | 0
            @warn "Code on L2 Channel invalid!"
        end

        # SV Accuracy, user range accuracy
        ura  = bin2dec(words[3][13:16])
        if ura <= 6
            dc.data.ura = 2^(1 + (ura / 2))
        elseif 6 < ura <= 14
            dc.data.ura = 2^(ura - 2)
        elseif ura == 15
            @warn "URA unsafe, no accuracy prediction available - use Satellite on own risk!"
            dc.data.ura = 99999
        end

        # Satellite Health
        svhealth = bitArray2Str(words[3][17:22])
        dc.data.svHealth = svhealth
        if words[3][17]
            @warn "Bad LNAV Data, SV-Health critical"
        end

        # Issue of Data Clock
        IODC = append!(words[3][23:24], words[8][1:8]) # 2 MSB in Word 2, LSB 8 in Word 8
        dc.data.IODC = bitArray2Str(IODC)


        # * Decoding Word 4
        # True: LNAV Datastream on PCode commanded OFF
        dc.data.l2pcode = words[4][1]


        # * Decoding Word 7
        # group time differential
        grp_delay_differential = bin2dec_twoscomp(words[7][17:24])
        dc.data.groupDelayDifferential = convert(Float64, grp_delay_differential) * 2^-31


        # *Decoding Word 8
        # IODC already in computed in word 3

        # Clock data reference
        t_oc = bin2dec(words[8][9:24]) << 4
        dc.data.toc = t_oc


        # * Decoding Word 9
        # clock correction parameter a_f2
        a_f2 = bin2dec_twoscomp(words[9][1:8]) * 2^-55
        dc.data.af2 = a_f2

        # clock correction parameter  a_f1
        a_f1 = bin2dec_twoscomp(words[9][9:24]) * 2^-43
        dc.data.af1 = a_f1


        # * Decoding Word 10
        # clock correction parameter a_f0
        a_f0 = bin2dec_twoscomp(words[10][1:22]) * 2^-31
        dc.data.af0 = a_f0

        # * Finish Decoding
        dc.subframes_decoded[1] = true

        


    elseif subfr_bits == [0,1,0]
        dc.subframe_count = 2
        println("Decoding subframe ", dc.subframe_count, "...")
        



        # * Decoding Word 3
        # Issue of ephemeris data
        IODE = words[3][1:8]
        dc.data.IODE = bitArray2Str(IODE)
        
        if dc.data.IODE != nothing && dc.data.IODC != nothing  
            if dc.data.IODE != dc.data.IODC[3:10]
                @info "new data required, IODC and IODE do not match"
                dc.new_data_needed = true
            end
        end

        # Amplitude of Sine Harmonic Correction Term to Orbit Radius
        c_rs = bin2dec_twoscomp(words[3][9:24]) * 2^-5
        dc.data.C_rs = c_rs



        # * Decoding Word 4
        # Mean motion difference from computed val
        Δn = bin2dec_twoscomp(words[4][1:16]) * 2^-43 * π
        dc.data.Δn = Δn

        # Mean anomaly at Reference Time (From word 4 and 5)
        M = bin2dec_twoscomp(append!(words[4][17:24], words[5][1:24])) * 2^-31 * π
        dc.data.M0 = M


        # * Decoding Word 5
        # Mean time anomaly computed in word 4


        # * Decoding Word 6
        # Amplitude of the Cosine Harmonic Correction Term to the Argument Latitude
        dc.data.C_uc = bin2dec_twoscomp(words[6][1:16]) * 2^-29 
        
        # Eccentricity
        dc.data.e = bin2dec(append!(words[6][17:24], words[7][1:24])) * 2^-33
    

        # * Decoding Word 7
        # Eccentricity already computed in word 6


        # * Decoding Word 8
        # Amplitude of the Sine Harmonic Correction Term to the Argument of Latitude
        dc.data.C_us = bin2dec_twoscomp(words[8][1:16]) * 2^-29

        # Square Root of Semi-Major Axis
        dc.data.sqrt_A = bin2dec(append!(words[8][17:24], words[9][1:24])) * 2^-19


        # * Decoding Word 9
        # square of A already computed in Word 8


        # * Decoding Word 10
        # Reference Time ephemeris
        dc.data.t_oe = bin2dec(words[10][1:16]) << 4

        # Curve ftir interval flag - (0: 4 hours| 1: greater than 4 hours)
        dc.data.fitinterval = words[10][17]

        # AODO Word
        dc.data.aodo = bin2dec(words[10][18:22])

        # * Finish Decoding
        dc.subframes_decoded[2] = true

    elseif subfr_bits == [0,1,1]
        dc.subframe_count = 3
        println("Decoding subframe ", dc.subframe_count, "...")

        # * Decoding Word 3
        # Amplitude of the Cosine Harmonic Correction to Angle of Inclination
        dc.data.C_ic = bin2dec_twoscomp(words[3][1:16]) * 2^-29

        # Longitude of Ascending Node of Orbit Plane at Weekly Epoch
        dc.data.𝛀_0 = bin2dec_twoscomp(append!(words[3][17:24], words[4][1:24])) * 2^-31 * π

        # * Decoding Word 4
        # Omega 0 already in word 3 computed

        # * Decoding Word 5
        # Amplitude of the sine harmonic correction term to angle of Inclination
        dc.data.C_is = bin2dec_twoscomp(words[5][1:16]) * 2^-29

        # inclination Angle at reference time
        dc.data.i_0 = bin2dec_twoscomp(append!(words[5][17:24], words[6][1:24])) * 2^-31 * π
        
        # * Decoding Word 6
        # i_0 already in Word 5 computed

        # * Decoding Word 7
        # Amplitude of the cosine harmonic correction term to orbit Radius
        dc.data.C_rc = bin2dec_twoscomp(words[7][1:16]) * 2^-5

        # Argument of Perigee
        dc.data.ω = bin2dec_twoscomp(append!(words[7][17:24], words[8][1:24])) * 2^-31 * π

        # * Decoding Word 8
        # Argument of Perigee already computed in word 7

        # * Decoding Word 9
        # Rate of Right Ascension
        dc.data.𝛀_dot = bin2dec_twoscomp(words[9][1:24]) * 2^-43 * π

        # * Decoding Word 10
        # Issue of Ephemeris Data
        IODE = words[10][1:8]
        dc.data.IODE = bitArray2Str(IODE)
        
        if dc.data.IODE != nothing && dc.data.IODC != nothing  
            if dc.data.IODE != dc.data.IODC[3:10]
                @info "new data required, IODC and IODE do not match"
                dc.new_data_needed = true
            end
        end

        # Rate of Inclination Angle
        dc.data.IDOT = bin2dec_twoscomp(words[10][9:22]) * 2^-43 * π

        # * Finish Decoding
        dc.subframes_decoded[3] = true



    elseif subfr_bits == [1,0,0]
        dc.subframe_count = 4
        println("Decoding subframe ", dc.subframe_count, "...")


        # * Finish Decoding
        dc.subframes_decoded[4] = true
    elseif subfr_bits == [1,0,1]
        dc.subframe_count = 5
        println("Decoding subframe ", dc.subframe_count, "...")


        # * Finish Decoding
        dc.subframes_decoded[5] = true
        if dc.subframes_decoded == [true, true, true, true, true]
            println("DECODING COMPLETED!")
        end
    end

    

    dc.num_bits_buffered = 8
    return dc

end
