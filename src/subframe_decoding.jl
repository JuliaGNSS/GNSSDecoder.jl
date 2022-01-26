
GPS_PI = 3.1415926535898

#include("gpsl1_structs.jl")

"""
    Decodes words of the first two Words of a subframe
    $(SIGNATURES)

    ´words´: buffer of 60 bits, sliced in 2 Words á 30 Bits


    # Details
    # Decodes the first two Words of the Substring and returns a custom struct for these values
"""
function decode_TLM_HOW(words)
    integrity_status = words[1][23]
    tow = convert(Float64, bin2dec(words[2][1:17]))
    alert_flag = words[2][18]
    if alert_flag
        @warn "Signal URA may be worse than indicated in subframe 1 - Use satellite at own risk!"
    end
    anti_spoof_flag = words[2][19]
    return TLM_HOW_Data_Struct(integrity_status, tow, alert_flag, anti_spoof_flag)
end


"""
    Decodes words of the first subframe
    $(SIGNATURES)

    words: buffer of 300 bits, sliced in 10 Words á 30 Bits


    # Details
    # Decodes the first substring of the Substring and returns a custom struct for these values alongside the computed HOW and TLM Data
"""
function decode_subframe_1(words)
    if PRINT_ON
        println("Decoding subframe 1...")
    end

    TLM_HOW = decode_TLM_HOW(words[1:2])
    # * Decoding Word 3
    # Transmission Week
    trans_week = convert(Float64, bin2dec(words[3][1:10]))

    # Codes on L2 Channel
    codeonl2 = convert(Float64, bin2dec(words[3][11:12]))

    if codeonl2 == 3 | 0
        @warn "Code on L2 Channel invalid!"
    end

    # SV Accuracy, user range accuracy
    ura  = bin2dec(words[3][13:16])
    if ura <= 6
        ura = 2^(1 + (ura / 2))
    elseif 6 < ura <= 14
        ura = 2^(ura - 2)
    elseif ura == 15
        @warn "URA unsafe, no accuracy prediction available - use Satellite on own risk!"
        ura = 99999
    end

    # Satellite Health
    svhealth = bitArray2Str(words[3][17:22])
    if words[3][17]
        @warn "Bad LNAV Data, SV-Health critical"
    end

    # Issue of Data Clock
    IODC = bitArray2Str(append!(words[3][23:24], words[8][1:8])) # 2 MSB in Word 2, LSB 8 in Word 8


    # * Decoding Word 4
    # True: LNAV Datastream on PCode commanded OFF
    l2pcode = words[4][1]


    # * Decoding Word 7
    # group time differential
    T_GD = convert(Float64, bin2dec_twoscomp(words[7][17:24])) * 2^-31


    # *Decoding Word 8
    # IODC already in computed in word 3

    # Clock data reference
    t_oc = bin2dec(words[8][9:24]) << 4


    # * Decoding Word 9
    # clock correction parameter a_f2
    a_f2 = bin2dec_twoscomp(words[9][1:8]) * 2^-55

    # clock correction parameter  a_f1
    a_f1 = bin2dec_twoscomp(words[9][9:24]) * 2^-43


    # * Decoding Word 10
    # clock correction parameter a_f0
    a_f0 = bin2dec_twoscomp(words[10][1:22]) * 2^-31

    # * Finish Decoding

    return TLM_HOW, Subframe_1_Data(trans_week, codeonl2, ura, svhealth, IODC, l2pcode, T_GD, t_oc, a_f2, a_f1, a_f0)
end


"""
    Decodes words of the second subframe
    $(SIGNATURES)

    ´words´: buffer of 300 bits, sliced in 10 Words á 30 Bits


    Details:
    # Decodes the second substring of the Substring and returns a custom struct for these values alongside the computed HOW and TLM Data
"""
function decode_subframe_2(words)
    if PRINT_ON
        println("Decoding subframe 2...")
    end

    TLM_HOW = decode_TLM_HOW(words[1:2])
    # * Decoding Word 3
    # Issue of ephemeris data
    IODE = bitArray2Str(words[3][1:8])
    
    
    # Amplitude of Sine Harmonic Correction Term to Orbit Radius
    C_rs = bin2dec_twoscomp(words[3][9:24]) * 2^-5



    # * Decoding Word 4
    # Mean motion difference from computed value
    Δn = bin2dec_twoscomp(words[4][1:16]) * GPS_PI * 2^-43 
    

    # Mean anomaly at Reference Time (From word 4 and 5)
    M_0 = bin2dec_twoscomp(append!(words[4][17:24], words[5][1:24])) * GPS_PI * 2^-31 
    


    # * Decoding Word 5
    # Mean time anomaly computed in word 4


    # * Decoding Word 6
    # Amplitude of the Cosine Harmonic Correction Term to the Argument Latitude
    C_uc = bin2dec_twoscomp(words[6][1:16]) * 2^-29 
    
    # Eccentricity
    e = bin2dec(append!(words[6][17:24], words[7][1:24])) * 2^-33


    # * Decoding Word 7
    # Eccentricity already computed in word 6


    # * Decoding Word 8
    # Amplitude of the Sine Harmonic Correction Term to the Argument of Latitude
    C_us = bin2dec_twoscomp(words[8][1:16]) * 2^-29

    # Square Root of Semi-Major Axis
    sqrt_A = bin2dec(append!(words[8][17:24], words[9][1:24])) * 2^-19
    

    # * Decoding Word 9
    # square of A already computed in Word 8


    # * Decoding Word 10
    # Reference Time ephemeris
    t_oe = bin2dec(words[10][1:16]) << 4

    # Curve ftir interval flag - (0: 4 hours | 1: greater than 4 hours)
    fitinterval = words[10][17]

    # AODO Word
    aodo = bin2dec(words[10][18:22])

    # * Finish Decoding
    return TLM_HOW, Subframe_2_Data(IODE, C_rs, Δn, M_0, C_uc, e, C_us, sqrt_A, t_oe, fitinterval, aodo)
end


"""
    Decodes words of the third subframe
    $(SIGNATURES)

    ´words´: buffer of 300 bits, sliced in 10 Words á 30 Bits


    # Details
    # Decodes the third substring of the Substring and returns a custom struct for these values alongside the computed HOW and TLM Data
"""
function decode_subframe_3(words)
    if PRINT_ON
        println("Decoding subframe 3...")
    end

    TLM_HOW = decode_TLM_HOW(words[1:2])
    # * Decoding Word 3
    # Amplitude of the Cosine Harmonic Correction to Angle of Inclination
    C_ic = bin2dec_twoscomp(words[3][1:16]) * 2^-29

    # Longitude of Ascending Node of Orbit Plane at Weekly Epoch
    Ω_0 = bin2dec_twoscomp(append!(words[3][17:24], words[4][1:24])) * GPS_PI * 2^-31 
    # * Decoding Word 4
    # Omega 0 already in word 3 computed

    # * Decoding Word 5
    # Amplitude of the sine harmonic correction term to angle of Inclination
    C_is = bin2dec_twoscomp(words[5][1:16]) * 2^-29

    # inclination Angle at reference time
    i_0 = bin2dec_twoscomp(append!(words[5][17:24], words[6][1:24])) * GPS_PI * 2^-31
    
    # * Decoding Word 6
    # i_0 already in Word 5 computed

    # * Decoding Word 7
    # Amplitude of the cosine harmonic correction term to orbit Radius
    C_rc = bin2dec_twoscomp(words[7][1:16]) * 2^-5

    # Argument of Perigee
    ω = bin2dec_twoscomp(append!(words[7][17:24], words[8][1:24])) * GPS_PI * 2^-31

    # * Decoding Word 8
    # Argument of Perigee already computed in word 7

    # * Decoding Word 9
    # Rate of Right Ascension
    Ω_dot = bin2dec_twoscomp(words[9][1:24]) * GPS_PI * 2^-43

    # * Decoding Word 10
    # Issue of Ephemeris Data
    IODE = bitArray2Str(words[10][1:8])

    # Rate of Inclination Angle
    IDOT = bin2dec_twoscomp(words[10][9:22]) * GPS_PI * 2^-43

    # * Finish Decoding
    return TLM_HOW, Subframe_3_Data(C_ic, Ω_0, C_is, i_0, C_rc, ω, Ω_dot, IODE, IDOT)
end

"""
    Decodes words of the fourth subframe
    $(SIGNATURES)

    ´words´: buffer of 300 bits, sliced in 10 Words á 30 Bits


    # Details
    # Does mot do anything at the Moment, except calling the TLM_HOW Decoder and returns its return
"""
function decode_subframe_4(words)
    if PRINT_ON
        println("Decoding subframe 4...")
    end

    TLM_HOW = decode_TLM_HOW(words[1:2])
    return TLM_HOW
end

"""
    Decodes words of the fifth subframe
    $(SIGNATURES)

    ´words´: buffer of 300 bits, sliced in 10 Words á 30 Bits


    # Details
    # Does mot do anything at the Moment, except calling the TLM_HOW Decoder and returns its return
"""
function decode_subframe_5(words)
    if PRINT_ON
        println("Decoding subframe 5...")
    end

    TLM_HOW = decode_TLM_HOW(words[1:2])
    return TLM_HOW
end