function decodeword(word,data,parameters,intermediate)

    #If the last bit of the previous word is set, the payload bits in the next word are inversed
    if parameters.prev_30
    word = [word[1:6]; map(!, word[7:30])]
    end

    parameters.data_integrity=paritycheck(word,parameters.prev_29,parameters.prev_30)

    #Actualize prev_29 and prev_30
    parameters.prev_29 = word[2]
    parameters.prev_30 = word[1]

    if parameters.word_count==2
    #Bits 20,21 and 22 of the HOW provide the ID of the subframe
    subframe_bits = reverse(word[9:11])
    subframe_bits==[0;0;1] ? parameters.subframe_count=1 : (subframe_bits==[0;1;0] ? parameters.subframe_count=2 : (subframe_bits==[0;1;1] ? parameters.subframe_count=3 : (subframe_bits==[1;0;0] ? parameters.subframe_count=4 : parameters.subframe_count=5)))
    parameters.subframe_count == 1 ? parameters.first_subframe = true : (parameters.subframe_count == 2 ? parameters.second_subframe = true : (parameters.subframe_count == 3 ? parameters.third_subframe = true : (parameters.subframe_count == 4 ? parameters.fourth_subframe = true : parameters.fifth_subframe = true)))
        if  parameters.first_subframe && parameters.second_subframe && parameters.third_subframe && parameters.fourth_subframe && parameters.fifth_subframe
            parameters.decoding_completed = true
        end
    @show parameters.subframe_count
    #decode tow-count only if it has not been read yet
        if data.tow == nothing
            tow = reverse(word[14:30])
            data.tow = convert(Float64,bin2dec(tow))
            data.alert = word[13]
            data.antiSpoof = word[14]
        end
    end

    if parameters.word_count==3

        if parameters.subframe_count==1
            IODC_1 = reverse(word[7:8])
            data.IODC=bitArray2Str(IODC_1)

            transmissionWeekNumber = reverse(word[21:30])
            data.transmissionWeekNumber = convert(Float64,bin2dec(transmissionWeekNumber))
            #codes on L2 channel
            l2codes = reverse(word[19:20])
            data.l2codes = convert(Float64,bin2dec(l2codes))
            #SV Accuracy / User Range Accuracy
            uraIndex = reverse(word[15:18])
            uraIndex = convert(Float64,bin2dec(uraIndex))
            if uraIndex <6
                data.ura = 2^(1+uraIndex/2)
            else
                data.ura = 99999
            end

            svHealth = reverse(word[9:14])
            data.svHealth = bitArray2Str(svHealth)
        end
        #Ephemeris data are only in subframes 2 and 3
        if parameters.subframe_count==2
            IODE = reverse(word[23:30])
            data.IODE=bitArray2Str(IODE)

            #IODE must be equal to the 8 LSB of the IODC. If they don't match new data must be collected
            if data.IODE !== nothing && data.IODC !== nothing
            parameters.new_data_requ = !(data.IODE == data.IODC[3:10])
            end
            #Crs sine correction of orbit radius
            _Crs= reverse(word[7:22])
            data.C_rs= convert(Float64,bin2dec_twoscomp(_Crs))*2^-5
        end

        if parameters.subframe_count==3
            #Cic cosine correction to the angle of inclination
            _Cic= reverse(word[15:30])
            data.C_ic= convert(Float64,bin2dec_twoscomp(_Cic))*2^-29

            #Omega0 Longitude of Ascenting Node of Orbit Plane
            intermediate._𝛀_0= reverse(word[7:14])
        end


    end

    if parameters.word_count==4

        #data flag for L2 P-Code
        if parameters.subframe_count==1
        data.l2pcode = word[30]
        end

        if parameters.subframe_count==2
            _deltan= reverse(word[15:30])
            data.Δn= convert(Float64,bin2dec_twoscomp(_deltan))*2^-43*π
            intermediate._M0 = reverse(word[7:14])
        end

        if parameters.subframe_count==3
            #Omega0 Longitude of Ascenting Node of Orbit Plane
            intermediate._𝛀_0=vcat(intermediate._𝛀_0,reverse(word[7:30]))
            data.𝛀_0= convert(Float64,bin2dec_twoscomp(intermediate._𝛀_0))*2^-31*π

        end

    end

    if parameters.word_count==5

        if parameters.subframe_count==2
            intermediate._M0 = vcat(intermediate._M0,reverse(word[7:30]))
            data.M0= convert(Float64,bin2dec_twoscomp(intermediate._M0))*2^-31*π
        end

        if parameters.subframe_count==3
            #Cis sine correction to the angle of inclination
            _C_is= reverse(word[15:30])
            data.C_is= convert(Float64,bin2dec_twoscomp(_C_is))*2^-29

            #i0 inclination angle at reference time
            intermediate._i_0=reverse(word[7:14])
        end

    end

    if parameters.word_count==6

        if parameters.subframe_count==2
            #Cuc Cosine Correction to the Argument of Latitude
            _Cuc= reverse(word[15:30])
            data.C_uc= convert(Float64,bin2dec_twoscomp(_Cuc))*2^-29
            #e Eccentricity
            intermediate._e =reverse(word[7:14])
        end

        if parameters.subframe_count==3
            #i0 inclination angle at reference time
            intermediate._i_0= vcat(intermediate._i_0,reverse(word[7:30]))
            data.i_0= convert(Float64,bin2dec_twoscomp(intermediate._i_0))*2^-31*π
        end

    end

    if parameters.word_count==7

        if parameters.subframe_count==1
        groupDelayDifferential = reverse(word[7:14])
        data.groupDelayDifferential= convert(Float64,bin2dec_twoscomp(groupDelayDifferential))*2^-31
        end

        if parameters.subframe_count==2
            #e Eccentricity
            intermediate._e = vcat(intermediate._e,reverse(word[7:30]))
            data.e= convert(Float64,bin2dec(intermediate._e))*2^-33

        end

        if parameters.subframe_count==3
            #Crc cosine correction of orbit radius
            _C_rc= reverse(word[15:30])
            data.C_rc= convert(Float64,bin2dec_twoscomp(_C_rc))*2^-5
            #omega Argument of perigee
            intermediate._ω=reverse(word[7:14])
        end

    end

    if parameters.word_count==8

        if parameters.subframe_count==1
            IODC_2 = reverse(word[23:30])
            data.IODC=data.IODC*bitArray2Str(IODC_2)
            #SV Clock correction toc
            toc = reverse(word[7:22])
            data.toc = convert(Float64,bin2dec(toc))*2^4
        end

        if parameters.subframe_count==2
            #Cus Sine Correction to the Argument of Latitude
            _Cus= reverse(word[15:30])
            data.C_us= convert(Float64,bin2dec_twoscomp(_Cus))*2^-29

            #sqrtA Square Root of the Semi-Major Axis
            intermediate._sqrt_A=reverse(word[7:14])
        end

        if parameters.subframe_count==3
            #omega Argument of perigee
            intermediate._ω= vcat(intermediate._ω,reverse(word[7:30]))
            data.ω= convert(Float64,bin2dec_twoscomp(intermediate._ω))*2^-31*π
        end

    end

    if parameters.word_count==9

        if parameters.subframe_count==1
        #SV Clock Correction af2
        af2 = reverse(word[23:30])
        data.af2 = convert(Float64,bin2dec_twoscomp(af2))*2^-55
        #SV Clock Correction af1
        af1 = reverse(word[7:22])
        data.af1 = convert(Float64,bin2dec_twoscomp(af1))*2^-43
        end

        if parameters.subframe_count==2
            #sqrtA Square Root of the Semi-Major Axis
            intermediate._sqrt_A=vcat(intermediate._sqrt_A,reverse(word[7:30]))
            data.sqrt_A= convert(Float64,bin2dec(intermediate._sqrt_A))*2^-19
        end

        if parameters.subframe_count==3
            #omega dot is the rate of right ascension
            𝛀_dot = reverse(word[7:30])
            data.𝛀_dot = convert(Float64,bin2dec_twoscomp(𝛀_dot))*2^-43*π
        end

    end

    if parameters.word_count==10

        if parameters.subframe_count==1
            #SV Clock correction af0
            af0 = reverse(word[9:30])
            data.af0 = convert(Float64,bin2dec_twoscomp(af0))*2^-31
        end

        if parameters.subframe_count==2
            data.fitinterval = word[14]
            #toe Reference Time Ephemeris
            _t_oe = reverse(word[15:30])
            data.t_oe= convert(Float64,bin2dec(_t_oe))*2^4
            #validity time for the NMCT data in subframe 4
            aodo = reverse(word[9:13])
            data.aodo = convert(Float64,bin2dec(aodo))
        end

        if parameters.subframe_count==3
            #GPSData.IODE = reverse(word[23:30])

            #IDOT Rate of Inclinations angle
            _IDOT= reverse(word[9:22])
            data.IDOT= convert(Float64,bin2dec_twoscomp(_IDOT))*2^-43*π
        end

    parameters.word_count = 0

    end #end  del if word_count==10
    #return data,parameters,intermediate
end
