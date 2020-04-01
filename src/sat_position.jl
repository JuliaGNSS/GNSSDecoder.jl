
struct MissingDataError <: Exception
    text
end




"""
Compute Satellite position
$(SIGNATURES)

´dc´: Decoder witch all data collected
´code_phase´: code_phase of code for precise time and position computing

# Details
This functions provide the satellite's antenna phase centre position in the WGS-84 ECEF 
#
""" 
function sat_position(dc::GNSSDecoderState, code_phase::Float64 = 0.0)
    
    can_get_sat_position(dc) || throw(MissingDataError("Subframes 1-3 are not completely decoded"))
    code_time = code_phase / GNSSSignals.get_code_frequency(GPSL1) * Hz
    t = get_sat_time(dc) + code_time
    # t is the satellite time of transmission  corrected for transit
    OmegadotE = 7.2921151467e-5

    A  = (dc.data.sqrt_A)^2


    dc.data.tot = t - dc.data.t_oe

    #= dc.data.tot shall be the actual total time difference between the time tc and the epoch time toe,
    and must account for beginning or end of week crossovers =#
    if dc.data.tot > 302400
        dc.data.tot = dc.data.tot - 604800
    elseif dc.data.tot < -302400
        dc.data.tot = dc.data.tot + 604800
    end

    # newton method for eccentric anomaly
    Ek = eccentric_anomaly(dc)
    dc.data.tot = correct_sat_time(dc.data, Ek)
    Ek = eccentric_anomaly(dc)
    


    vk = atan(sqrt(1 - dc.data.e^2) * sin(Ek), (cos(Ek) - dc.data.e))
    # Argument of latitude
    Phik = vk + dc.data.ω

    # Second Harmonic Pertubations
    dUk = dc.data.C_us * sin(2 * Phik) + dc.data.C_uc * cos(2 * Phik)
    dRk = dc.data.C_rs * sin(2 * Phik) + dc.data.C_rc * cos(2 * Phik)
    dIk = dc.data.C_is * sin(2 * Phik) + dc.data.C_ic * cos(2 * Phik)

    # Corrected argument of latitude, radius and inclination
    Uk = Phik + dUk
    rk = A * (1 - dc.data.e * cos(Ek)) + dRk
    ik = dc.data.i_0 + dIk + dc.data.IDOT * dc.data.tot

    # Orbital Positions including corrections
    xk0 = rk * cos(Uk)
    yk0 = rk * sin(Uk)

    # Corrected longitude of ascending node
    Omegak = dc.data.𝛀_0 + (dc.data.𝛀_dot - OmegadotE) * dc.data.tot - OmegadotE * dc.data.t_oe

    # ECEF coordinates including earth rotation
    x = xk0 * cos(Omegak) - yk0 * cos(ik) * sin(Omegak)
    y = xk0 * sin(Omegak) + yk0 * cos(ik) * cos(Omegak)
    z = yk0 * sin(ik)

    position = SVector(x, y, z)
    return position

end

"""
Computes the eccentric anomaly
$SIGNATURES

´dc´ : Decoder with completely decoded data

# Details
Computes the eccentric anomaly using the decoded data, time of transmission and Newton-Raphson method
#
"""
function eccentric_anomaly(dc::GNSSDecoderState)
    A = dc.data.sqrt_A^2
    mu = 3.986005e14
    n0 = sqrt(mu / A^3)


    # Corrected mean motion and Mean anomaly at observation time
    n  = n0 + dc.data.Δn
    Mk = dc.data.M0 + n * dc.data.tot

    # newton method for eccentric anomaly
    Ek = newtonRaphson(dc.data.e, Mk)
    return Ek
end

"""
Checks if satellite position computation is possible
$(SIGNATURES)


´dc´ : Decoder

# Details
Checks if subframes 1-3 are decoded properly, Those contain all neccessary ephemeris data
#
"""
function can_get_sat_position(dc::GNSSDecoderState)
    if dc.subframes_decoded[1:3] == trues(3)
        return true
    else
        return false
    end
end

function newtonRaphson(ecc, Mk)
    E0 = Mk
    E1 = Mk + ecc * sin(E0)
    while abs(E1 - E0) > 10.0^-12
        E0 = E1
        E2 = Mk + ecc * sin(E1)
        E1 = E2
    end
    Ek = E1
    return Ek
end

"""
Compute corrected Satellite Time
$(SIGNATURES)


´tTX´ : satellite time provided in GPS-Seconds, not corrected! (Time of Data)
´eph´: ephemeris data
´Ek´: Eccentric Anomaly (Computed with ephemeris data in ´sat_position()´)

# Details
Corrects satellite time bias with provided satellite data 
#
"""
function correct_sat_time(data::GPSData, Ek)

    #Parameters
    mu = 3.986005e14    
    c  = 299792458
    F  = -2*sqrt(mu)/c^2
   
    #calculate satellite clock offset dtTX
    dtTX = data.tot - data.toc
    
    #relativistic correction term 
    dtr = F * data.e * data.sqrt_A * sin(Ek)
   
   #satellite clock error
    tTXerror = data.af0 + data.af1*dtTX + data.af2*dtTX^2 + dtr - data.groupDelayDifferential
    
   return totc = data.tot - tTXerror
end


"""
Computes Satellite Time
$(SIGNATURES)


´dc´ : Decoder, decoded at least first 3 subframes

# Details
Computes time of data, using the tow and bit count. t = tow*6 + num_bits*time_per_bit
#
"""
function get_sat_time(dc::GNSSDecoderState)
    tow_count = dc.data.tow * 6

    t_bits = dc.num_bits_buffered / get_data_frequency(GPSL1)*Hz
    return phase = tow_count + t_bits
end
    
