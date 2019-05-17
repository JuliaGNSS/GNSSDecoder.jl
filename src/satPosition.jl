#=
This functions provide the satellite's antenna phase centre position in the WGS-84 ECEF
=#
function satPosition(data, tc)

    mu = 3.986005e14
    OmegadotE = 7.2921151467e-5

    A  = (data.sqrt_A)^2
    n0 = sqrt(mu/A^3)

    tk = tc - data.t_oe

    #=tk shall be the actual total time difference between the time tc and the epoch time toe,
    and must account for beginning or end of week crossovers=#
    if tk > 302400
        tk = tk-604800
    elseif tk < -302400
        tk = tk+604800
    end

    #Corrected mean motion and Mean anomaly at observation time
    n  = n0 + data.Δn
    Mk = data.M0 + n*tk

    #newton method for eccentric anomaly
    Ek = newtonRaphson(data.e,Mk)


    vk = atan(sqrt(1-data.e^2) * sin(Ek) , (cos(Ek)-data.e))
    #Argument of latitude
    Phik = vk + data.ω

    #Second Harmonic Pertubations
    dUk = data.C_us*sin(2 * Phik) + data.C_uc*cos(2 * Phik)
    dRk = data.C_rs*sin(2 * Phik) + data.C_rc*cos(2 * Phik)
    dIk = data.C_is*sin(2 * Phik) + data.C_ic*cos(2 * Phik)

    #Corrected argument of latitude, radius and inclination
    Uk = Phik + dUk
    rk = A*(1-data.e*cos(Ek)) + dRk
    ik = data.i_0 + dIk + data.IDOT*tk

    #Orbital Positions including corrections
    xk0 = rk*cos(Uk)
    yk0 = rk*sin(Uk)

    #Corrected longitude of ascending node
    Omegak = data.𝛀_0 + (data.𝛀_dot-OmegadotE)*tk - OmegadotE*data.t_oe

    #ECEF coordinates including earth rotation
    x = xk0*cos(Omegak) - yk0*cos(ik)*sin(Omegak)
    y = xk0*sin(Omegak) + yk0*cos(ik)*cos(Omegak)
    z = yk0*sin(ik)
    return x,y,z,Ek
end


function newtonRaphson(ecc, Mk)
    E0 = Mk
    E1 = Mk + ecc*sin(E0)
    while abs(E1-E0)>10.0^-12
       E0 = E1
       E2 = Mk + ecc*sin(E1)
       E1 = E2
   end
   Ek = E1
   return Ek
end
