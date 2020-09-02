
    @with_kw struct GPSData
        integrity_status_flag::Union{Nothing,Bool} = nothing
        TOW::Union{Nothing,Int64} = nothing
        alert_flag::Union{Nothing,Bool} = nothing
        anti_spoof_flag::Union{Nothing,Bool} = nothing

        trans_week::Union{Nothing,Int64} = nothing
        codeonl2::Union{Nothing,Int64} = nothing
        ura::Union{Nothing,Float64} = nothing
        svhealth::Union{Nothing,String} = nothing
        IODC::Union{Nothing,String} = nothing
        l2pcode::Union{Nothing,Bool} = nothing
        T_GD::Union{Nothing,Float64} = nothing
        t_oc::Union{Nothing,Int64} = nothing
        a_f2::Union{Nothing,Float64} = nothing
        a_f1::Union{Nothing,Float64} = nothing
        a_f0::Union{Nothing,Float64} = nothing

        IODE_Sub_2::Union{Nothing,String} = nothing
        C_rs::Union{Nothing,Float64} = nothing
        Δn::Union{Nothing,Float64} = nothing
        M_0::Union{Nothing,Float64} = nothing
        C_uc::Union{Nothing,Float64} = nothing
        e::Union{Nothing,Float64} = nothing
        C_us::Union{Nothing,Float64} = nothing
        sqrt_A::Union{Nothing,Float64} = nothing
        t_oe::Union{Nothing,Int64} = nothing
        fit_interval::Union{Nothing,Bool} = nothing
        AODO::Union{Nothing,Int64} = nothing

        C_ic::Union{Nothing,Float64} = nothing
        Ω_0::Union{Nothing,Float64} = nothing
        C_is::Union{Nothing,Float64} = nothing
        i_0::Union{Nothing,Float64} = nothing
        C_rc::Union{Nothing,Float64} = nothing
        ω::Union{Nothing,Float64} = nothing
        Ω_dot::Union{Nothing,Float64} = nothing
        IODE_Sub_3::Union{Nothing,String} = nothing
        IDOT::Union{Nothing,Float64} = nothing
    end 


    @with_kw struct GPSL1Constants
        PI = 3.1415926535898
        Ω_dot_e = 7.2921151467e-5 
        c = 2.99792458e8
        μ = 3.986005e14    
        F =  -4.442807633e-10
    end


    @with_kw mutable struct GNSSDecoderState
        PRN::Int64
        buffer::BitArray{1} = falses(BUFFER_LENGTH)
        data::GPSData = GPSData()
        constants::GPSL1Constants = GPSL1Constants()
        preamble_found::Bool = false
        subframe_count::Int = 0
        prev_30::Bool = false
        prev_29::Bool = false
        data_integrity::Bool = true # pass of parity check
        new_data_needed::Bool = false # IODE must be equal to the 8 LSB of the IODC. If they don't match new data must be collected
        subframes_decoded::MVector{5,Bool} = MVector{5,Bool}(false, false, false, false, false)
        nb_prev::Int = 0
        num_bits_buffered::Int = 0
    end