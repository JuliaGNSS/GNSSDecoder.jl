
struct TLM_HOW_Data_Struct
    integrity_status_flag::Bool
    TOW::Int64
    alert_flag::Bool
    anti_spoof_flag::Bool
end


struct Subframe_1_Data
    trans_week::Int64
    codeonl2::Int64
    ura::Float64
    svhealth::String
    IODC::String
    l2pcode::Bool
    T_GD::Float64
    t_oc::Int64
    a_f2::Float64
    a_f1::Float64
    a_f0::Float64
end

struct Subframe_2_Data
    IODE::String
    C_rs::Float64
    Δn::Float64
    M_0::Float64
    C_uc::Float64
    e::Float64
    C_us::Float64
    sqrt_A::Float64
    t_oe::Int64
    fit_interval::Bool
    AODO::Int64
end

struct Subframe_3_Data
    C_ic::Float64
    Ω_0::Float64
    C_is::Float64
    i_0::Float64
    C_rc::Float64
    ω::Float64
    Ω_dot::Float64
    IODE::String
    IDOT::Float64
end

struct Subframe_4_Data
end

struct Subframe_5_Data
end




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

function GPSData(
    prev_data::GPSData,
    subfr_1_data::Subframe_1_Data)
    

    data = GPSData(
        prev_data,

        trans_week = subfr_1_data.trans_week,
        codeonl2 = subfr_1_data.codeonl2,
        ura = subfr_1_data.ura,
        svhealth = subfr_1_data.svhealth,
        IODC = subfr_1_data.IODC,
        l2pcode = subfr_1_data.l2pcode,
        T_GD = subfr_1_data.T_GD,
        t_oc = subfr_1_data.t_oc,
        a_f2 = subfr_1_data.a_f2,
        a_f1 = subfr_1_data.a_f1,
        a_f0 = subfr_1_data.a_f0
    )
    return data
end

function GPSData(
    prev_data::GPSData,
    subfr_2_data::Subframe_2_Data)
    

    data = GPSData(
        prev_data,

        IODE_Sub_2 = subfr_2_data.IODE,
        C_rs = subfr_2_data.C_rs,
        Δn = subfr_2_data.Δn,
        M_0  = subfr_2_data.M_0,
        C_uc = subfr_2_data.C_uc,
        e = subfr_2_data.e,
        C_us = subfr_2_data.C_us,
        sqrt_A = subfr_2_data.sqrt_A,
        t_oe = subfr_2_data.t_oe,
        fit_interval = subfr_2_data.fit_interval,
        AODO = subfr_2_data.AODO
    )
    return data
end

function GPSData(
    prev_data::GPSData,
    subfr_3_data::Subframe_3_Data)
    

    data = GPSData(
        prev_data,

        C_ic = subfr_3_data.C_ic,
        Ω_0 = subfr_3_data.Ω_0,
        C_is = subfr_3_data.C_is,
        i_0 = subfr_3_data.i_0,
        C_rc = subfr_3_data.C_rc,
        ω = subfr_3_data.ω,
        Ω_dot = subfr_3_data.Ω_dot,
        IODE_Sub_3 = subfr_3_data.IODE,
        IDOT = subfr_3_data.IDOT
    )
    return data
end

function GPSData(
    prev_data::GPSData,
    TLM_HOW_Data::TLM_HOW_Data_Struct)
    
    data = GPSData(
        prev_data,
        
        integrity_status_flag = TLM_HOW_Data.integrity_status_flag,
        TOW = TLM_HOW_Data.TOW,
        alert_flag = TLM_HOW_Data.alert_flag,
        anti_spoof_flag = TLM_HOW_Data.anti_spoof_flag,
    )
    return data
end



@with_kw struct GPSL1Constants
    PI::Float64 = 3.1415926535898
    Ω_dot_e::Float64 = 7.2921151467e-5 
    c::Float64 = 2.99792458e8
    μ::Float64 = 3.986005e14    
    F::Float64 = -4.442807633e-10
end


@with_kw mutable struct GNSSDecoderState
    PRN::Int64
    buffer::BitArray{1} = falses(BUFFER_LENGTH)
    data::GPSData = GPSData()
    constants::GPSL1Constants = GPSL1Constants()
    prev_30::Bool = false
    prev_29::Bool = false
    data_integrity::Bool = true # pass of parity check
    subframes_decoded::MVector{5,Bool} = MVector{5,Bool}(false, false, false, false, false)
    num_bits_buffered::Int = 0
end