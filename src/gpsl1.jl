# UInt320 buffer for GPS L1 (more efficient than UInt312)
# which holds at least a complete GPS L1 subframe
BitIntegers.@define_integers 320

Base.@kwdef struct GPSL1Constants <: AbstractGNSSConstants
    subframe_length::Int = 300
    preamble::UInt8 = 0b10001011
    preamble_length::Int = 8
    word_length::Int = 30
    PI::Float64 = 3.1415926535898
    Ω_dot_e::Float64 = 7.2921151467e-5 
    c::Float64 = 2.99792458e8
    μ::Float64 = 3.986005e14    
    F::Float64 = -4.442807633e-10
end

Base.@kwdef struct GPSL1Data <: AbstractGNSSData
    last_subframe_id::Int = 0
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

function GPSL1Data(
    data::GPSL1Data;
    last_subframe_id = data.last_subframe_id,
    integrity_status_flag = data.integrity_status_flag,
    TOW = data.TOW,
    alert_flag = data.alert_flag,
    anti_spoof_flag = data.anti_spoof_flag,
    trans_week = data.trans_week,
    codeonl2 = data.codeonl2,
    ura = data.ura,
    svhealth = data.svhealth,
    IODC = data.IODC,
    l2pcode = data.l2pcode,
    T_GD = data.T_GD,
    t_oc = data.t_oc,
    a_f2 = data.a_f2,
    a_f1 = data.a_f1,
    a_f0 = data.a_f0,
    IODE_Sub_2 = data.IODE_Sub_2,
    C_rs = data.C_rs,
    Δn = data.Δn,
    M_0 = data.M_0,
    C_uc = data.C_uc,
    e = data.e,
    C_us = data.C_us,
    sqrt_A = data.sqrt_A,
    t_oe = data.t_oe,
    fit_interval = data.fit_interval,
    AODO = data.AODO,
    C_ic = data.C_ic,
    Ω_0 = data.Ω_0,
    C_is = data.C_is,
    i_0 = data.i_0,
    C_rc = data.C_rc,
    ω = data.ω,
    Ω_dot = data.Ω_dot,
    IODE_Sub_3 = data.IODE_Sub_3,
    IDOT = data.IDOT,
)
    GPSL1Data(
        last_subframe_id,
        integrity_status_flag,
        TOW,
        alert_flag,
        anti_spoof_flag,
        trans_week,
        codeonl2,
        ura,
        svhealth,
        IODC,
        l2pcode,
        T_GD,
        t_oc,
        a_f2,
        a_f1,
        a_f0,
        IODE_Sub_2,
        C_rs,
        Δn,
        M_0,
        C_uc,
        e,
        C_us,
        sqrt_A,
        t_oe,
        fit_interval,
        AODO,
        C_ic,
        Ω_0,
        C_is,
        i_0,
        C_rc,
        ω,
        Ω_dot,
        IODE_Sub_3,
        IDOT,
    )
end

function is_subframe1_decoded(data::GPSL1Data)
    !isnothing(data.trans_week) &&
        !isnothing(data.codeonl2) &&
        !isnothing(data.ura) &&
        !isnothing(data.svhealth) &&
        !isnothing(data.IODC) &&
        !isnothing(data.l2pcode) &&
        !isnothing(data.T_GD) &&
        !isnothing(data.t_oc) &&
        !isnothing(data.a_f2) &&
        !isnothing(data.a_f1) &&
        !isnothing(data.a_f0)
end

function is_subframe2_decoded(data::GPSL1Data)
    !isnothing(data.IODE_Sub_2) &&
        !isnothing(data.C_rs) &&
        !isnothing(data.Δn) &&
        !isnothing(data.M_0) &&
        !isnothing(data.C_uc) &&
        !isnothing(data.e) &&
        !isnothing(data.C_us) &&
        !isnothing(data.sqrt_A) &&
        !isnothing(data.t_oe) &&
        !isnothing(data.fit_interval) &&
        !isnothing(data.AODO)
end

function is_subframe3_decoded(data::GPSL1Data)
    !isnothing(data.C_ic) &&
        !isnothing(data.Ω_0) &&
        !isnothing(data.C_is) &&
        !isnothing(data.i_0) &&
        !isnothing(data.C_rc) &&
        !isnothing(data.ω) &&
        !isnothing(data.Ω_dot) &&
        !isnothing(data.IODE_Sub_3) &&
        !isnothing(data.IDOT)
end

function is_subframe4_decoded(data::GPSL1Data)
    false
end

function is_subframe5_decoded(data::GPSL1Data)
    false
end

function is_decoding_completed_for_positioning(data::GPSL1Data)
    !isnothing(data.integrity_status_flag) &&
        !isnothing(data.TOW) &&
        !isnothing(data.alert_flag) &&
        !isnothing(data.anti_spoof_flag) &&
        is_subframe1_decoded(data) &&
        is_subframe2_decoded(data) &&
        is_subframe3_decoded(data)
end

function GPSL1DecoderState(prn)
    GNSSDecoderState(prn, UInt320(0), UInt320(0), GPSL1Data(), GPSL1Data(), GPSL1Constants(), 0, nothing, false)
end

function check_gpsl1_parity(word::Unsigned, prev_29 = false, prev_30 = false)
    function bit(bit_number)
        cbit = get_bit(word, 30, bit_number)
        prev_30 ? !cbit : cbit
    end
    # Parity check to verify the data integrity:
    D_25 = prev_29 ⊻ bit(1) ⊻ bit(2) ⊻ bit(3) ⊻ bit(5) ⊻ bit(6) ⊻ bit(10) ⊻ bit(11) ⊻ bit(12) ⊻ bit(13) ⊻ bit(14) ⊻ bit(17) ⊻ bit(18) ⊻ bit(20) ⊻ bit(23)
    D_26 = prev_30 ⊻ bit(2) ⊻ bit(3) ⊻ bit(4) ⊻ bit(6) ⊻ bit(7) ⊻ bit(11) ⊻ bit(12) ⊻ bit(13) ⊻ bit(14) ⊻ bit(15) ⊻ bit(18) ⊻ bit(19) ⊻ bit(21) ⊻ bit(24)
    D_27 = prev_29 ⊻ bit(1) ⊻ bit(3) ⊻ bit(4) ⊻ bit(5) ⊻ bit(7) ⊻ bit(8) ⊻ bit(12) ⊻ bit(13) ⊻ bit(14) ⊻ bit(15) ⊻ bit(16) ⊻ bit(19) ⊻ bit(20) ⊻ bit(22)
    D_28 = prev_30 ⊻ bit(2) ⊻ bit(4) ⊻ bit(5) ⊻ bit(6) ⊻ bit(8) ⊻ bit(9) ⊻ bit(13) ⊻ bit(14) ⊻ bit(15) ⊻ bit(16) ⊻ bit(17) ⊻ bit(20) ⊻ bit(21) ⊻ bit(23)
    D_29 = prev_30 ⊻ bit(1) ⊻ bit(3) ⊻ bit(5) ⊻ bit(6) ⊻ bit(7) ⊻ bit(9) ⊻ bit(10) ⊻ bit(14) ⊻ bit(15) ⊻ bit(16) ⊻ bit(17) ⊻ bit(18) ⊻ bit(21) ⊻ bit(22) ⊻ bit(24)
    D_30 = prev_29 ⊻ bit(3) ⊻ bit(5) ⊻ bit(6) ⊻ bit(8) ⊻ bit(9) ⊻ bit(10) ⊻ bit(11) ⊻ bit(13) ⊻ bit(15) ⊻ bit(19) ⊻ bit(22) ⊻ bit(23) ⊻ bit(24)
    computed_parity_bits = ((((D_25 << UInt(1) + D_26) << UInt(1) + D_27) << UInt(1) + D_28) << UInt(1) + D_29) << UInt(1) + D_30
    computed_parity_bits == get_bits(word, 30, 25, 6)
end

function can_decode_word(decode_bits::Function, state::GNSSDecoderState, word_number::Int)
    word = get_word(state, word_number)
    prev_word = get_word(state, word_number - 1)
    prev_word_bit_30 = get_bit(prev_word, 30, 30)
    is_checked_word = word_number != 1 && word_number != 3
    if check_gpsl1_parity(word, is_checked_word * get_bit(prev_word, 30, 29), is_checked_word * prev_word_bit_30)
        word_comp = (is_checked_word * prev_word_bit_30) ? ~word : word
        data = decode_bits(word_comp, state)
        state = GNSSDecoderState(state, raw_data = data)
    end
    return state
end

function can_decode_two_words(decode_bits::Function, state::GNSSDecoderState, word_number1::Int, word_number2::Int)
    word1 = get_word(state, word_number1)
    prev_word1 = get_word(state, word_number1 - 1)
    prev_word1_bit_30 = get_bit(prev_word1, 30, 30)
    is_checked_word1 = word_number1 != 1 && word_number1 != 3
    word2 = get_word(state, word_number2)
    prev_word2 = get_word(state, word_number2 - 1)
    prev_word2_bit_30 = get_bit(prev_word2, 30, 30)
    is_checked_word2 = word_number2 != 1 && word_number2 != 3
    if check_gpsl1_parity(word1, is_checked_word1 * get_bit(prev_word1, 30, 29), is_checked_word1 * prev_word1_bit_30) &&
        check_gpsl1_parity(word2, is_checked_word2 * get_bit(prev_word2, 30, 29), is_checked_word2 * prev_word2_bit_30)
        word1_comp = (is_checked_word1 * prev_word1_bit_30) ? ~word1 : word1
        word2_comp = (is_checked_word2 * prev_word2_bit_30) ? ~word2 : word2
        data = decode_bits(word1_comp, word2_comp, state)
        state = GNSSDecoderState(state, raw_data = data)
    end
    return state
end

function read_tlm_and_how_words(state)
    state = can_decode_word(state, 1) do tlm_word, state
        integrity_status_flag = get_bit(tlm_word, 30, 23)
        GPSL1Data(state.raw_data; integrity_status_flag)
    end
    prev_TOW = state.raw_data.TOW
    state = can_decode_word(state, 2) do how_word, state
        TOW = get_bits(how_word, 30, 1, 17)
        alert_flag = get_bit(how_word, 30, 18)
        anti_spoof_flag = get_bit(how_word, 30, 19)
        last_subframe_id = get_bits(how_word, 30, 20, 3)
        GPSL1Data(
            state.raw_data;
            last_subframe_id,
            TOW,
            alert_flag,
            anti_spoof_flag
        )
    end
    if !isnothing(prev_TOW) && prev_TOW + 1 != state.raw_data.TOW
        # Time of week must be decodable
        state = GNSSDecoderState(state, raw_data = GPSL1Data(
            state.raw_data;
            TOW = nothing
        ))
    end
    state
end

function decode_frame(state::GNSSDecoderState{<:GPSL1Data})
    state = read_tlm_and_how_words(state)
    subframe_id = state.raw_data.last_subframe_id

    if subframe_id == 1
        state = can_decode_word(state, 3) do word3, state
            trans_week = get_bits(word3, 30, 1, 10)

            # Codes on L2 Channel
            codeonl2 = get_bits(word3, 30, 11, 2)

            # SV Accuracy, user range accuracy
            ura  = get_bits(word3, 30, 13, 4)
            if ura <= 6
                ura = 2^(1 + (ura / 2))
            elseif 6 < ura <= 14
                ura = 2^(ura - 2)
            elseif ura == 15
                ura = nothing
            end

            
            # Satellite Health
            svhealth = bitstring(get_bits(word3, 30, 17, 6))[end - 5:end]
            if get_bit(word3, 30, 17)
                @warn "Bad LNAV Data, SV-Health critical", svhealth
            end

            GPSL1Data(state.raw_data; trans_week, codeonl2, ura, svhealth)
        end


        state = can_decode_two_words(state, 3, 8) do word3, word8, state
            # Issue of Data Clock
            # 2 MSB in Word 2, LSB 8 in Word 8
            IODC = bitstring(get_bits(word3, 30, 23, 2))[end - 1:end] * bitstring(get_bits(word8, 30, 1, 8))[end - 7:end]
            GPSL1Data(state.raw_data; IODC)
        end

        state = can_decode_word(state, 4) do word4, state
            # True: LNAV Datastream on PCode commanded OFF
            l2pcode = get_bit(word4, 30, 1)
            GPSL1Data(state.raw_data; l2pcode)
        end

        state = can_decode_word(state, 7) do word7, state
            # group time differential
            T_GD = get_two_complement_num(word7, 30, 17, 8) / 1 << 31
            GPSL1Data(state.raw_data; T_GD)
        end

        state = can_decode_word(state, 8) do word8, state
            # Clock data reference
            t_oc = get_bits(word8, 30, 9, 16) << 4
            GPSL1Data(state.raw_data; t_oc)
        end

        state = can_decode_word(state, 9) do word9, state
            # clock correction parameter a_f2
            a_f2 = get_two_complement_num(word9, 30, 1, 8) / 1 << 55
            # clock correction parameter a_f1
            a_f1 = get_two_complement_num(word9, 30, 9, 16) / 1 << 43
            GPSL1Data(state.raw_data; a_f2, a_f1)
        end

        state = can_decode_word(state, 10) do word10, state
            # Clock data reference
            a_f0 = get_two_complement_num(word10, 30, 1, 22) / 1 << 31
            GPSL1Data(state.raw_data; a_f0)
        end
    elseif subframe_id == 2
        state = can_decode_word(state, 3) do word3, state
            # Issue of ephemeris data
            IODE_Sub_2 = bitstring(get_bits(word3, 30, 1, 8))[end-7:end]
            C_rs = get_two_complement_num(word3, 30, 9, 16) / 1 << 5
            GPSL1Data(state.raw_data; IODE_Sub_2, C_rs)
        end

        state = can_decode_word(state, 4) do word4, state
            # Mean motion difference from computed value
            Δn = get_two_complement_num(word4, 30, 1, 16) * state.constants.PI / 1 << 43
            GPSL1Data(state.raw_data; Δn)
        end

        state = can_decode_two_words(state, 4, 5) do word4, word5, state
            # Mean motion difference from computed value
            combined_word = UInt(get_bits(word4, 30, 17, 8) << 24 + get_bits(word5, 30, 1, 24))
            M_0 = get_two_complement_num(combined_word, 32, 1, 32) * state.constants.PI / 1 << 31
            GPSL1Data(state.raw_data; M_0)
        end

        state = can_decode_word(state, 6) do word6, state
            # Amplitude of the Cosine Harmonic Correction Term to the Argument Latitude
            C_uc = get_two_complement_num(word6, 30, 1, 16) / 1 << 29
            GPSL1Data(state.raw_data; C_uc)
        end
        
        state = can_decode_two_words(state, 6, 7) do word6, word7, state
            # Eccentricity
            e = (get_bits(word6, 30, 17, 8) << 24 + get_bits(word7, 30, 1, 24)) / 1 << 33
            GPSL1Data(state.raw_data; e)
        end

        state = can_decode_word(state, 8) do word8, state
            # Amplitude of the Sine Harmonic Correction Term to the Argument of Latitude
            C_us = get_two_complement_num(word8, 30, 1, 16) / 1 << 29
            GPSL1Data(state.raw_data; C_us)
        end

        state = can_decode_two_words(state, 8, 9) do word8, word9, state
            # Square Root of Semi-Major Axis
            sqrt_A = (get_bits(word8, 30, 17, 8) << 24 + get_bits(word9, 30, 1, 24)) / 1 << 19
            GPSL1Data(state.raw_data; sqrt_A)
        end

        state = can_decode_word(state, 10) do word10, state
            # Reference Time ephemeris
            t_oe = get_bits(word10, 30, 1, 16) << 4
            fit_interval = get_bit(word10, 30, 17)
            AODO = get_bits(word10, 30, 18, 5)
            GPSL1Data(state.raw_data; t_oe, fit_interval, AODO)
        end
    elseif subframe_id == 3
        state = can_decode_word(state, 3) do word3, state
            # Amplitude of the Cosine Harmonic Correction to Angle of Inclination
            C_ic = get_two_complement_num(word3, 30, 1, 16) / 1 << 29
            GPSL1Data(state.raw_data; C_ic)
        end

        state = can_decode_two_words(state, 3, 4) do word3, word4, state
            # Longitude of Ascending Node of Orbit Plane at Weekly Epoch
            combined_word = UInt(get_bits(word3, 30, 17, 8) << 24 + get_bits(word4, 30, 1, 24))
            Ω_0 = get_two_complement_num(combined_word, 32, 1, 32) * state.constants.PI / 1 << 31
            GPSL1Data(state.raw_data; Ω_0)
        end

        state = can_decode_word(state, 5) do word5, state
            # Amplitude of the sine harmonic correction term to angle of Inclination
            C_is = get_two_complement_num(word5, 30, 1, 16) / 1 << 29
            GPSL1Data(state.raw_data; C_is)
        end

        state = can_decode_two_words(state, 5, 6) do word5, word6, state
            # inclination Angle at reference time
            combined_word = UInt(get_bits(word5, 30, 17, 8) << 24 + get_bits(word6, 30, 1, 24))
            i_0 = get_two_complement_num(combined_word, 32, 1, 32) * state.constants.PI / 1 << 31
            GPSL1Data(state.raw_data; i_0)
        end

        state = can_decode_word(state, 7) do word7, state
            # Amplitude of the cosine harmonic correction term to orbit Radius
            C_rc = get_two_complement_num(word7, 30, 1, 16) / 1 << 5
            GPSL1Data(state.raw_data; C_rc)
        end

        state = can_decode_two_words(state, 7, 8) do word7, word8, state
            # Argument of Perigee
            combined_word = UInt(get_bits(word7, 30, 17, 8) << 24 + get_bits(word8, 30, 1, 24))
            ω = get_two_complement_num(combined_word, 32, 1, 32) * state.constants.PI / 1 << 31
            GPSL1Data(state.raw_data; ω)
        end

        state = can_decode_word(state, 9) do word9, state
            # Amplitude of the cosine harmonic correction term to orbit Radius
            Ω_dot = get_two_complement_num(word9, 30, 1, 24) * state.constants.PI / 1 << 43
            GPSL1Data(state.raw_data; Ω_dot)
        end

        state = can_decode_word(state, 10) do word10, state
            # Issue of Ephemeris Data
            IODE_Sub_3 = bitstring(get_bits(word10, 30, 1, 8))[end-7:end]
            # Rate of Inclination Angle
            IDOT = get_two_complement_num(word10, 30, 9, 14) * state.constants.PI / 1 << 43
            GPSL1Data(state.raw_data; IODE_Sub_3, IDOT)
        end
    end

    return state
end

function validate_data(state::GNSSDecoderState{<:GPSL1Data})
    if is_decoding_completed_for_positioning(state.raw_data) &&
        state.raw_data.IODC[3:10] == state.raw_data.IODE_Sub_2 == state.raw_data.IODE_Sub_3
        state = GNSSDecoderState(state, data = state.raw_data, num_bits_after_valid_subframe = 8)
    end
    return state
end 