# UInt320 buffer for GPS L1 (more efficient than UInt312)
# which holds at least a complete GPS L1 subframe plus
# 8 extra syncronization bíts
BitIntegers.@define_integers 320

Base.@kwdef struct GPSL1Constants <: AbstractGNSSConstants
    syncro_sequence_length::Int = 300
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
    t_0c::Union{Nothing,Int64} = nothing
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
    t_0e::Union{Nothing,Int64} = nothing
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
    i_dot::Union{Nothing,Float64} = nothing

    # Subframe 4 page 18: Ionospheric parameters
    α_0::Union{Nothing,Float64} = nothing
    α_1::Union{Nothing,Float64} = nothing
    α_2::Union{Nothing,Float64} = nothing
    α_3::Union{Nothing,Float64} = nothing
    β_0::Union{Nothing,Float64} = nothing
    β_1::Union{Nothing,Float64} = nothing
    β_2::Union{Nothing,Float64} = nothing
    β_3::Union{Nothing,Float64} = nothing

    # Subframe 4 page 18: UTC parameters
    A_0::Union{Nothing,Float64} = nothing
    A_1::Union{Nothing,Float64} = nothing
    Δt_LS::Union{Nothing,Int64} = nothing
    t_ot::Union{Nothing,Int64} = nothing
    WN_t::Union{Nothing,Int64} = nothing
    WN_LSF::Union{Nothing,Int64} = nothing
    DN::Union{Nothing,Int64} = nothing
    Δt_LSF::Union{Nothing,Int64} = nothing

    # Subframe 4 page 25: A-S flags and SV configurations (32 SVs)
    sv_config::Union{Nothing,Vector{Int64}} = nothing

    # Subframe 4 page 25: SV health for SV 25-32 (6-bit health words)
    sv_health_sf4_25::Union{Nothing,Vector{String}} = nothing

    # Subframe 5 pages 1-24: Almanac data (stored per SV ID)
    almanac::Union{Nothing,Dict{Int64,NamedTuple}} = nothing

    # Subframe 5 page 25: SV health for SV 1-24 (6-bit health words)
    sv_health_sf5_25::Union{Nothing,Vector{String}} = nothing

    # Subframe 5 page 25: Almanac reference time and week
    t_oa::Union{Nothing,Int64} = nothing
    WN_a::Union{Nothing,Int64} = nothing
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
    t_0c = data.t_0c,
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
    t_0e = data.t_0e,
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
    i_dot = data.i_dot,
    α_0 = data.α_0,
    α_1 = data.α_1,
    α_2 = data.α_2,
    α_3 = data.α_3,
    β_0 = data.β_0,
    β_1 = data.β_1,
    β_2 = data.β_2,
    β_3 = data.β_3,
    A_0 = data.A_0,
    A_1 = data.A_1,
    Δt_LS = data.Δt_LS,
    t_ot = data.t_ot,
    WN_t = data.WN_t,
    WN_LSF = data.WN_LSF,
    DN = data.DN,
    Δt_LSF = data.Δt_LSF,
    sv_config = data.sv_config,
    sv_health_sf4_25 = data.sv_health_sf4_25,
    almanac = data.almanac,
    sv_health_sf5_25 = data.sv_health_sf5_25,
    t_oa = data.t_oa,
    WN_a = data.WN_a,
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
        t_0c,
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
        t_0e,
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
        i_dot,
        α_0,
        α_1,
        α_2,
        α_3,
        β_0,
        β_1,
        β_2,
        β_3,
        A_0,
        A_1,
        Δt_LS,
        t_ot,
        WN_t,
        WN_LSF,
        DN,
        Δt_LSF,
        sv_config,
        sv_health_sf4_25,
        almanac,
        sv_health_sf5_25,
        t_oa,
        WN_a,
    )
end

struct VotedGPSL1Data
    vote::Int
    data::GPSL1Data
end

struct GPSL1Cache <: AbstractGNSSCache
    old_data::Vector{VotedGPSL1Data}
end

function GPSL1Cache()
    GPSL1Cache(Vector{VotedGPSL1Data}())
end

function is_subframe1_decoded(data::GPSL1Data)
    !isnothing(data.trans_week) &&
        !isnothing(data.codeonl2) &&
        !isnothing(data.ura) &&
        !isnothing(data.svhealth) &&
        !isnothing(data.IODC) &&
        !isnothing(data.l2pcode) &&
        !isnothing(data.T_GD) &&
        !isnothing(data.t_0c) &&
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
        !isnothing(data.t_0e) &&
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
        !isnothing(data.i_dot)
end

function is_subframe4_decoded(data::GPSL1Data)
    # Subframe 4 is considered decoded when we have ionospheric and UTC parameters
    # (page 18) and SV configurations/health (page 25)
    !isnothing(data.α_0) &&
        !isnothing(data.α_1) &&
        !isnothing(data.α_2) &&
        !isnothing(data.α_3) &&
        !isnothing(data.β_0) &&
        !isnothing(data.β_1) &&
        !isnothing(data.β_2) &&
        !isnothing(data.β_3) &&
        !isnothing(data.A_0) &&
        !isnothing(data.A_1) &&
        !isnothing(data.Δt_LS) &&
        !isnothing(data.t_ot) &&
        !isnothing(data.WN_t) &&
        !isnothing(data.WN_LSF) &&
        !isnothing(data.DN) &&
        !isnothing(data.Δt_LSF) &&
        !isnothing(data.sv_config) &&
        !isnothing(data.sv_health_sf4_25)
end

function is_subframe5_decoded(data::GPSL1Data)
    # Subframe 5 is considered decoded when we have SV health (page 25)
    # and almanac reference time/week
    !isnothing(data.sv_health_sf5_25) &&
        !isnothing(data.t_oa) &&
        !isnothing(data.WN_a)
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
    GNSSDecoderState(
        prn,
        UInt320(0),
        UInt320(0),
        GPSL1Data(),
        GPSL1Data(),
        GPSL1Constants(),
        GPSL1Cache(),
        0,
        nothing,
        false,
    )
end

function GNSSDecoderState(system::GPSL1, prn)
    GNSSDecoderState(
        prn,
        UInt320(0),
        UInt320(0),
        GPSL1Data(),
        GPSL1Data(),
        GPSL1Constants(),
        GPSL1Cache(),
        0,
        nothing,
        false,
    )
end

function reset_decoder_state(state::GNSSDecoderState{<:GPSL1Data})
    # Reset bit buffers and TOW data field, while keeping the
    # remaining parameters in raw_data. This allows a GNSSReceiver
    # to use a satellite after a reacquisition without waiting for
    # the decoding of all data fields.
    # Note: trans_week is currently not reset as it is only
    # broadcast in subframe 1 and thus may increase the time until
    # the decoder is available again after an outage. This will
    # lead to erroneous decoder information for a few seconds after
    # reacquisition when a new GPS week started during a signal outage.
    GNSSDecoderState(
        state;
        raw_buffer = UInt320(0),
        buffer = UInt320(0),
        raw_data = GPSL1Data(state.raw_data; TOW = nothing),
        data = GPSL1Data(),
        num_bits_buffered = 0,
        num_bits_after_valid_syncro_sequence = nothing,
    )
end

function check_gpsl1_parity(word::Unsigned, prev_29 = false, prev_30 = false)
    function bit(bit_number)
        cbit = get_bit(word, 30, bit_number)
        prev_30 ? !cbit : cbit
    end
    # Parity check to verify the data integrity:
    D_25 =
        prev_29 ⊻ bit(1) ⊻ bit(2) ⊻ bit(3) ⊻ bit(5) ⊻ bit(6) ⊻ bit(10) ⊻ bit(11) ⊻ bit(12) ⊻
        bit(13) ⊻ bit(14) ⊻ bit(17) ⊻ bit(18) ⊻ bit(20) ⊻ bit(23)
    D_26 =
        prev_30 ⊻ bit(2) ⊻ bit(3) ⊻ bit(4) ⊻ bit(6) ⊻ bit(7) ⊻ bit(11) ⊻ bit(12) ⊻ bit(13) ⊻
        bit(14) ⊻ bit(15) ⊻ bit(18) ⊻ bit(19) ⊻ bit(21) ⊻ bit(24)
    D_27 =
        prev_29 ⊻ bit(1) ⊻ bit(3) ⊻ bit(4) ⊻ bit(5) ⊻ bit(7) ⊻ bit(8) ⊻ bit(12) ⊻ bit(13) ⊻
        bit(14) ⊻ bit(15) ⊻ bit(16) ⊻ bit(19) ⊻ bit(20) ⊻ bit(22)
    D_28 =
        prev_30 ⊻ bit(2) ⊻ bit(4) ⊻ bit(5) ⊻ bit(6) ⊻ bit(8) ⊻ bit(9) ⊻ bit(13) ⊻ bit(14) ⊻
        bit(15) ⊻ bit(16) ⊻ bit(17) ⊻ bit(20) ⊻ bit(21) ⊻ bit(23)
    D_29 =
        prev_30 ⊻ bit(1) ⊻ bit(3) ⊻ bit(5) ⊻ bit(6) ⊻ bit(7) ⊻ bit(9) ⊻ bit(10) ⊻ bit(14) ⊻
        bit(15) ⊻ bit(16) ⊻ bit(17) ⊻ bit(18) ⊻ bit(21) ⊻ bit(22) ⊻ bit(24)
    D_30 =
        prev_29 ⊻ bit(3) ⊻ bit(5) ⊻ bit(6) ⊻ bit(8) ⊻ bit(9) ⊻ bit(10) ⊻ bit(11) ⊻ bit(13) ⊻
        bit(15) ⊻ bit(19) ⊻ bit(22) ⊻ bit(23) ⊻ bit(24)
    computed_parity_bits =
        (
            (((D_25 << UInt(1) + D_26) << UInt(1) + D_27) << UInt(1) + D_28) << UInt(1) +
            D_29
        ) << UInt(1) + D_30
    computed_parity_bits == get_bits(word, 30, 25, 6)
end

function get_word(state::GNSSDecoderState{<:GPSL1Data}, word_number::Int)
    num_words = Int(state.constants.syncro_sequence_length / state.constants.word_length)
    word =
        state.buffer >> UInt(
            state.constants.word_length * (num_words - word_number) +
            state.constants.preamble_length,
        )
    UInt(word & (UInt(1) << UInt(state.constants.word_length) - UInt(1)))
end

function can_decode_word(decode_bits::Function, state::GNSSDecoderState, word_number::Int)
    word = get_word(state, word_number)
    prev_word = get_word(state, word_number - 1)
    prev_word_bit_30 = get_bit(prev_word, 30, 30)
    is_checked_word = word_number != 1 && word_number != 3
    if check_gpsl1_parity(
        word,
        is_checked_word * get_bit(prev_word, 30, 29),
        is_checked_word * prev_word_bit_30,
    )
        word_comp = (is_checked_word * prev_word_bit_30) ? ~word : word
        data = decode_bits(word_comp, state)
        state = GNSSDecoderState(state; raw_data = data)
    end
    return state
end

function can_decode_two_words(
    decode_bits::Function,
    state::GNSSDecoderState,
    word_number1::Int,
    word_number2::Int,
)
    word1 = get_word(state, word_number1)
    prev_word1 = get_word(state, word_number1 - 1)
    prev_word1_bit_30 = get_bit(prev_word1, 30, 30)
    is_checked_word1 = word_number1 != 1 && word_number1 != 3
    word2 = get_word(state, word_number2)
    prev_word2 = get_word(state, word_number2 - 1)
    prev_word2_bit_30 = get_bit(prev_word2, 30, 30)
    is_checked_word2 = word_number2 != 1 && word_number2 != 3
    if check_gpsl1_parity(
        word1,
        is_checked_word1 * get_bit(prev_word1, 30, 29),
        is_checked_word1 * prev_word1_bit_30,
    ) && check_gpsl1_parity(
        word2,
        is_checked_word2 * get_bit(prev_word2, 30, 29),
        is_checked_word2 * prev_word2_bit_30,
    )
        word1_comp = (is_checked_word1 * prev_word1_bit_30) ? ~word1 : word1
        word2_comp = (is_checked_word2 * prev_word2_bit_30) ? ~word2 : word2
        data = decode_bits(word1_comp, word2_comp, state)
        state = GNSSDecoderState(state; raw_data = data)
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
        TOW = get_bits(how_word, 30, 1, 17) * 6
        alert_flag = get_bit(how_word, 30, 18)
        anti_spoof_flag = get_bit(how_word, 30, 19)
        last_subframe_id = get_bits(how_word, 30, 20, 3)
        GPSL1Data(state.raw_data; last_subframe_id, TOW, alert_flag, anti_spoof_flag)
    end
    if !isnothing(prev_TOW) && prev_TOW + 1 != state.raw_data.TOW
        # Time of week must be decodable
        state = GNSSDecoderState(state; raw_data = GPSL1Data(state.raw_data; TOW = nothing))
    end
    state
end

function decode_syncro_sequence(state::GNSSDecoderState{<:GPSL1Data})
    state = read_tlm_and_how_words(state)
    subframe_id = state.raw_data.last_subframe_id

    if subframe_id == 1
        state = can_decode_word(state, 3) do word3, state
            trans_week = get_bits(word3, 30, 1, 10)

            # Codes on L2 Channel
            codeonl2 = get_bits(word3, 30, 11, 2)

            # SV Accuracy, user range accuracy
            ura = get_bits(word3, 30, 13, 4)
            if ura <= 6
                ura = 2^(1 + (ura / 2))
            elseif 6 < ura <= 14
                ura = 2^(ura - 2)
            elseif ura == 15
                ura = nothing
            end

            # Satellite Health
            svhealth = bitstring(get_bits(word3, 30, 17, 6))[end-5:end]
            if get_bit(word3, 30, 17)
                @warn "Bad LNAV Data, SV-Health critical", svhealth
            end

            GPSL1Data(state.raw_data; trans_week, codeonl2, ura, svhealth)
        end

        state = can_decode_two_words(state, 3, 8) do word3, word8, state
            # Issue of Data Clock
            # 2 MSB in Word 2, LSB 8 in Word 8
            IODC =
                bitstring(get_bits(word3, 30, 23, 2))[end-1:end] *
                bitstring(get_bits(word8, 30, 1, 8))[end-7:end]
            GPSL1Data(state.raw_data; IODC)
        end

        state = can_decode_word(state, 4) do word4, state
            # True: LNAV Datastream on PCode commanded OFF
            l2pcode = get_bit(word4, 30, 1)
            GPSL1Data(state.raw_data; l2pcode)
        end

        state = can_decode_word(state, 7) do word7, state
            # group time differential
            T_GD = get_twos_complement_num(word7, 30, 17, 8) / 1 << 31
            GPSL1Data(state.raw_data; T_GD)
        end

        state = can_decode_word(state, 8) do word8, state
            # Clock data reference
            t_0c = get_bits(word8, 30, 9, 16) << 4
            GPSL1Data(state.raw_data; t_0c)
        end

        state = can_decode_word(state, 9) do word9, state
            # clock correction parameter a_f2
            a_f2 = get_twos_complement_num(word9, 30, 1, 8) / 1 << 55
            # clock correction parameter a_f1
            a_f1 = get_twos_complement_num(word9, 30, 9, 16) / 1 << 43
            GPSL1Data(state.raw_data; a_f2, a_f1)
        end

        state = can_decode_word(state, 10) do word10, state
            # Clock data reference
            a_f0 = get_twos_complement_num(word10, 30, 1, 22) / 1 << 31
            GPSL1Data(state.raw_data; a_f0)
        end
    elseif subframe_id == 2
        state = can_decode_word(state, 3) do word3, state
            # Issue of ephemeris data
            IODE_Sub_2 = bitstring(get_bits(word3, 30, 1, 8))[end-7:end]
            C_rs = get_twos_complement_num(word3, 30, 9, 16) / 1 << 5
            GPSL1Data(state.raw_data; IODE_Sub_2, C_rs)
        end

        state = can_decode_word(state, 4) do word4, state
            # Mean motion difference from computed value
            Δn = get_twos_complement_num(word4, 30, 1, 16) * state.constants.PI / 1 << 43
            GPSL1Data(state.raw_data; Δn)
        end

        state = can_decode_two_words(state, 4, 5) do word4, word5, state
            # Mean motion difference from computed value
            combined_word =
                UInt(get_bits(word4, 30, 17, 8) << 24 + get_bits(word5, 30, 1, 24))
            M_0 =
                get_twos_complement_num(combined_word, 32, 1, 32) * state.constants.PI /
                1 << 31
            GPSL1Data(state.raw_data; M_0)
        end

        state = can_decode_word(state, 6) do word6, state
            # Amplitude of the Cosine Harmonic Correction Term to the Argument Latitude
            C_uc = get_twos_complement_num(word6, 30, 1, 16) / 1 << 29
            GPSL1Data(state.raw_data; C_uc)
        end

        state = can_decode_two_words(state, 6, 7) do word6, word7, state
            # Eccentricity
            e = (get_bits(word6, 30, 17, 8) << 24 + get_bits(word7, 30, 1, 24)) / 1 << 33
            GPSL1Data(state.raw_data; e)
        end

        state = can_decode_word(state, 8) do word8, state
            # Amplitude of the Sine Harmonic Correction Term to the Argument of Latitude
            C_us = get_twos_complement_num(word8, 30, 1, 16) / 1 << 29
            GPSL1Data(state.raw_data; C_us)
        end

        state = can_decode_two_words(state, 8, 9) do word8, word9, state
            # Square Root of Semi-Major Axis
            sqrt_A =
                (get_bits(word8, 30, 17, 8) << 24 + get_bits(word9, 30, 1, 24)) / 1 << 19
            GPSL1Data(state.raw_data; sqrt_A)
        end

        state = can_decode_word(state, 10) do word10, state
            # Reference Time ephemeris
            t_0e = get_bits(word10, 30, 1, 16) << 4
            fit_interval = get_bit(word10, 30, 17)
            AODO = get_bits(word10, 30, 18, 5)
            GPSL1Data(state.raw_data; t_0e, fit_interval, AODO)
        end
    elseif subframe_id == 3
        state = can_decode_word(state, 3) do word3, state
            # Amplitude of the Cosine Harmonic Correction to Angle of Inclination
            C_ic = get_twos_complement_num(word3, 30, 1, 16) / 1 << 29
            GPSL1Data(state.raw_data; C_ic)
        end

        state = can_decode_two_words(state, 3, 4) do word3, word4, state
            # Longitude of Ascending Node of Orbit Plane at Weekly Epoch
            combined_word =
                UInt(get_bits(word3, 30, 17, 8) << 24 + get_bits(word4, 30, 1, 24))
            Ω_0 =
                get_twos_complement_num(combined_word, 32, 1, 32) * state.constants.PI /
                1 << 31
            GPSL1Data(state.raw_data; Ω_0)
        end

        state = can_decode_word(state, 5) do word5, state
            # Amplitude of the sine harmonic correction term to angle of Inclination
            C_is = get_twos_complement_num(word5, 30, 1, 16) / 1 << 29
            GPSL1Data(state.raw_data; C_is)
        end

        state = can_decode_two_words(state, 5, 6) do word5, word6, state
            # inclination Angle at reference time
            combined_word =
                UInt(get_bits(word5, 30, 17, 8) << 24 + get_bits(word6, 30, 1, 24))
            i_0 =
                get_twos_complement_num(combined_word, 32, 1, 32) * state.constants.PI /
                1 << 31
            GPSL1Data(state.raw_data; i_0)
        end

        state = can_decode_word(state, 7) do word7, state
            # Amplitude of the cosine harmonic correction term to orbit Radius
            C_rc = get_twos_complement_num(word7, 30, 1, 16) / 1 << 5
            GPSL1Data(state.raw_data; C_rc)
        end

        state = can_decode_two_words(state, 7, 8) do word7, word8, state
            # Argument of Perigee
            combined_word =
                UInt(get_bits(word7, 30, 17, 8) << 24 + get_bits(word8, 30, 1, 24))
            ω =
                get_twos_complement_num(combined_word, 32, 1, 32) * state.constants.PI /
                1 << 31
            GPSL1Data(state.raw_data; ω)
        end

        state = can_decode_word(state, 9) do word9, state
            # Amplitude of the cosine harmonic correction term to orbit Radius
            Ω_dot = get_twos_complement_num(word9, 30, 1, 24) * state.constants.PI / 1 << 43
            GPSL1Data(state.raw_data; Ω_dot)
        end

        state = can_decode_word(state, 10) do word10, state
            # Issue of Ephemeris Data
            IODE_Sub_3 = bitstring(get_bits(word10, 30, 1, 8))[end-7:end]
            # Rate of Inclination Angle
            i_dot =
                get_twos_complement_num(word10, 30, 9, 14) * state.constants.PI / 1 << 43
            GPSL1Data(state.raw_data; IODE_Sub_3, i_dot)
        end
    elseif subframe_id == 4
        # Get page ID (SV ID) from word 3 bits 3-8
        state = can_decode_word(state, 3) do word3, state
            sv_page_id = get_bits(word3, 30, 3, 6)
            GPSL1Data(state.raw_data; last_subframe_id = 4 + sv_page_id * 100) # encode page in subframe_id temporarily
        end
        sv_page_id = (state.raw_data.last_subframe_id - 4) ÷ 100
        state = GNSSDecoderState(state; raw_data = GPSL1Data(state.raw_data; last_subframe_id = 4))

        if sv_page_id == 56 # Page 18: Ionospheric and UTC data
            state = decode_subframe4_page18(state)
        elseif sv_page_id == 63 # Page 25: A-S flags and SV health
            state = decode_subframe4_page25(state)
        elseif sv_page_id in 25:32 # Pages 2-5, 7-10: Almanac for SV 25-32
            state = decode_almanac_page(state, sv_page_id)
        end
    elseif subframe_id == 5
        # Get SV ID from word 3 bits 3-8
        state = can_decode_word(state, 3) do word3, state
            sv_id = get_bits(word3, 30, 3, 6)
            GPSL1Data(state.raw_data; last_subframe_id = 5 + sv_id * 100) # encode SV ID temporarily
        end
        sv_id = (state.raw_data.last_subframe_id - 5) ÷ 100
        state = GNSSDecoderState(state; raw_data = GPSL1Data(state.raw_data; last_subframe_id = 5))

        if sv_id == 51 # Page 25: SV health and almanac reference
            state = decode_subframe5_page25(state)
        elseif sv_id in 1:24 # Pages 1-24: Almanac for SV 1-24
            state = decode_almanac_page(state, sv_id)
        end
    end

    return state
end

function decode_subframe4_page18(state::GNSSDecoderState{<:GPSL1Data})
    # Page 18 contains ionospheric parameters and UTC parameters
    # Word 3: bits 9-16 = α0, bits 17-24 = α1
    # Word 4: bits 1-8 = α2, bits 9-16 = α3, bits 17-24 = β0
    # Word 5: bits 1-8 = β1, bits 9-16 = β2, bits 17-24 = β3
    # Word 6: bits 1-24 = A1 (24 bits)
    # Word 7: bits 1-24 = A0 MSBs (24 bits)
    # Word 8: bits 1-8 = A0 LSBs (8 bits), bits 9-16 = tot, bits 17-24 = WNt
    # Word 9: bits 1-8 = ΔtLS, bits 9-16 = WNLSF, bits 17-24 = DN
    # Word 10: bits 1-8 = ΔtLSF

    state = can_decode_word(state, 3) do word3, state
        α_0 = get_twos_complement_num(word3, 30, 9, 8) / 1 << 30
        α_1 = get_twos_complement_num(word3, 30, 17, 8) / 1 << 27
        GPSL1Data(state.raw_data; α_0, α_1)
    end

    state = can_decode_word(state, 4) do word4, state
        α_2 = get_twos_complement_num(word4, 30, 1, 8) / 1 << 24
        α_3 = get_twos_complement_num(word4, 30, 9, 8) / 1 << 24
        β_0 = get_twos_complement_num(word4, 30, 17, 8) * (1 << 11)
        GPSL1Data(state.raw_data; α_2, α_3, β_0)
    end

    state = can_decode_word(state, 5) do word5, state
        β_1 = get_twos_complement_num(word5, 30, 1, 8) * (1 << 14)
        β_2 = get_twos_complement_num(word5, 30, 9, 8) * (1 << 16)
        β_3 = get_twos_complement_num(word5, 30, 17, 8) * (1 << 16)
        GPSL1Data(state.raw_data; β_1, β_2, β_3)
    end

    state = can_decode_word(state, 6) do word6, state
        A_1 = get_twos_complement_num(word6, 30, 1, 24) / 1 << 50
        GPSL1Data(state.raw_data; A_1)
    end

    state = can_decode_two_words(state, 7, 8) do word7, word8, state
        # A0 is 32 bits: 24 MSBs in word 7, 8 LSBs in word 8
        combined_word = UInt(get_bits(word7, 30, 1, 24) << 8 + get_bits(word8, 30, 1, 8))
        A_0 = get_twos_complement_num(combined_word, 32, 1, 32) / 1 << 30
        GPSL1Data(state.raw_data; A_0)
    end

    state = can_decode_word(state, 8) do word8, state
        t_ot = get_bits(word8, 30, 9, 8) << 12
        WN_t = get_bits(word8, 30, 17, 8)
        GPSL1Data(state.raw_data; t_ot, WN_t)
    end

    state = can_decode_word(state, 9) do word9, state
        Δt_LS = get_twos_complement_num(word9, 30, 1, 8)
        WN_LSF = get_bits(word9, 30, 9, 8)
        DN = get_bits(word9, 30, 17, 8)
        GPSL1Data(state.raw_data; Δt_LS, WN_LSF, DN)
    end

    state = can_decode_word(state, 10) do word10, state
        Δt_LSF = get_twos_complement_num(word10, 30, 1, 8)
        GPSL1Data(state.raw_data; Δt_LSF)
    end

    return state
end

function decode_subframe4_page25(state::GNSSDecoderState{<:GPSL1Data})
    # Page 25 contains A-S flags and SV configurations for 32 SVs
    # and SV health for SV 25-32
    # Word 3: bits 9-24 = 4 SVs config (4 bits each)
    # Words 4-7: 24 MSBs = 6 SVs config each (4 bits each)
    # Word 8: bits 1-16 = 4 SVs config, bits 19-24 = SV25 health (6 bits)
    # Word 9: bits 1-24 = SV26-29 health (6 bits each)
    # Word 10: bits 1-18 = SV30-32 health (6 bits each)

    # Decode SV configurations (32 x 4-bit values)
    sv_config = Vector{Int64}(undef, 32)

    state = can_decode_word(state, 3) do word3, state
        for i in 1:4
            sv_config[i] = get_bits(word3, 30, 9 + (i - 1) * 4, 4)
        end
        GPSL1Data(state.raw_data; sv_config)
    end

    state = can_decode_word(state, 4) do word4, state
        cfg = something(state.raw_data.sv_config, Vector{Int64}(undef, 32))
        for i in 1:6
            cfg[4 + i] = get_bits(word4, 30, 1 + (i - 1) * 4, 4)
        end
        GPSL1Data(state.raw_data; sv_config = cfg)
    end

    state = can_decode_word(state, 5) do word5, state
        cfg = something(state.raw_data.sv_config, Vector{Int64}(undef, 32))
        for i in 1:6
            cfg[10 + i] = get_bits(word5, 30, 1 + (i - 1) * 4, 4)
        end
        GPSL1Data(state.raw_data; sv_config = cfg)
    end

    state = can_decode_word(state, 6) do word6, state
        cfg = something(state.raw_data.sv_config, Vector{Int64}(undef, 32))
        for i in 1:6
            cfg[16 + i] = get_bits(word6, 30, 1 + (i - 1) * 4, 4)
        end
        GPSL1Data(state.raw_data; sv_config = cfg)
    end

    state = can_decode_word(state, 7) do word7, state
        cfg = something(state.raw_data.sv_config, Vector{Int64}(undef, 32))
        for i in 1:6
            cfg[22 + i] = get_bits(word7, 30, 1 + (i - 1) * 4, 4)
        end
        GPSL1Data(state.raw_data; sv_config = cfg)
    end

    state = can_decode_word(state, 8) do word8, state
        cfg = something(state.raw_data.sv_config, Vector{Int64}(undef, 32))
        for i in 1:4
            cfg[28 + i] = get_bits(word8, 30, 1 + (i - 1) * 4, 4)
        end
        # SV 25 health (6 bits) at bits 19-24
        sv_health_sf4_25 = Vector{String}(undef, 8)
        sv_health_sf4_25[1] = bitstring(get_bits(word8, 30, 19, 6))[end-5:end]
        GPSL1Data(state.raw_data; sv_config = cfg, sv_health_sf4_25)
    end

    state = can_decode_word(state, 9) do word9, state
        health = something(state.raw_data.sv_health_sf4_25, Vector{String}(undef, 8))
        for i in 1:4
            health[1 + i] = bitstring(get_bits(word9, 30, 1 + (i - 1) * 6, 6))[end-5:end]
        end
        GPSL1Data(state.raw_data; sv_health_sf4_25 = health)
    end

    state = can_decode_word(state, 10) do word10, state
        health = something(state.raw_data.sv_health_sf4_25, Vector{String}(undef, 8))
        for i in 1:3
            health[5 + i] = bitstring(get_bits(word10, 30, 1 + (i - 1) * 6, 6))[end-5:end]
        end
        GPSL1Data(state.raw_data; sv_health_sf4_25 = health)
    end

    return state
end

function decode_subframe5_page25(state::GNSSDecoderState{<:GPSL1Data})
    # Page 25 contains SV health for SV 1-24 and almanac reference time/week
    # Word 3: bits 9-16 = toa, bits 17-24 = WNa
    # Word 4: bits 1-24 = SV1-4 health (6 bits each)
    # Word 5-8: bits 1-24 = SV health (6 bits each, 4 SVs per word)
    # Word 9: bits 1-24 = SV21-24 health (6 bits each)

    state = can_decode_word(state, 3) do word3, state
        t_oa = get_bits(word3, 30, 9, 8) << 12
        WN_a = get_bits(word3, 30, 17, 8)
        GPSL1Data(state.raw_data; t_oa, WN_a)
    end

    sv_health_sf5_25 = Vector{String}(undef, 24)

    state = can_decode_word(state, 4) do word4, state
        for i in 1:4
            sv_health_sf5_25[i] = bitstring(get_bits(word4, 30, 1 + (i - 1) * 6, 6))[end-5:end]
        end
        GPSL1Data(state.raw_data; sv_health_sf5_25)
    end

    state = can_decode_word(state, 5) do word5, state
        health = something(state.raw_data.sv_health_sf5_25, Vector{String}(undef, 24))
        for i in 1:4
            health[4 + i] = bitstring(get_bits(word5, 30, 1 + (i - 1) * 6, 6))[end-5:end]
        end
        GPSL1Data(state.raw_data; sv_health_sf5_25 = health)
    end

    state = can_decode_word(state, 6) do word6, state
        health = something(state.raw_data.sv_health_sf5_25, Vector{String}(undef, 24))
        for i in 1:4
            health[8 + i] = bitstring(get_bits(word6, 30, 1 + (i - 1) * 6, 6))[end-5:end]
        end
        GPSL1Data(state.raw_data; sv_health_sf5_25 = health)
    end

    state = can_decode_word(state, 7) do word7, state
        health = something(state.raw_data.sv_health_sf5_25, Vector{String}(undef, 24))
        for i in 1:4
            health[12 + i] = bitstring(get_bits(word7, 30, 1 + (i - 1) * 6, 6))[end-5:end]
        end
        GPSL1Data(state.raw_data; sv_health_sf5_25 = health)
    end

    state = can_decode_word(state, 8) do word8, state
        health = something(state.raw_data.sv_health_sf5_25, Vector{String}(undef, 24))
        for i in 1:4
            health[16 + i] = bitstring(get_bits(word8, 30, 1 + (i - 1) * 6, 6))[end-5:end]
        end
        GPSL1Data(state.raw_data; sv_health_sf5_25 = health)
    end

    state = can_decode_word(state, 9) do word9, state
        health = something(state.raw_data.sv_health_sf5_25, Vector{String}(undef, 24))
        for i in 1:4
            health[20 + i] = bitstring(get_bits(word9, 30, 1 + (i - 1) * 6, 6))[end-5:end]
        end
        GPSL1Data(state.raw_data; sv_health_sf5_25 = health)
    end

    return state
end

function decode_almanac_page(state::GNSSDecoderState{<:GPSL1Data}, sv_id::Int)
    # Almanac pages (subframe 4 pages 2-5, 7-10 for SV 25-32;
    #                subframe 5 pages 1-24 for SV 1-24)
    # Word 3: bits 9-24 = e (16 bits)
    # Word 4: bits 1-8 = toa, bits 9-24 = δi (16 bits)
    # Word 5: bits 1-16 = Ω_dot, bits 17-24 = SV health (8 bits)
    # Word 6: bits 1-24 = √A (24 bits)
    # Word 7: bits 1-24 = Ω0 (24 bits)
    # Word 8: bits 1-24 = ω (24 bits)
    # Word 9: bits 1-24 = M0 (24 bits)
    # Word 10: bits 1-8 = af0 MSBs (8 bits), bits 9-19 = af1 (11 bits), bits 20-22 = af0 LSBs (3 bits)

    if sv_id == 0
        # Dummy SV, skip decoding
        return state
    end

    # Decode all almanac words and build the almanac entry
    # We need to decode all words and check parity for each
    word3 = get_word(state, 3)
    word4 = get_word(state, 4)
    word5 = get_word(state, 5)
    word6 = get_word(state, 6)
    word7 = get_word(state, 7)
    word8 = get_word(state, 8)
    word9 = get_word(state, 9)
    word10 = get_word(state, 10)

    # Check parity for all words
    prev_word2 = get_word(state, 2)
    prev_word3 = get_word(state, 3)
    prev_word4 = get_word(state, 4)
    prev_word5 = get_word(state, 5)
    prev_word6 = get_word(state, 6)
    prev_word7 = get_word(state, 7)
    prev_word8 = get_word(state, 8)
    prev_word9 = get_word(state, 9)

    parity_ok =
        check_gpsl1_parity(word3, get_bit(prev_word2, 30, 29), get_bit(prev_word2, 30, 30)) &&
        check_gpsl1_parity(word4, get_bit(prev_word3, 30, 29), get_bit(prev_word3, 30, 30)) &&
        check_gpsl1_parity(word5, get_bit(prev_word4, 30, 29), get_bit(prev_word4, 30, 30)) &&
        check_gpsl1_parity(word6, get_bit(prev_word5, 30, 29), get_bit(prev_word5, 30, 30)) &&
        check_gpsl1_parity(word7, get_bit(prev_word6, 30, 29), get_bit(prev_word6, 30, 30)) &&
        check_gpsl1_parity(word8, get_bit(prev_word7, 30, 29), get_bit(prev_word7, 30, 30)) &&
        check_gpsl1_parity(word9, get_bit(prev_word8, 30, 29), get_bit(prev_word8, 30, 30)) &&
        check_gpsl1_parity(word10, get_bit(prev_word9, 30, 29), get_bit(prev_word9, 30, 30))

    if !parity_ok
        return state
    end

    # Apply bit complement if needed based on prev word bit 30
    word3_comp = get_bit(prev_word2, 30, 30) ? ~word3 : word3
    word4_comp = get_bit(prev_word3, 30, 30) ? ~word4 : word4
    word5_comp = get_bit(prev_word4, 30, 30) ? ~word5 : word5
    word6_comp = get_bit(prev_word5, 30, 30) ? ~word6 : word6
    word7_comp = get_bit(prev_word6, 30, 30) ? ~word7 : word7
    word8_comp = get_bit(prev_word7, 30, 30) ? ~word8 : word8
    word9_comp = get_bit(prev_word8, 30, 30) ? ~word9 : word9
    word10_comp = get_bit(prev_word9, 30, 30) ? ~word10 : word10

    # Decode almanac parameters
    alm_e = get_bits(word3_comp, 30, 9, 16) / 1 << 21
    alm_toa = get_bits(word4_comp, 30, 1, 8) << 12
    alm_δi = get_twos_complement_num(word4_comp, 30, 9, 16) * state.constants.PI / 1 << 19
    alm_Ω_dot = get_twos_complement_num(word5_comp, 30, 1, 16) * state.constants.PI / 1 << 38
    alm_sv_health = bitstring(get_bits(word5_comp, 30, 17, 8))[end-7:end]
    alm_sqrt_A = get_bits(word6_comp, 30, 1, 24) / 1 << 11
    alm_Ω_0 = get_twos_complement_num(word7_comp, 30, 1, 24) * state.constants.PI / 1 << 23
    alm_ω = get_twos_complement_num(word8_comp, 30, 1, 24) * state.constants.PI / 1 << 23
    alm_M_0 = get_twos_complement_num(word9_comp, 30, 1, 24) * state.constants.PI / 1 << 23

    # af0: 8 MSBs at bits 1-8, 3 LSBs at bits 20-22 = 11 bits total
    af0_msbs = get_bits(word10_comp, 30, 1, 8)
    af0_lsbs = get_bits(word10_comp, 30, 20, 3)
    af0_combined = UInt(af0_msbs << 3 + af0_lsbs)
    alm_af0 = get_twos_complement_num(af0_combined, 11, 1, 11) / 1 << 20
    alm_af1 = get_twos_complement_num(word10_comp, 30, 9, 11) / 1 << 38

    # Store almanac data for this SV
    almanac_entry = (
        e = alm_e,
        t_oa = alm_toa,
        δi = alm_δi,
        Ω_dot = alm_Ω_dot,
        sv_health = alm_sv_health,
        sqrt_A = alm_sqrt_A,
        Ω_0 = alm_Ω_0,
        ω = alm_ω,
        M_0 = alm_M_0,
        af0 = alm_af0,
        af1 = alm_af1,
    )

    almanac = something(state.raw_data.almanac, Dict{Int64,NamedTuple}())
    almanac[sv_id] = almanac_entry
    state = GNSSDecoderState(state; raw_data = GPSL1Data(state.raw_data; almanac))

    return state
end

function compare_data(data::GPSL1Data, new_data::GPSL1Data)
    data.IODC == new_data.IODC && # IODE_Sub_2 and IODE_Sub_3 is already checked for in validate_data
        (
            data.TOW == new_data.TOW ||
            (data.TOW > new_data.TOW && data.trans_week < new_data.trans_week) ||
            (data.TOW < new_data.TOW && data.trans_week >= new_data.trans_week)
        ) &&
        data.T_GD == new_data.T_GD &&
        data.t_0c == new_data.t_0c &&
        data.a_f0 == new_data.a_f0 &&
        data.a_f1 == new_data.a_f1 &&
        data.a_f2 == new_data.a_f2 &&
        data.C_rs == new_data.C_rs &&
        data.Δn == new_data.Δn &&
        data.M_0 == new_data.M_0 &&
        data.C_uc == new_data.C_uc &&
        data.e == new_data.e &&
        data.C_us == new_data.C_us &&
        data.sqrt_A == new_data.sqrt_A &&
        data.t_0e == new_data.t_0e &&
        data.fit_interval == new_data.fit_interval &&
        data.AODO == new_data.AODO &&
        data.C_ic == new_data.C_ic &&
        data.Ω_0 == new_data.Ω_0 &&
        data.C_is == new_data.C_is &&
        data.i_0 == new_data.i_0 &&
        data.C_rc == new_data.C_rc &&
        data.ω == new_data.ω &&
        data.Ω_dot == new_data.Ω_dot &&
        data.i_dot == new_data.i_dot
end

function increment_voting(old_vote, max_vote)
    min(max_vote, old_vote + 1)
end

function update_vote_at(old_data, idx, new_vote, new_data)
    [i == idx ? VotedGPSL1Data(new_vote, new_data) : entry for (i, entry) in enumerate(old_data)]
end

function confirm_data(state, max_vote = 20)
    old_data = state.cache.old_data

    # Check if any entry has same IODC
    has_same_iodc = any(e.data.IODC == state.raw_data.IODC for e in old_data)

    # Find matching entry (same IODC and same data)
    matching_idx = findfirst(old_data) do entry
        entry.data.IODC == state.raw_data.IODC && compare_data(entry.data, state.raw_data)
    end

    if isnothing(matching_idx)
        if has_same_iodc
            # Same IODC exists but data differs - add as new entry, don't use data yet
            return GNSSDecoderState(
                state;
                cache = GPSL1Cache([old_data; VotedGPSL1Data(0, state.raw_data)]),
                raw_data = GPSL1Data(),
            )
        else
            # New IODC entirely
            return if state.data == GPSL1Data() # no data yet - add to cache and use data
                GNSSDecoderState(
                    state;
                    cache = GPSL1Cache([VotedGPSL1Data(0, state.raw_data)]),
                    data = state.raw_data,
                    num_bits_after_valid_syncro_sequence = state.constants.preamble_length,
                )
            else # add as new entry, don't use data yet
                GNSSDecoderState(
                    state;
                    cache = GPSL1Cache([old_data; VotedGPSL1Data(0, state.raw_data)]),
                    raw_data = GPSL1Data(),
                )
            end
        end
    end

    # Found matching entry - upvote it
    curr_score = old_data[matching_idx].vote
    new_vote = increment_voting(curr_score, max_vote)

    # Find best score among entries with same IODC
    best_score = maximum(e.vote for e in old_data if e.data.IODC == state.raw_data.IODC)

    if best_score > curr_score
        # Another entry has higher score - reject this data
        new_cache = update_vote_at(old_data, matching_idx, new_vote, state.raw_data)
        return GNSSDecoderState(state; cache = GPSL1Cache(new_cache), raw_data = GPSL1Data())
    end

    # This entry has the best (or tied best) score - use the data
    new_cache = if new_vote == max_vote && length(old_data) > 1
        # Max votes reached - keep only this entry
        [VotedGPSL1Data(new_vote, state.raw_data)]
    else
        update_vote_at(old_data, matching_idx, new_vote, state.raw_data)
    end

    GNSSDecoderState(
        state;
        cache = GPSL1Cache(new_cache),
        data = state.raw_data,
        num_bits_after_valid_syncro_sequence = state.constants.preamble_length,
    )
end

function validate_data(state::GNSSDecoderState{<:GPSL1Data})
    if is_decoding_completed_for_positioning(state.raw_data) &&
       state.raw_data.IODC[3:10] == state.raw_data.IODE_Sub_2 == state.raw_data.IODE_Sub_3
        state = confirm_data(state)
    end
    return state
end

function is_sat_healthy(state::GNSSDecoderState{<:GPSL1Data})
    state.data.svhealth == "000000"
end