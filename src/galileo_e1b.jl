# UInt288 buffer for Galileo E1B
# which holds at least a complete Galileo E1B page
# plus 10 extra syncronization bits
BitIntegers.@define_integers 288

Base.@kwdef struct GalileoE1BConstants <: AbstractGNSSConstants
    syncro_sequence_length::Int = 250
    preamble::UInt16 = 0b0101100000
    preamble_length::Int = 10
    PI::Float64 = 3.1415926535898
    Ω_dot_e::Float64 = 7.2921151467e-5
    c::Float64 = 2.99792458e8
    μ::Float64 = 3.986004418e14
    F::Float64 = -4.442807309e-10
end

# Page is splitted in even and odd parts
# Cache even part and decode after odd part
# Page contains 120 bits
struct GalileoE1BCache <: AbstractGNSSCache
    even_page_part_bits::UInt128
end

GalileoE1BCache() = GalileoE1BCache(UInt128(0))

@enum SignalHealth begin
    signal_ok
    signal_out_of_service
    signal_will_be_out_of_service
    signal_component_currently_in_test
end

@enum DataValidityStatus begin
    navigation_data_valid
    working_without_guarantee
end

Base.@kwdef struct GalileoE1BData <: AbstractGNSSData
    WN::Union{Nothing,Int64} = nothing
    TOW::Union{Nothing,Int64} = nothing

    t_0e::Union{Nothing,Float64} = nothing
    M_0::Union{Nothing,Float64} = nothing
    e::Union{Nothing,Float64} = nothing
    sqrt_A::Union{Nothing,Float64} = nothing
    Ω_0::Union{Nothing,Float64} = nothing
    i_0::Union{Nothing,Float64} = nothing
    ω::Union{Nothing,Float64} = nothing
    i_dot::Union{Nothing,Float64} = nothing
    Ω_dot::Union{Nothing,Float64} = nothing
    Δn::Union{Nothing,Float64} = nothing
    C_uc::Union{Nothing,Float64} = nothing
    C_us::Union{Nothing,Float64} = nothing
    C_rc::Union{Nothing,Float64} = nothing
    C_rs::Union{Nothing,Float64} = nothing
    C_ic::Union{Nothing,Float64} = nothing
    C_is::Union{Nothing,Float64} = nothing

    t_0c::Union{Nothing,Float64} = nothing
    a_f0::Union{Nothing,Float64} = nothing
    a_f1::Union{Nothing,Float64} = nothing
    a_f2::Union{Nothing,Float64} = nothing

    IOD_nav1::Union{Nothing,UInt} = nothing
    IOD_nav2::Union{Nothing,UInt} = nothing
    IOD_nav3::Union{Nothing,UInt} = nothing
    IOD_nav4::Union{Nothing,UInt} = nothing
    num_pages_after_last_TOW::Int = 0
    num_bits_after_valid_syncro_sequence_after_last_TOW::Union{Nothing,Int} = nothing

    signal_health_e1b::Union{Nothing,SignalHealth} = nothing
    signal_health_e5b::Union{Nothing,SignalHealth} = nothing
    data_validity_status_e1b::Union{Nothing,DataValidityStatus} = nothing
    data_validity_status_e5b::Union{Nothing,DataValidityStatus} = nothing

    broadcast_group_delay_e1_e5a::Union{Nothing,Float64} = nothing
    broadcast_group_delay_e1_e5b::Union{Nothing,Float64} = nothing
end

function GalileoE1BData(
    data::GalileoE1BData;
    WN = data.WN,
    TOW = data.TOW,
    t_0e = data.t_0e,
    M_0 = data.M_0,
    e = data.e,
    sqrt_A = data.sqrt_A,
    Ω_0 = data.Ω_0,
    i_0 = data.i_0,
    ω = data.ω,
    i_dot = data.i_dot,
    Ω_dot = data.Ω_dot,
    Δn = data.Δn,
    C_uc = data.C_uc,
    C_us = data.C_us,
    C_rc = data.C_rc,
    C_rs = data.C_rs,
    C_ic = data.C_ic,
    C_is = data.C_is,
    t_0c = data.t_0c,
    a_f0 = data.a_f0,
    a_f1 = data.a_f1,
    a_f2 = data.a_f2,
    IOD_nav1 = data.IOD_nav1,
    IOD_nav2 = data.IOD_nav2,
    IOD_nav3 = data.IOD_nav3,
    IOD_nav4 = data.IOD_nav4,
    num_pages_after_last_TOW = data.num_pages_after_last_TOW,
    num_bits_after_valid_syncro_sequence_after_last_TOW = data.num_bits_after_valid_syncro_sequence_after_last_TOW,
    signal_health_e1b = data.signal_health_e1b,
    signal_health_e5b = data.signal_health_e5b,
    data_validity_status_e1b = data.data_validity_status_e1b,
    data_validity_status_e5b = data.data_validity_status_e5b,
    broadcast_group_delay_e1_e5a = data.broadcast_group_delay_e1_e5a,
    broadcast_group_delay_e1_e5b = data.broadcast_group_delay_e1_e5b,
)
    GalileoE1BData(
        WN,
        TOW,
        t_0e,
        M_0,
        e,
        sqrt_A,
        Ω_0,
        i_0,
        ω,
        i_dot,
        Ω_dot,
        Δn,
        C_uc,
        C_us,
        C_rc,
        C_rs,
        C_ic,
        C_is,
        t_0c,
        a_f0,
        a_f1,
        a_f2,
        IOD_nav1,
        IOD_nav2,
        IOD_nav3,
        IOD_nav4,
        num_pages_after_last_TOW,
        num_bits_after_valid_syncro_sequence_after_last_TOW,
        signal_health_e1b,
        signal_health_e5b,
        data_validity_status_e1b,
        data_validity_status_e5b,
        broadcast_group_delay_e1_e5a,
        broadcast_group_delay_e1_e5b,
    )
end

function is_ephemeris_decoded(data::GalileoE1BData)
    !isnothing(data.t_0e) &&
        !isnothing(data.M_0) &&
        !isnothing(data.e) &&
        !isnothing(data.sqrt_A) &&
        !isnothing(data.Ω_0) &&
        !isnothing(data.i_0) &&
        !isnothing(data.ω) &&
        !isnothing(data.i_dot) &&
        !isnothing(data.Ω_dot) &&
        !isnothing(data.Δn) &&
        !isnothing(data.C_uc) &&
        !isnothing(data.C_us) &&
        !isnothing(data.C_rc) &&
        !isnothing(data.C_rs) &&
        !isnothing(data.C_ic) &&
        !isnothing(data.C_is)
end

function is_clock_correction_decoded(data::GalileoE1BData)
    !isnothing(data.t_0c) &&
        !isnothing(data.a_f0) &&
        !isnothing(data.a_f1) &&
        !isnothing(data.a_f2)
end

function is_health_status_decoded(data::GalileoE1BData)
    !isnothing(data.signal_health_e1b) &&
        !isnothing(data.signal_health_e5b) &&
        !isnothing(data.data_validity_status_e1b) &&
        !isnothing(data.data_validity_status_e5b)
end

function is_decoding_completed_for_positioning(data::GalileoE1BData)
    !isnothing(data.TOW) &&
        !isnothing(data.WN) &&
        !isnothing(data.broadcast_group_delay_e1_e5a) &&
        !isnothing(data.broadcast_group_delay_e1_e5b) &&
        is_ephemeris_decoded(data) &&
        is_clock_correction_decoded(data) &&
        is_health_status_decoded(data)
end

function GalileoE1BDecoderState(prn)
    GNSSDecoderState(
        prn,
        UInt288(0),
        UInt288(0),
        GalileoE1BData(),
        GalileoE1BData(),
        GalileoE1BConstants(),
        GalileoE1BCache(),
        0,
        nothing,
        false,
    )
end

function GNSSDecoderState(system::GalileoE1B, prn)
    GNSSDecoderState(
        prn,
        UInt288(0),
        UInt288(0),
        GalileoE1BData(),
        GalileoE1BData(),
        GalileoE1BConstants(),
        GalileoE1BCache(),
        0,
        nothing,
        false,
    )
end

function reset_decoder_state(state::GNSSDecoderState{<:GalileoE1BData})
    # Reset bit buffers and TOW data field, while keeping the
    # remaining parameters in raw_data. This allows a GNSSReceiver
    # to use a satellite after a reacquisition without waiting for
    # the decoding of all data fields.
    # Note: WN is currently not reset as it is broadcast not as
    # frequently as the TOW and thus may increase the time until
    # the decoder is available again after an outage. This will
    # lead to erroneous decoder information for a few seconds after
    # reacquisition when a new week started during a signal outage.
    GNSSDecoderState(
        state;
        raw_buffer = UInt288(0),
        buffer = UInt288(0),
        raw_data = GalileoE1BData(
            state.raw_data;
            TOW = nothing,
            num_bits_after_valid_syncro_sequence_after_last_TOW = nothing,
        ),
        data = GalileoE1BData(),
        num_bits_buffered = 0,
        num_bits_after_valid_syncro_sequence = nothing
    )
end

function decode_syncro_sequence(state::GNSSDecoderState{<:GalileoE1BData})
    encoded_bits = bitstring(state.buffer >> state.constants.preamble_length)[sizeof(
        state.buffer,
    )*8-state.constants.syncro_sequence_length+state.constants.preamble_length+1:end]
    deinterleaved_encoded_bits = deinterleave(encoded_bits, 30, 8)
    inv_deinterleaved_encoded_bits = invert_every_second_bit(deinterleaved_encoded_bits)
    decoded_bits = viterbi_decode(7, [79, 109], inv_deinterleaved_encoded_bits)
    bits = parse(UInt128, decoded_bits; base = 2)
    is_even = !get_bit(bits, 114, 1)
    is_nominal_page = !get_bit(bits, 114, 2)
    state = GNSSDecoderState(
        state;
        raw_data = GalileoE1BData(
            state.raw_data;
            num_pages_after_last_TOW = state.raw_data.num_pages_after_last_TOW + 1,
        ),
    )
    if is_even
        state = GNSSDecoderState(
            state;
            cache = GalileoE1BCache(is_nominal_page ? bits : UInt128(0)),
        )
    elseif state.cache.even_page_part_bits != 0 && is_nominal_page
        data =
            get_bits(state.cache.even_page_part_bits, 114, 3, 112) << 16 +
            get_bits(bits, 114, 3, 16)
        bits_to_check_CRC =
            UInt288(state.cache.even_page_part_bits) << 106 + get_bits(bits, 114, 1, 106)
        if galCRC24(reverse(digits(UInt8, bits_to_check_CRC; base = 256))) == 0
            data_type = get_bits(data, 128, 1, 6)
            if data_type == 0
                if get_bits(data, 128, 7, 2) == 2 # '10'
                    WN = get_bits(data, 128, 97, 12)
                    TOW = get_bits(data, 128, 109, 20)
                    state = GNSSDecoderState(
                        state;
                        raw_data = GalileoE1BData(
                            state.raw_data;
                            WN,
                            TOW,
                            num_pages_after_last_TOW = 1,
                            num_bits_after_valid_syncro_sequence_after_last_TOW = state.num_bits_after_valid_syncro_sequence,
                        ),
                    )
                end
            elseif data_type == 1
                IOD_nav1 = get_bits(data, 128, 7, 10)
                t_0e = get_bits(data, 128, 17, 14) * 60
                M_0 =
                    get_twos_complement_num(data, 128, 31, 32) * state.constants.PI /
                    1 << 31
                e = get_bits(data, 128, 63, 32) / 1 << 33
                sqrt_A = get_bits(data, 128, 95, 32) / 1 << 19
                state = GNSSDecoderState(
                    state;
                    raw_data = GalileoE1BData(
                        state.raw_data;
                        IOD_nav1,
                        t_0e,
                        M_0,
                        e,
                        sqrt_A,
                    ),
                )
            elseif data_type == 2
                IOD_nav2 = get_bits(data, 128, 7, 10)
                Ω_0 =
                    get_twos_complement_num(data, 128, 17, 32) * state.constants.PI /
                    1 << 31
                i_0 =
                    get_twos_complement_num(data, 128, 49, 32) * state.constants.PI /
                    1 << 31
                ω =
                    get_twos_complement_num(data, 128, 81, 32) * state.constants.PI /
                    1 << 31
                i_dot =
                    get_twos_complement_num(data, 128, 113, 14) * state.constants.PI /
                    1 << 43
                state = GNSSDecoderState(
                    state;
                    raw_data = GalileoE1BData(state.raw_data; IOD_nav2, Ω_0, i_0, ω, i_dot),
                )
            elseif data_type == 3
                IOD_nav3 = get_bits(data, 128, 7, 10)
                Ω_dot =
                    get_twos_complement_num(data, 128, 17, 24) * state.constants.PI /
                    1 << 43
                Δn =
                    get_twos_complement_num(data, 128, 41, 16) * state.constants.PI /
                    1 << 43
                C_uc = get_twos_complement_num(data, 128, 57, 16) / 1 << 29
                C_us = get_twos_complement_num(data, 128, 73, 16) / 1 << 29
                C_rc = get_twos_complement_num(data, 128, 89, 16) / 1 << 5
                C_rs = get_twos_complement_num(data, 128, 105, 16) / 1 << 5
                state = GNSSDecoderState(
                    state;
                    raw_data = GalileoE1BData(
                        state.raw_data;
                        IOD_nav3,
                        Ω_dot,
                        Δn,
                        C_uc,
                        C_us,
                        C_rc,
                        C_rs,
                    ),
                )
            elseif data_type == 4
                IOD_nav4 = get_bits(data, 128, 7, 10)
                C_ic = get_twos_complement_num(data, 128, 23, 16) / 1 << 29
                C_is = get_twos_complement_num(data, 128, 39, 16) / 1 << 29
                t_0c = get_bits(data, 128, 55, 14) * 60
                a_f0 = get_twos_complement_num(data, 128, 69, 31) / 1 << 34
                a_f1 = get_twos_complement_num(data, 128, 100, 21) / 1 << 46
                a_f2 = get_twos_complement_num(data, 128, 121, 6) / 1 << 59
                state = GNSSDecoderState(
                    state;
                    raw_data = GalileoE1BData(
                        state.raw_data;
                        IOD_nav4,
                        C_ic,
                        C_is,
                        t_0c,
                        a_f0,
                        a_f1,
                        a_f2,
                    ),
                )
            elseif data_type == 5
                broadcast_group_delay_e1_e5a =
                    get_twos_complement_num(data, 128, 48, 10) / 1 << 32
                broadcast_group_delay_e1_e5b =
                    get_twos_complement_num(data, 128, 58, 10) / 1 << 32
                signal_health_e5b = SignalHealth(get_bits(data, 128, 68, 2))
                signal_health_e1b = SignalHealth(get_bits(data, 128, 70, 2))
                data_validity_status_e5b = DataValidityStatus(get_bit(data, 128, 72))
                data_validity_status_e1b = DataValidityStatus(get_bit(data, 128, 73))
                WN = get_bits(data, 128, 74, 12)
                TOW = get_bits(data, 128, 86, 20)
                state = GNSSDecoderState(
                    state;
                    raw_data = GalileoE1BData(
                        state.raw_data;
                        broadcast_group_delay_e1_e5a,
                        broadcast_group_delay_e1_e5b,
                        signal_health_e5b,
                        signal_health_e1b,
                        data_validity_status_e5b,
                        data_validity_status_e1b,
                        WN,
                        TOW,
                        num_pages_after_last_TOW = 1,
                        num_bits_after_valid_syncro_sequence_after_last_TOW = state.num_bits_after_valid_syncro_sequence,
                    ),
                )
            elseif data_type == 6
                TOW = get_bits(data, 128, 106, 20)
                state = GNSSDecoderState(
                    state;
                    raw_data = GalileoE1BData(
                        state.raw_data;
                        TOW,
                        num_pages_after_last_TOW = 1,
                        num_bits_after_valid_syncro_sequence_after_last_TOW = state.num_bits_after_valid_syncro_sequence,
                    ),
                )
            end
        end
    end
    return state
end

function validate_data(state::GNSSDecoderState{<:GalileoE1BData})
    if is_decoding_completed_for_positioning(state.raw_data) &&
       state.raw_data.IOD_nav1 ==
       state.raw_data.IOD_nav2 ==
       state.raw_data.IOD_nav3 ==
       state.raw_data.IOD_nav4
        num_bits_after_valid_syncro_sequence = 0
        if state.data.TOW == state.raw_data.TOW
            num_bits_after_valid_syncro_sequence = state.num_bits_after_valid_syncro_sequence
        elseif !isnothing(state.raw_data.num_bits_after_valid_syncro_sequence_after_last_TOW)
            num_bits_after_valid_syncro_sequence = state.num_bits_after_valid_syncro_sequence -
                (state.raw_data.num_bits_after_valid_syncro_sequence_after_last_TOW -
                2 * state.constants.syncro_sequence_length - state.constants.preamble_length)
        else # first succesful decoding
            num_bits_after_valid_syncro_sequence = state.constants.preamble_length +
                (
                    state.raw_data.num_pages_after_last_TOW + 1
                ) * 250
        end
        state = GNSSDecoderState(
            state;
            data = state.raw_data,
            num_bits_after_valid_syncro_sequence,
        )
    end
    return state
end

function is_sat_healthy(state::GNSSDecoderState{<:GalileoE1BData})
    state.data.signal_health_e1b == signal_ok &&
        state.data.data_validity_status_e1b == navigation_data_valid
end
