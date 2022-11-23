abstract type AbstractGNSSConstants end
abstract type AbstractGNSSData end
abstract type AbstractGNSSCache end

Base.@kwdef struct GNSSDecoderState{
    D<:AbstractGNSSData,
    C<:AbstractGNSSConstants,
    CA<:AbstractGNSSCache,
    B<:Unsigned,
}
    prn::Int
    raw_buffer::B
    buffer::B
    raw_data::D
    data::D
    constants::C
    cache::CA
    num_bits_buffered::Int = 0
    num_bits_after_valid_syncro_sequence::Union{Nothing,Int} = 0
    is_shifted_by_180_degrees = false
end

function GNSSDecoderState(
    state::GNSSDecoderState;
    raw_buffer = state.raw_buffer,
    buffer = state.buffer,
    raw_data = state.raw_data,
    data = state.data,
    cache = state.cache,
    num_bits_buffered = state.num_bits_buffered,
    num_bits_after_valid_syncro_sequence = state.num_bits_after_valid_syncro_sequence,
    is_shifted_by_180_degrees = state.is_shifted_by_180_degrees,
)
    GNSSDecoderState(
        state.prn,
        raw_buffer,
        buffer,
        raw_data,
        data,
        state.constants,
        cache,
        num_bits_buffered,
        num_bits_after_valid_syncro_sequence,
        is_shifted_by_180_degrees,
    )
end

function push_bit(state::GNSSDecoderState, current_bit)
    GNSSDecoderState(
        state;
        raw_buffer = state.raw_buffer << true + (current_bit > 0),
        num_bits_buffered = min(state.num_bits_buffered + 1, sizeof(state.raw_buffer) * 8),
    )
end

function is_enough_buffered_bits_to_decode(state::GNSSDecoderState)
    state.num_bits_buffered >=
    state.constants.syncro_sequence_length + state.constants.preamble_length
end

calc_preamble_mask(state::GNSSDecoderState) =
    UInt(1) << UInt(state.constants.preamble_length) - UInt(1)

function find_preamble(state::GNSSDecoderState)
    state.raw_buffer & calc_preamble_mask(state) == state.constants.preamble &&
        (state.raw_buffer >> state.constants.syncro_sequence_length) &
        calc_preamble_mask(state) == state.constants.preamble ||
        state.raw_buffer & calc_preamble_mask(state) ==
        ~state.constants.preamble & calc_preamble_mask(state) &&
            (state.raw_buffer >> state.constants.syncro_sequence_length) &
            calc_preamble_mask(state) ==
            ~state.constants.preamble & calc_preamble_mask(state)
end

function complement_buffer_if_necessary(state::GNSSDecoderState)
    if state.raw_buffer & calc_preamble_mask(state) ==
       ~state.constants.preamble & calc_preamble_mask(state)
        return GNSSDecoderState(
            state;
            buffer = ~state.raw_buffer,
            is_shifted_by_180_degrees = true,
        )
    else
        return GNSSDecoderState(
            state;
            buffer = state.raw_buffer,
            is_shifted_by_180_degrees = false,
        )
    end
end

function decode(
    state::GNSSDecoderState,
    bits::T,
    num_bits::Int;
    decode_once::Bool = false,
) where {T<:Unsigned}
    num_bits > sizeof(bits) * 8 ||
        ArgumentError("Number of bits is too large to fit type of bits")
    for i = num_bits-1:-1:0
        current_bit = bits & (T(1) << i)
        state = push_bit(state, current_bit)
        if !isnothing(state.num_bits_after_valid_syncro_sequence)
            state = GNSSDecoderState(
                state;
                num_bits_after_valid_syncro_sequence = state.num_bits_after_valid_syncro_sequence +
                                                       1,
            )
        end

        if is_enough_buffered_bits_to_decode(state) && find_preamble(state)
            state = decode_syncro_sequence(complement_buffer_if_necessary(state))
            if !decode_once || !is_decoding_completed_for_positioning(state.data)
                state = validate_data(state)
            end
            state =
                GNSSDecoderState(state; num_bits_buffered = state.constants.preamble_length)
        end
    end
    return state
end