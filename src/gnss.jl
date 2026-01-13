abstract type AbstractGNSSConstants end
abstract type AbstractGNSSData end
abstract type AbstractGNSSCache end

"""
$(TYPEDEF)

Generic decoder state for GNSS signal decoding. This parametric struct holds all state
required for decoding navigation messages from GNSS satellites.

# Type Parameters
- `D<:AbstractGNSSData`: The data type holding decoded navigation message fields
- `C<:AbstractGNSSConstants`: Constants specific to the GNSS system (e.g., preamble, timing)
- `CA<:AbstractGNSSCache`: Cache for intermediate decoding state
- `B<:Unsigned`: Buffer type for storing raw bits (sized for the specific GNSS system)

# Fields
$(TYPEDFIELDS)

# See Also
- [`GPSL1DecoderState`](@ref): Constructor for GPS L1 C/A decoder state
- [`GalileoE1BDecoderState`](@ref): Constructor for Galileo E1B decoder state
- [`decode`](@ref): Main function to decode bits using this state
- [`reset_decoder_state`](@ref): Reset decoder state after signal loss
"""
Base.@kwdef struct GNSSDecoderState{
    D<:AbstractGNSSData,
    C<:AbstractGNSSConstants,
    CA<:AbstractGNSSCache,
    B<:Unsigned,
}
    "Pseudo-Random Noise code identifier for the satellite"
    prn::Int
    "Raw bit buffer before phase correction"
    raw_buffer::B
    "Bit buffer after phase correction"
    buffer::B
    "Partially decoded navigation data (not yet validated)"
    raw_data::D
    "Validated navigation data ready for use"
    data::D
    "System-specific constants (preamble, timing parameters)"
    constants::C
    "Cache for intermediate decoding state"
    cache::CA
    "Number of bits currently stored in buffer"
    num_bits_buffered::Int = 0
    "Number of bits received after last valid synchronization sequence, or `nothing` if not yet synchronized"
    num_bits_after_valid_syncro_sequence::Union{Nothing,Int} = 0
    "Whether the signal phase is inverted by 180 degrees"
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

"""
$(TYPEDSIGNATURES)

Decode GNSS navigation message bits and update the decoder state.

Processes incoming bits from a GNSS signal, detecting preambles and decoding
synchronization sequences to extract navigation data. The function handles
both normal and 180-degree phase-shifted signals automatically.

# Arguments
- `state::GNSSDecoderState`: Current decoder state
- `bits::T`: Unsigned integer containing the bits to decode (MSB first)
- `num_bits::Int`: Number of valid bits in `bits` to process

# Keywords
- `decode_once::Bool=false`: If `true`, stops decoding after all required
  positioning data has been decoded (subframes 1-3 for GPS, pages 1-5 for Galileo)

# Returns
- `GNSSDecoderState`: Updated decoder state with newly decoded data

# Example
```julia
state = GPSL1DecoderState(1)  # PRN 1
state = decode(state, UInt8(0b10110010), 8)
```

# See Also
- [`GNSSDecoderState`](@ref): The state structure being updated
- [`is_sat_healthy`](@ref): Check if satellite health is good after decoding
"""
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