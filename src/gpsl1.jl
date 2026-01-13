# UInt320 buffer for GPS L1 (more efficient than UInt312)
# which holds at least a complete GPS L1 subframe plus
# 8 extra syncronization bíts
BitIntegers.@define_integers 320

"""
    GPSL1Constants

WGS 84 constants and LNAV message structure parameters for GPS L1 C/A signal decoding.

The physical constants are defined in IS-GPS-200 (Interface Specification) and are used
for computing satellite positions and clock corrections from broadcast ephemeris data.

# Fields
- `syncro_sequence_length::Int`: Length of synchronization sequence in bits (300 bits = 10 words × 30 bits)
- `preamble::UInt8`: TLM word preamble pattern (10001011 binary, 0x8B)
- `preamble_length::Int`: Length of preamble in bits (8)
- `word_length::Int`: Length of each LNAV word in bits (30)
- `PI::Float64`: Mathematical constant π = 3.1415926535898 (IS-GPS-200 Table 20-IV)
- `Ω_dot_e::Float64`: WGS 84 Earth rotation rate = 7.2921151467×10⁻⁵ rad/s
- `c::Float64`: Speed of light = 2.99792458×10⁸ m/s
- `μ::Float64`: WGS 84 Earth gravitational parameter = 3.986005×10¹⁴ m³/s²
- `F::Float64`: Relativistic correction constant = -4.442807633×10⁻¹⁰ s/√m

# Reference
IS-GPS-200N, Section 20.3.3 and Table 20-IV
"""
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

"""
    GPSL1Data

Decoded GPS L1 C/A LNAV navigation message data.

Contains ephemeris, clock correction, and satellite health parameters decoded from
subframes 1, 2, and 3 of the GPS LNAV message. All parameters conform to IS-GPS-200N.

# Telemetry and Handover Word (TLM/HOW) Fields
- `last_subframe_id::Int`: ID of the last decoded subframe (1-5)
- `integrity_status_flag::Bool`: LNAV data integrity status (0=OK, 1=bad)
- `TOW::Int64`: Time of Week at message transmission (seconds, 0-604784)
- `alert_flag::Bool`: URA may be worse than indicated (0=OK, 1=alert)
- `anti_spoof_flag::Bool`: Anti-spoofing mode (0=off, 1=on)

# Subframe 1 - Clock Correction Parameters
- `trans_week::Int64`: GPS week number (modulo 1024)
- `codeonl2::Int64`: Code on L2 channel (0=invalid, 1=P-code, 2=C/A-code, 3=invalid)
- `ura::Float64`: User Range Accuracy (meters), derived from URA index
- `svhealth::String`: 6-bit satellite health status ("000000" = healthy)
- `IODC::String`: Issue of Data, Clock (10-bit binary string)
- `l2pcode::Bool`: L2 P-code data flag (1=LNAV OFF on P-code)
- `T_GD::Float64`: L1-L2 group delay correction (seconds)
- `t_0c::Int64`: Clock reference time (seconds)
- `a_f0::Float64`: Clock bias correction coefficient (seconds)
- `a_f1::Float64`: Clock drift correction coefficient (s/s)
- `a_f2::Float64`: Clock drift rate correction coefficient (s/s²)

# Subframe 2 - Ephemeris Parameters (Part 1)
- `IODE_Sub_2::String`: Issue of Data, Ephemeris from subframe 2 (8-bit binary string)
- `C_rs::Float64`: Sine harmonic correction to orbit radius (meters)
- `Δn::Float64`: Mean motion difference from computed value (semi-circles/s)
- `M_0::Float64`: Mean anomaly at reference time (semi-circles)
- `C_uc::Float64`: Cosine harmonic correction to argument of latitude (rad)
- `e::Float64`: Eccentricity (dimensionless, range 0-0.03)
- `C_us::Float64`: Sine harmonic correction to argument of latitude (rad)
- `sqrt_A::Float64`: Square root of semi-major axis (√m)
- `t_0e::Int64`: Ephemeris reference time (seconds)
- `fit_interval::Bool`: Curve fit interval flag (0=4h, 1=>4h)
- `AODO::Int64`: Age of Data Offset for NMCT (seconds)

# Subframe 3 - Ephemeris Parameters (Part 2)
- `C_ic::Float64`: Cosine harmonic correction to inclination (rad)
- `Ω_0::Float64`: Longitude of ascending node at weekly epoch (semi-circles)
- `C_is::Float64`: Sine harmonic correction to inclination (rad)
- `i_0::Float64`: Inclination angle at reference time (semi-circles)
- `C_rc::Float64`: Cosine harmonic correction to orbit radius (meters)
- `ω::Float64`: Argument of perigee (semi-circles)
- `Ω_dot::Float64`: Rate of right ascension (semi-circles/s)
- `IODE_Sub_3::String`: Issue of Data, Ephemeris from subframe 3 (8-bit binary string)
- `i_dot::Float64`: Rate of inclination angle (semi-circles/s)

# Reference
IS-GPS-200N, Tables 20-I, 20-II, 20-III, Sections 20.3.3.3-20.3.3.4
"""
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

"""
$(TYPEDSIGNATURES)

Create a decoder state for GPS L1 C/A navigation messages.

Initializes a [`GNSSDecoderState`](@ref) configured for decoding GPS L1 C/A
(Coarse/Acquisition) civil navigation messages. The decoder extracts ephemeris,
clock correction, and health data from the 50 bps LNAV data stream.

# Arguments
- `prn::Int`: Pseudo-Random Noise code identifier (1-32 for GPS satellites)

# Returns
- `GNSSDecoderState{GPSL1Data}`: Initialized decoder state for GPS L1

# Example
```julia
state = GPSL1DecoderState(1)  # Create decoder for PRN 1
state = decode(state, bits, num_bits)
if is_sat_healthy(state)
    # Use state.data for positioning
end
```

# See Also
- [`GNSSDecoderState`](@ref): The underlying state structure
- [`decode`](@ref): Decode bits using this state
- [`reset_decoder_state`](@ref): Reset after signal loss
- [`is_sat_healthy`](@ref): Check satellite health status
"""
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

"""
$(TYPEDSIGNATURES)

Reset the GPS L1 decoder state after a signal loss or reacquisition.

Clears the bit buffers and time-of-week (TOW) field while preserving other
decoded ephemeris and clock data in `raw_data`. This allows faster recovery
after brief signal outages without requiring a full re-decode of all subframes.

!!! note
    The `trans_week` field is intentionally not reset as it is only broadcast
    in subframe 1. This may cause brief errors if a GPS week rollover occurs
    during a signal outage.

# Arguments
- `state::GNSSDecoderState{<:GPSL1Data}`: Current GPS L1 decoder state

# Returns
- `GNSSDecoderState{<:GPSL1Data}`: Reset decoder state with cleared buffers

# Example
```julia
# After detecting signal loss
state = reset_decoder_state(state)
# Continue decoding with preserved ephemeris
state = decode(state, new_bits, num_bits)
```

# See Also
- [`GPSL1DecoderState`](@ref): Create a fresh decoder state
- [`decode`](@ref): Continue decoding after reset
"""
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
    end

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

"""
$(TYPEDSIGNATURES)

Check if the GPS satellite is healthy and usable for positioning.

Examines the 6-bit satellite health field (`svhealth`) from subframe 1. A satellite
is considered healthy only if all health bits are zero (`"000000"`).

!!! warning
    This function requires that subframe 1 has been successfully decoded.
    Check that `state.data.svhealth` is not `nothing` before relying on this result.

# Arguments
- `state::GNSSDecoderState{<:GPSL1Data}`: GPS L1 decoder state with decoded data

# Returns
- `Bool`: `true` if satellite health status indicates normal operation

# Example
```julia
state = GPSL1DecoderState(1)
state = decode(state, bits, num_bits)
if is_sat_healthy(state)
    # Safe to use for positioning
end
```

# See Also
- [`GPSL1DecoderState`](@ref): Create decoder state
- [`decode`](@ref): Decode navigation data
"""
function is_sat_healthy(state::GNSSDecoderState{<:GPSL1Data})
    state.data.svhealth == "000000"
end