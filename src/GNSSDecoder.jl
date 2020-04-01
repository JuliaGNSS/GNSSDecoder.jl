module GNSSDecoder

    using DocStringExtensions, GNSSSignals, Parameters, FixedPointNumbers, StaticArrays, Logging
    using Unitful: Hz


    export decode,
        sat_position,
        GNSSDecoderState,
        correct_sat_time,
        get_sat_time,
        can_get_sat_position,
        pseudo_range

    @with_kw mutable struct GPSData
    IODC::Union{Nothing,String} = nothing
    IODE::Union{Nothing,String} = nothing
    C_rs::Union{Nothing,Float64} = nothing
    Δn::Union{Nothing,Float64} = nothing
    M0::Union{Nothing,Float64} = nothing
    C_uc::Union{Nothing,Float64} = nothing
    e::Union{Nothing,Float64} = nothing
    C_us::Union{Nothing,Float64} = nothing
    sqrt_A::Union{Nothing,Float64} = nothing
    t_oe::Union{Nothing,Float64} = nothing
    C_ic::Union{Nothing,Float64} = nothing
    𝛀_0::Union{Nothing,Float64} = nothing
    C_is::Union{Nothing,Float64} = nothing
    i_0::Union{Nothing,Float64} = nothing
    C_rc::Union{Nothing,Float64} = nothing
    ω::Union{Nothing,Float64} = nothing
    𝛀_dot::Union{Nothing,Float64} = nothing
    IDOT::Union{Nothing,Float64} = nothing
    tow::Union{Nothing,Float64} = nothing # all frames, word 2, bits 1-17
    alert::Union{Nothing,Bool} = nothing  # word 2, bit 18
    antiSpoof::Union{Nothing,Bool} = nothing  # word 2, bit 19
    transmissionWeekNumber::Union{Nothing,Float64} = nothing  # subframe 1, word 3, bit 1-10
    l2codes::Union{Nothing,Float64} = nothing  # subframe 1, word 3, bit 11-12
    ura::Union{Nothing,Float64} = nothing # subframe 1, word 3, bit 13-16
    svHealth::Union{Nothing,String} = nothing  # subframe 1, word 3, bit 17-22
    l2pcode::Union{Nothing,Bool} = nothing  # subframe 1, word 4, bit 1
    groupDelayDifferential::Union{Nothing,Float64} = nothing  # subframe 1, word 7, bit 17-24
    toc::Union{Nothing,Float64} = nothing # subframe 1, word 8, bit 9-24
    af2::Union{Nothing,Float64} = nothing # subframe 1, word 9, bit 1-8
    af1::Union{Nothing,Float64} = nothing # subframe 1, word 9, bit 9-24
    af0::Union{Nothing,Float64} = nothing # subframe 1, word 10, bit 1-22
    fitinterval::Union{Nothing,Bool} = nothing  # subframe 2, word 10, bit 17
    aodo::Union{Nothing,Float64} = nothing  # subframe 2, word 10, bit 18-22
    tow_computed::Union{Nothing, Int64} = nothing
    tot::Union{Nothing, Float64} = nothing 
 
end

    @with_kw mutable struct GNSSDecoderState
    buffer::BitArray{1} = falses(310)
    data::GPSData = GPSData()
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
    """
        stores element in a bitbuffer

        $(SIGNATURES)
        
        # Buffer must be of Type BitArray{N}
    """
    function push_buffer(buffer::BitArray{1}, val::Bool)
        buffer = buffer >> 1
        buffer.chunks[1] = buffer.chunks[1] + val
        return buffer
    end


    """
        Decodes the bits inserted using a bit buffer
        $(SIGNATURES)

        # Details
        checks, if the buffer has the expected size of 308 bits (length of subframe + preamble),
        then tries to find a preamble.
        
    """
    function decode(dc, data_bits::UInt64, num_bits)
    
        
        a = num_bits - dc.nb_prev
        dc.nb_prev = num_bits
        
        if a !== 0
            num_bits = a
        end

        
        for i = num_bits - 1:-1:0
            current_bit = data_bits & (1 << i)
            dc.buffer = push_buffer(dc.buffer, current_bit > 0)
            
            if dc.num_bits_buffered != 308
                dc.num_bits_buffered += 1
            end


            # * Find first preamble
            if (!dc.preamble_found) && (dc.num_bits_buffered == 308)
                dc.preamble_found = find_preamble(dc.buffer)
            end
        
            # * Begin decoding after preamble is confirmed
            if dc.preamble_found
                rev_buf = reverse(dc.buffer)
                words = map(i->rev_buf[(i * 30) + 3:((i + 1) * 30 ) + 2], 0:9) # gets the words of the subframe (word length = 30)
            dc.prev_29 = rev_buf[1]
                dc.prev_30 = rev_buf[2]
                dc = decode_words(dc, words)
                dc.preamble_found = false
            end

            

        end # end of for-loop 
    end # end of _decode()

    include("pseudo_range.jl")
    include("find_preamble.jl")
    include("decode_words.jl")
    include("parity_check.jl")
    include("bin2dec.jl")
    include("sat_position.jl")
end# end of module