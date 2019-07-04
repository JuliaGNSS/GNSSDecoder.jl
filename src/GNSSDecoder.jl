module GNSSDecoder

    using DocStringExtensions, Parameters, FixedPointNumbers

    export init_decode,
    satPosition

    abstract type GNSSData end

    @with_kw mutable struct GPSData <: GNSSData
        IODC::Union{Nothing, String} = nothing
        IODE::Union{Nothing, String} = nothing
        C_rs::Union{Nothing, Float64} = nothing
        Δn::Union{Nothing, Float64} = nothing
        M0::Union{Nothing, Float64} = nothing
        C_uc::Union{Nothing, Float64} = nothing
        e::Union{Nothing, Float64} = nothing
        C_us::Union{Nothing, Float64} = nothing
        sqrt_A::Union{Nothing, Float64} = nothing
        t_oe::Union{Nothing, Float64} = nothing
        C_ic::Union{Nothing, Float64} = nothing
        𝛀_0::Union{Nothing, Float64} = nothing
        C_is::Union{Nothing, Float64} = nothing
        i_0::Union{Nothing, Float64} = nothing
        C_rc::Union{Nothing, Float64} = nothing
        ω::Union{Nothing, Float64} = nothing
        𝛀_dot::Union{Nothing, Float64} = nothing
        IDOT::Union{Nothing, Float64} = nothing
        ##
        tow::Union{Nothing, Float64} = nothing # all frames, word 2, bits 1-17
        alert::Union{Nothing, Bool} = nothing  # word 2, bit 18
        antiSpoof::Union{Nothing, Bool} = nothing  # word 2, bit 19
        transmissionWeekNumber::Union{Nothing, Float64} = nothing  # subframe 1, word 3, bit 1-10
        l2codes::Union{Nothing, Float64} = nothing  # subframe 1, word 3, bit 11-12
        ura::Union{Nothing, Float64} = nothing # subframe 1, word 3, bit 13-16
        svHealth::Union{Nothing, String} = nothing  # subframe 1, word 3, bit 17-22
        l2pcode::Union{Nothing, Bool} = nothing  # subframe 1, word 4, bit 1
        groupDelayDifferential::Union{Nothing, Float64} = nothing  # subframe 1, word 7, bit 17-24
        toc::Union{Nothing, Float64} = nothing # subframe 1, word 8, bit 9-24
        af2::Union{Nothing, Float64} = nothing # subframe 1, word 9, bit 1-8
        af1::Union{Nothing, Float64} = nothing # subframe 1, word 9, bit 9-24
        af0::Union{Nothing, Float64} = nothing # subframe 1, word 10, bit 1-22
        fitinterval::Union{Nothing, Bool} = nothing  # subframe 2, word 10, bit 17
        aodo::Union{Nothing, Float64} = nothing  # subframe 2, word 10, bit 18-22

    end

    @with_kw mutable struct preambles
        found_preamble::Bool = false
        found_inverted_preamble::Bool = false
        preamble_pos::Int = 0
    end


    @with_kw mutable struct parameters
        first_TLM::Bool = false
        #Intermediate variables for decoding
        subframe_count::Int = 0
        word_count::Int = 0
        prev_30::Bool = false
        prev_29::Bool = false
        word_window::Int = 1
        data_integrity::Bool = true
        new_data_requ::Bool = false #IODE must be equal to the 8 LSB of the IODC. If they don't match new data must be collected
        first_subframe::Bool = false
        second_subframe::Bool = false
        third_subframe::Bool = false
        fourth_subframe::Bool = false
        fifth_subframe::Bool = false
        decoding_completed::Bool = false #flag to indicate that a whole frame (1500 bits) has been decoded
    end

    @with_kw mutable struct GPSData_interm
        _M0::Union{Nothing, BitArray{1}} = nothing
        _e::Union{Nothing, BitArray{1}} = nothing
        _sqrt_A::Union{Nothing, BitArray{1}} = nothing
        _𝛀_0::Union{Nothing, BitArray{1}} = nothing
        _i_0::Union{Nothing, BitArray{1}} = nothing
        _ω::Union{Nothing, BitArray{1}} = nothing
    end

    @with_kw mutable struct nb
        num_bits_prev::Int = 0
    end

    function init_decode()
        buffer = BitArray(undef, 0)
        data = GPSData()
        found_preambles = preambles()
        parameters1 = parameters()
        n_b = nb()

        (data_bits, num_bits) ->_decode(buffer, data, data_bits, num_bits, found_preambles, parameters1, n_b)

    end

    function _decode(buffer, data, data_bits::UInt64, num_bits, found_preambles, parameters1, n_b)

        a = num_bits - n_b.num_bits_prev
        n_b.num_bits_prev = num_bits

        if a !== 0
        num_bits = a
        end

        for i = num_bits-1:-1:0
            current_bit = data_bits & (1 << i)
            pushfirst!(buffer, current_bit > 0)

            if !found_preambles.found_preamble && length(buffer) >= 308

                found_preambles.found_preamble,found_preambles.found_inverted_preamble, found_preambles.preamble_pos=find_preamble(buffer, found_preambles.preamble_pos)

                if found_preambles.found_preamble
                    parameters1.first_TLM = true
                end

            end


            if found_preambles.found_preamble && length(buffer) < 1500 + found_preambles.preamble_pos #if less than 5 subframes are found the signal cannot be decoded

                if parameters1.first_TLM
                    if found_preambles.found_inverted_preamble #if we found the inverted form of the preamble, the whole datastream needs to be inverted
                        buffer[1:end] = map(!,buffer[1:end])
                    end

                    if found_preambles.preamble_pos > 1
                        parameters1.prev_30 = buffer[length(buffer) - found_preambles.preamble_pos + 1]
                        parameters1.prev_29 = buffer[length(buffer) - found_preambles.preamble_pos + 2]
                    end

                    parameters1.first_TLM = false

                else
                    if found_preambles.found_inverted_preamble
                        buffer[1] = map(!,buffer[1])
                    end
                end
            end

            if found_preambles.found_preamble && length(buffer) == 1500 + found_preambles.preamble_pos
                println("DECODING...")
                if found_preambles.found_inverted_preamble
                    buffer[1] = map(!,buffer[1])
                end

                _interm = GPSData_interm() #Intermediate variables for GPSData
                buffer = buffer[1:(length(buffer)-found_preambles.preamble_pos)]

                for i=1:1:50
                    word,parameters1.word_count,parameters1.word_window=getword(buffer,parameters1.word_count,parameters1.word_window)
                    decodeword(word, data, parameters1,_interm)
                end

                parameters1.word_window = 1

                println("DECODING COMPLETED")

            end #end if found_preamble

            if (found_preambles.found_preamble) && (length(buffer) > 1500 + found_preambles.preamble_pos) && (length(buffer)-found_preambles.preamble_pos % 300 != 0)
                println("entras en la funcion 1")

                if found_preambles.found_inverted_preamble
                    buffer[1] = map(!,buffer[1])
                end

            end

            if (found_preambles.found_preamble) && (length(buffer) > 1500 + found_preambles.preamble_pos) && (length(buffer)-found_preambles.preamble_pos % 300 == 0)
                println("entras en la funcion 2")

                if found_preambles.found_inverted_preamble
                    buffer[1] = map(!,buffer[1])
                end

                _interm = GPSData_interm() #Intermediate variables for GPSData
                buffer = buffer[1:300]

                for i=1:1:10
                    word,parameters1.word_count,parameters1.word_window=getword(buffer,parameters1.word_count,parameters1.word_window)
                    decodeword(word, data, parameters1,_interm)
                end

            end

        end #end for

    end #end function _decode

    include("bin2dec.jl")
    include("getword.jl")
    include("paritycheck.jl")
    include("find_preamble.jl")
    include("decodeword.jl")
    include("satPosition.jl")

end #end module
