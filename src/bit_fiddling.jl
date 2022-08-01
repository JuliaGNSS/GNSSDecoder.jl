function get_bit(word::Unsigned, word_length::Int, bit_number::Int)
    Bool(get_bits(word, word_length, bit_number, 1))
end

function get_bits(word::T, word_length::Int, start::Int, length::Int) where T <: Unsigned
    T((word >> (word_length - start - length + 1)) & (T(1) << length - T(1)))
end

function get_two_complement_num(word::Unsigned, word_length::Int, start::Int, length::Int)
    sign = get_bit(word, word_length, start)
    value = Int(get_bits(word, word_length, start, length))
    if sign
        return -1 << length + value
    else
        return value
    end
end

function deinterleave(bits::String, columns, rows)
    String(vec(permutedims(reshape(collect(bits), columns, rows))))
end

function invert_every_second_bit(bits::String)
    reshaped_bits = reshape(collect(bits), 2, length(bits) >> 1)
    reshaped_bits[2,:] .= ifelse.(reshaped_bits[2,:] .== '1', '0', '1')
    String(vec(reshaped_bits))
end