function get_bit(word::Unsigned, word_length::Int, bit_number::Int)
    Bool(get_bits(word, word_length, bit_number, 1))
end

function get_bits(word::Unsigned, word_length::Int, start::Int, length::Int)
    Int((word >> (word_length - start - length + 1)) & (UInt(1) << UInt(length) - UInt(1)))
end

function get_two_complement_num(word::Unsigned, word_length::Int, start::Int, length::Int)
    sign = get_bit(word, word_length, start)
    value = get_bits(word, word_length, start, length)
    if sign
        return -1 << length + value
    else
        return value
    end
end