function get_bit(word::Unsigned, word_length::Int, bit_number::Int)
    Bool(get_bits(word, word_length, bit_number, 1))
end

function get_bits(word::T, word_length::Int, start::Int, length::Int) where {T<:Unsigned}
    T((word >> (word_length - start - length + 1)) & (T(1) << length - T(1)))
end

function get_twos_complement_num(word::Unsigned, word_length::Int, start::Int, length::Int)
    value = UInt(get_bits(word, word_length, start, length))
    num_shift_bits = sizeof(value) * 8 - length
    signed(value << num_shift_bits) >> num_shift_bits
end