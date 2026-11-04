function get_bit(word::Unsigned, word_length::Int, bit_number::Int)
    Bool(get_bits(word, word_length, bit_number, 1))
end

function get_bits(word::T, word_length::Int, start::Int, length::Int) where {T<:Unsigned}
    T((word >> (word_length - start - length + 1)) & (T(1) << length - T(1)))
end

# Sign-extend through UInt64, not the platform UInt: several ICD fields are wider than 32
# bits (e.g. 33-bit BeiDou B-CNAV and GPS CNAV/CNAV-2 angles). `length` must be ≤ 64.
#
# The same 32-bit trap applies to scale factors: every scale of 2^-31 or finer is written
# `* 2.0^-n` rather than `/ 1 << n`, because on 32-bit Julia `1 << 31` is `typemin(Int32)`
# and `1 << 32` or wider is `0`. Scales up to 2^-30 are safe and left as shifts.
function get_twos_complement_num(word::Unsigned, word_length::Int, start::Int, length::Int)
    value = UInt64(get_bits(word, word_length, start, length))
    num_shift_bits = 64 - length
    signed(value << num_shift_bits) >> num_shift_bits
end
