function get_bit(word::Unsigned, word_length::Int, bit_number::Int)
    Bool(get_bits(word, word_length, bit_number, 1))
end

function get_bits(word::T, word_length::Int, start::Int, length::Int) where {T<:Unsigned}
    T((word >> (word_length - start - length + 1)) & (T(1) << length - T(1)))
end

# Sign-extend through UInt64 rather than the platform UInt: several ICD fields
# are wider than 32 bits (e.g. the 33-bit BeiDou B-CNAV and GPS CNAV/CNAV-2
# angle fields), which would overflow a 32-bit UInt. `length` must be ≤ 64.
#
# The same 32-bit trap applies to the *scale factors* a parser applies to the
# extracted field, which is why every scale of 2^-31 or finer is written
# `* 2.0^-n` rather than `/ 1 << n`: `1 << n` is an `Int`, so on 32-bit Julia
# `1 << 31` is `typemin(Int32)` (a sign flip) and `1 << 32` and wider are `0`
# (a division by zero). Anything up to 2^-30 is safe either way, and the
# narrow ones are left as shifts.
function get_twos_complement_num(word::Unsigned, word_length::Int, start::Int, length::Int)
    value = UInt64(get_bits(word, word_length, start, length))
    num_shift_bits = 64 - length
    signed(value << num_shift_bits) >> num_shift_bits
end
