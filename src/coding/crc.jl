# CRC-24Q, used by Galileo I/NAV pages, GPS L1C-D subframes 2 and 3, RTCM messages and
# several other CNAV variants.
#
# Polynomial: g(x) = x^24 + x^23 + x^18 + x^17 + x^14 + x^11 + x^10 + x^7
#                  + x^6  + x^5  + x^4  + x^3  + x  + 1
# Wire form: 0x1864cfb, with the leading x^24 implicit -> 0x864cfb.
# Init: 0, no reflection, xor-out 0. Reference: RTCM 10403.3 §3.1.5,
# IS-GPS-800G §3.5.3.5, IS-GPS-200N §3.5.3.5.
#
# The CRC field ends the message, so `crc24q(whole_message_with_crc) == 0` when received
# correctly.

const CRC24Q_POLY = UInt32(0x01864cfb)
const CRC24Q_MASK = UInt32(0x00ffffff)

"""
    crc24q(bytes::AbstractVector{UInt8}) -> UInt32

CRC-24Q checksum (polynomial `0x1864cfb`, init `0`, no reflection, xor-out `0`) over
`bytes`, in the low 24 bits of the returned `UInt32`.

For a message followed by its big-endian 24-bit checksum, `crc24q(message_with_crc)`
is `0` iff the checksum matches.
"""
function crc24q(bytes::AbstractVector{UInt8})
    crc = UInt32(0)
    for b in bytes
        crc ⊻= UInt32(b) << 16
        for _ = 1:8
            crc <<= 1
            if (crc & 0x01000000) != 0
                crc ⊻= CRC24Q_POLY
            end
        end
    end
    return crc & CRC24Q_MASK
end

"""
    crc24q(bits::AbstractVector{Bool}) -> UInt32

Bit-stream variant. The first bit of `bits` enters the CRC register first (wire order).
`length(bits)` need not be a multiple of 8. An appended CRC field must be 24 bits,
MSB-first.
"""
function crc24q(bits::AbstractVector{Bool})
    # `CRC24Q_POLY` carries the implicit x^24; mask it off, the register holds 24 bits.
    poly_low = CRC24Q_POLY & CRC24Q_MASK
    crc = UInt32(0)
    @inbounds for b in bits
        feedback = ((crc >> 23) & UInt32(0x01)) ⊻ UInt32(b)
        crc = (crc << 1) & CRC24Q_MASK
        if feedback != 0
            crc ⊻= poly_low
        end
    end
    return crc & CRC24Q_MASK
end

"""
    crc24q(word::Unsigned, num_bytes::Integer) -> UInt32

Packed-word variant: CRC over the low `num_bytes` octets of `word`, most significant
octet first. Equals `crc24q(bytes)` for the same octets, without materialising them.

Used by the Galileo decoders, where a CRC-protected page is packed into one wide
`BitIntegers` word (`UInt288` for an I/NAV page pair, `UInt512` for an E6-B C/NAV page)
whose bit length is not a whole number of octets (220 and 486). Rounding `num_bytes` up is
safe: the register starts at zero, so the leading zero bits are no-ops.

Prefer this over `crc24q(reverse(digits(UInt8, word; base = 256)))`: that allocates and
does `num_bytes` software divisions on a wide integer (`BitIntegers` division is generic;
133 µs for a 61-octet `UInt512` against 0.33 µs here).
"""
function crc24q(word::Unsigned, num_bytes::Integer)
    crc = UInt32(0)
    @inbounds for i = (num_bytes-1):-1:0
        crc ⊻= UInt32(UInt8((word >> (8 * i)) & 0xff)) << 16
        for _ = 1:8
            crc <<= 1
            if (crc & 0x01000000) != 0
                crc ⊻= CRC24Q_POLY
            end
        end
    end
    return crc & CRC24Q_MASK
end

# CRC-32Q, the ISM CRC of the GPS ISM Packet (IS-GPS-200 §30.3.3.10.1.2.2.13 as amended
# by IRN-IS-200N-004), the polynomial of the RTCA DO-229 / DO-246 32-bit CRCs.
#
# Polynomial: g(X) = X^32 + X^31 + X^24 + X^22 + X^16 + X^14 + X^8 + X^7 + X^5 + X^3
#                  + X + 1
# Wire form: 0x1814141ab, with the leading X^32 implicit -> 0x814141ab.
# Init: 0, no reflection, xor-out 0 (the ICD's r(X) = m(X)·X^32 mod g(X)).

const CRC32Q_POLY = UInt32(0x814141ab)

"""
    crc32q(word::Unsigned, word_length::Int, start::Int, len::Int) -> UInt32

CRC-32Q remainder (polynomial `0x1814141ab`, init `0`, no reflection, xor-out `0`) of
the `len` bits starting at 1-based bit `start` of `word`, an MSB-first `word_length`-bit
field. The first bit enters the register first (wire order), so for a field followed by
its 32-bit checksum `crc32q(word, word_length, start, len + 32) == 0` iff they match.
"""
function crc32q(word::Unsigned, word_length::Int, start::Int, len::Int)
    crc = UInt32(0)
    for i = start:(start+len-1)
        feedback = (crc >> 31) ⊻ UInt32(get_bit(word, word_length, i))
        crc <<= 1
        if feedback != 0
            crc ⊻= CRC32Q_POLY
        end
    end
    return crc
end
