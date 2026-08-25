# CRC-24Q — the Q polynomial used by Galileo I/NAV pages, GPS L1C-D
# subframes 2 and 3, RTCM messages, and several other CNAV variants.
#
# Polynomial: g(x) = x^24 + x^23 + x^18 + x^17 + x^14 + x^11 + x^10 + x^7
#                  + x^6  + x^5  + x^4  + x^3  + x  + 1
# Wire form: 0x1864cfb, with the leading x^24 implicit -> 0x864cfb.
# Init: 0, no reflection, xor-out 0. Reference: RTCM 10403.3 §3.1.5,
# IS-GPS-800G §3.5.3.5, IS-GPS-200N §3.5.3.5.
#
# The CRC field appears at the end of the message, so a correctly received
# message satisfies `crc24q(whole_message_with_crc) == 0`. This replaced the
# legacy inline `galCRC24` (built from the `CRC` package) when Galileo E1B
# migrated to the soft-symbol decode path in issue #37; the `CRC` dependency
# is no longer required.

const CRC24Q_POLY = UInt32(0x01864cfb)
const CRC24Q_MASK = UInt32(0x00ffffff)

"""
    crc24q(bytes::AbstractVector{UInt8}) -> UInt32

Compute the CRC-24Q checksum (polynomial `0x1864cfb`, init `0`, no input or
output reflection, xor-out `0`) over `bytes`. The result is right-aligned in
the low 24 bits of the returned `UInt32`; bits 24..31 are always zero.

For a complete CRC-protected message — i.e. a message followed by its
big-endian 24-bit checksum — `crc24q(message_with_crc)` returns `0` iff
the checksum matches.
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

Bit-stream variant. `bits` is interpreted MSB-first in the same direction as
the wire (i.e. the first bit of `bits` enters the CRC register first).
`length(bits)` need not be a multiple of 8 — any leftover bits at the tail
are processed bit-by-bit. The CRC field, if appended, must therefore also
appear MSB-first as 24 individual bits.
"""
function crc24q(bits::AbstractVector{Bool})
    # Standard non-reflected bit-serial CRC. `CRC24Q_POLY` carries the
    # implicit x^24; we mask it away here because the register only holds
    # the 24 low bits.
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

Packed-word variant: compute the CRC over the low `num_bytes` octets of `word`,
most significant octet first. Equivalent to `crc24q(bytes)` for the same octets,
without materialising them.

This is the form the Galileo decoders need. A CRC-protected page arrives packed
into one wide `BitIntegers` word (`UInt288` for an I/NAV page pair, `UInt512` for
an E6-B C/NAV page), and its bit length is not a whole number of octets — 106 and
486 respectively. Rounding `num_bytes` up is safe and is what both callers do:
the register initialises to zero, so the few leading zero bits of the rounded-up
representation feed in as no-ops and leave the checksum unchanged.

Prefer this over `crc24q(reverse(digits(UInt8, word; base = 256)))`, which is
correct but allocates two vectors and, more importantly, drives `num_bytes`
software divisions on a wide integer — measured at 133 µs for a 61-octet
`UInt512` against 0.33 µs here, because `BitIntegers` division falls back to a
generic implementation while shifts and masks compile down.
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
