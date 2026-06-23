using Test
using GNSSDecoder
using GNSSDecoder: crc24q

@testset "CRC-24Q" begin
    @testset "Empty input" begin
        @test crc24q(UInt8[]) == 0
    end

    @testset "RTCM 10403.3 / IS-GPS-800 reference vector" begin
        # Standard CRC-24Q "check" vector: CRC over the ASCII byte sequence
        # "123456789" is 0xCDE703 (see CRC.jl's `spec(24, 0x864cfb, ...,
        # 0xcde703)` and Greg Cook's CRC catalogue under CRC-24/OPENPGP).
        msg = UInt8.(collect("123456789"))
        @test crc24q(msg) == 0x00CDE703
    end

    @testset "Single-byte sanity" begin
        # Single byte 0x00 should give CRC 0; single byte 0x01 should give
        # exactly the polynomial low-24 (8 shift steps with a 1 in MSB).
        @test crc24q(UInt8[0x00]) == 0
        # Compute manually: register = 0x010000, then 8 shift-and-condxor steps.
        crc = UInt32(0x010000)
        for _ = 1:8
            crc <<= 1
            if (crc & 0x01000000) != 0
                crc ⊻= UInt32(0x01864cfb)
            end
        end
        @test crc24q(UInt8[0x01]) == crc & UInt32(0x00ffffff)
    end

    @testset "Equivalence: byte API vs bit API" begin
        # The two APIs must produce identical CRCs when fed identical bit
        # streams. Convert "123456789" to bits MSB-first and check.
        msg = UInt8.(collect("123456789"))
        bits = Bool[]
        for byte in msg
            for bitpos = 7:-1:0
                push!(bits, (byte >> bitpos) & 0x1 == 0x1)
            end
        end
        @test crc24q(bits) == crc24q(msg)
    end

    @testset "Self-consistency: message + appended CRC checks to zero" begin
        # Appending the CRC to the message should make the whole thing CRC
        # to zero — that's the receiver's check.
        for raw in (
            UInt8.(collect("123456789")),
            UInt8[0x55, 0xaa, 0xff, 0x00, 0x12, 0x34],
            UInt8[i for i = 0x00:0x10],
        )
            c = crc24q(raw)
            tail = UInt8[(c>>16)&0xff, (c>>8)&0xff, c&0xff]
            @test crc24q(vcat(raw, tail)) == 0
        end
    end

    @testset "Frozen golden values (legacy galCRC24 parity)" begin
        # `crc24q` replaced the legacy `CRC.jl`-based `galCRC24` when Galileo
        # E1B migrated to the soft-symbol decode path (issue #37). These golden
        # values were captured from `galCRC24` before its removal, so they pin
        # bit-identical behaviour across the migration. The "123456789" entry
        # is the canonical CRC-24Q check value (0xCDE703).
        for (raw, expected) in (
            (UInt8.(collect("123456789")), UInt32(0x00cde703)),
            (UInt8[0x55, 0xaa, 0xff, 0x00, 0x12, 0x34], UInt32(0x005126f1)),
            (UInt8[i for i = 0x00:0x20], UInt32(0x00bf660c)),
        )
            @test crc24q(raw) == expected
        end
    end
end
