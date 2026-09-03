# Precompile workload (PrecompileTools). The first `decode` of a real navigation
# stream costs 12.5 s of compilation on a workstation and three to four times
# that on an embedded ARM host — paid by a live receiver at the moment its first
# satellite has bit sync, tens of seconds into a run, while every tracking loop
# waits (GNSSReceiver.jl#107). Every decoder this package has runs here, so a
# receiver on any signal finds it compiled: real captured streams where the
# test suite has them (GPS L1 C/A and Galileo E1B as hex — the latter also
# drives E5b, whose I/NAV pages are the same; GPS L2C and L5I CNAV messages and
# Galileo E5a F/NAV pages from the test suite's own `test/data/`, encoded to
# symbols the way the transmitter would), structurally valid pages for Galileo E6-B, and random
# symbols for the decoders without a capture (GPS L1C-D, BeiDou
# B1I/B1C/B2a/B2b/B3I), which
# still compiles their sync search, FEC and CRC paths, just not the field
# parsing that only a valid frame reaches.
using PrecompileTools: @setup_workload, @compile_workload

const _PRECOMPILE_GPS_L1CA_SYMBOLS =
    "8baaa8beae7de8c0e9000f75555540aaaaaaf2aaaaabcaaaa0257fafbd4593dc273ffb9e8f48baaa" *
    "8beae7bd1c0100c0660003861bf59b8e80864cfbcccccc5010ca84df264bc510429fa48baaa8beae" *
    "79c94064b410c9688d79f79cd40d52e77cdff5f32eba317a7ff29406c004370148baaa8beae77b64" *
    "79aaaaa5555556aaaaaaa55555556aaaaaaa55555556aaaaaaa55555559c8baaa8beae75aec50000" *
    "0a1f02d85ffffffffd7bc9904a35103943ee7972a0f43ab7fffff208baaa8beae73ee80e9000f755" *
    "55540aaaaaaf2aaaaabcaaaa0257fafbd4593dc273ffb9e8f48baaa8beae71d800100c0660003861" *
    "bf59b8e80864cfbcccccc5010ca84df264bc510429fa48baaa8beae6fc74064b410c9688d79f79cd" *
    "40d52e77cdff5f32eba317a7ff29406c004370148baaa8beae6db307753509edab6e87b1abdfd6f2" *
    "c2cb84d554e01250d045abb6b001394eab08baaa8beae6ba40510000960fd2789ffffff5a84366d2" *
    "cdc02e9dc02f43b8b4b1557fffff208baaa8beae69ebc0e9000f75555540aaaaaaf2aaaaabcaaaa0" *
    "257fafbd4593dc273ffb9e8f48baaa8beae67d600100c0660003861bf59b8e80864cfbcccccc5010" *
    "ca84df264bc510429fa48baaa8beae65ce8064b410c9688d79f79cd40d52e77cdff5f32eba317a7f" *
    "f29406c004370148baaa8beae63b54780502a400"

const _PRECOMPILE_GALILEO_E1B_SYMBOLS =
    "580ba78311d174abc085cca0ec39253e13c2ff0fae80ff9d381d939d486fb29609400f2e5cffc233" *
    "a00002ae1ffbd94e002fda08ff43ce8c020bd5cfffc5a580f10b21cc68ef3ad013f196a7c8f6692a" *
    "750cfca2db1a9885e2b89c03f6dd60fc00f1c5bffc283c400032e6ffbc504002f3b20ff4052ec020" *
    "4574fffcea5804243d38ef0fac918402d459eff55640602d1b9fff7b74870786bbe9e7729609400f" *
    "5f56ffc383a00000dfaffbd680402f583eff407608020294bfffc925808a00162e77ffef510001a0" *
    "ae7ffefcca002e6affff6e9ae00067f87ff27960cc00f1543ffc327d400026e2ffbce4f002fc606f" *
    "f43c6e40206752ffff1a5810010000bffbffcd000000adff3fff1807000f9fe3fffd0060015bffff" *
    "fb160e000f4347ffc387c000046ffffbd404002f872fff42ba00020d74ffffdea5813fbb0b5e04fd" *
    "c05e2d4f89b1959ee136ec032b3bd3af3ff92d31176a1fd960b400f4c51ffc22f8c0004de5ffbd0c" *
    "fc02f0d0cff40c2a8020935ffffe5258143b6667ebc63c1971394c6b98cbf40af9d331fbf5befac4" *
    "0037a0cf1d57d60dc00f2e55ffc2178c0006decffbc0ca802f9d0dff43c220020cb65fffde258114" *
    "5b160daf5a411a06a19992fc2f151cd5701b7c0af9a040aa5e72d1cb4560d000f1344ffc34be4000" *
    "6ff9ffbd281c02f923dff43863c020235efffe0a"

# Hex digits → ±1 soft symbols (a `0` bit is +1), most significant bit first.
function _precompile_soft_symbols(hex::AbstractString)
    out = Vector{Float32}(undef, 4length(hex))
    for (i, c) in enumerate(hex)
        nibble = parse(UInt8, c; base = 16)
        for b = 0:3
            out[4(i-1)+b+1] = (nibble >> (3 - b)) & 0x1 == 0 ? 1.0f0 : -1.0f0
        end
    end
    out
end

# 26 CNAV messages of 300 bits, stored as 38-byte blocks.
function _precompile_cnav_messages(path)
    raw = read(path)
    map(0:25) do i
        block = raw[(38i+1):(38i+38)]
        bits = [(byte >> (7 - j)) & 0x01 == 0x01 for byte in block for j = 0:7]
        bits[1:300]
    end
end

# The continuous K=7 rate-1/2 FEC (G1 = 0o171, G2 = 0o133) a CNAV transmitter
# applies, as ±1 symbols (bit 0 ⇒ +1).
function _precompile_cnav_symbols(messages)
    register = 0x00
    soft = Float32[]
    for bit in Iterators.flatten(messages)
        u = UInt8(bit)
        s1 = (register >> 5) & 0x01
        s2 = (register >> 4) & 0x01
        s3 = (register >> 3) & 0x01
        s5 = (register >> 1) & 0x01
        s6 = register & 0x01
        push!(soft, (u ⊻ s1 ⊻ s2 ⊻ s3 ⊻ s6) == 0x01 ? -1.0f0 : 1.0f0)
        push!(soft, (u ⊻ s2 ⊻ s3 ⊻ s5 ⊻ s6) == 0x01 ? -1.0f0 : 1.0f0)
        register = ((u << 5) | (register >> 1)) & 0x3f
    end
    soft
end

# One Galileo page on air: the sync pattern, then the information bits
# convolutionally encoded (K=7, G1 = 0o171, G2 = 0o133 inverted) and block
# interleaved over `GALILEO_INTERLEAVER_ROWS` rows. `columns` is the only shape
# that differs between the channels — 61 for E5a's F/NAV, 123 for E6-B's C/NAV.
function _precompile_galileo_page_symbols(info::Vector{Bool}, sync, columns::Int)
    register = zeros(Bool, 6)
    encoded = Float32[]
    for b in vcat(info, zeros(Bool, 6))
        taps = vcat([b], register)
        for (k, poly) in enumerate((0o171, 0o133))
            acc = false
            for j = 1:7
                ((poly >> (7 - j)) & 1 == 1) && (acc ⊻= taps[j])
            end
            push!(encoded, (k == 2 ? !acc : acc) ? -1.0f0 : 1.0f0)
        end
        register = vcat([b], register[1:(end-1)])
    end
    vcat(
        Float32[b == 1 ? -1.0f0 : 1.0f0 for b in sync],
        interleave(encoded, columns, GALILEO_INTERLEAVER_ROWS),
    )
end

# A stream of such pages, closed by a trailing sync so the last page's
# both-ends preamble match succeeds.
function _precompile_galileo_stream(pages::Vector{Vector{Bool}}, sync, columns::Int)
    stream = Float32[]
    for page in pages
        append!(stream, _precompile_galileo_page_symbols(page, sync, columns))
    end
    append!(stream, Float32[b == 1 ? -1.0f0 : 1.0f0 for b in sync])
    stream
end

# Galileo E5a F/NAV pages (244 bits each, stored as 31-byte blocks) → the 238
# information bits each carries.
function _precompile_e5a_pages(path)
    raw = read(path)
    map(1:(length(raw)÷31)) do p
        page = UInt256(0)
        for b = 1:31
            page = (page << 8) | UInt256(raw[(p-1)*31+b])
        end
        info = (page >> 4) >> 6
        Bool[get_bit(info, 238, k) for k = 1:238]
    end
end

# The captured fixtures are the ones the test suite already ships in
# `test/data/`, read from there rather than copied into `src/`: a package
# tarball carries both directories, and one copy cannot drift from the other.
# A tree without them still precompiles — `_precompile_fixture_symbols` falls
# back to random symbols, which compiles the affected decoder's sync search,
# FEC and CRC paths, just not its field parsing.
_precompile_fixture(name) = joinpath(@__DIR__, "..", "test", "data", name)

function _precompile_fixture_symbols(build, name, fallback)
    path = _precompile_fixture(name)
    isfile(path) ? build(path) : fallback
end

@setup_workload begin
    random_symbols = Float32[rand(Bool) ? 1.0f0 : -1.0f0 for _ = 1:8000]
    gps_l1ca = _precompile_soft_symbols(_PRECOMPILE_GPS_L1CA_SYMBOLS)
    galileo_e1b = _precompile_soft_symbols(_PRECOMPILE_GALILEO_E1B_SYMBOLS)
    gps_l2c =
        _precompile_fixture_symbols("gps_l2c_prn25_nav_bits.bin", random_symbols) do path
            _precompile_cnav_symbols(_precompile_cnav_messages(path))
        end
    gps_l5i =
        _precompile_fixture_symbols("gps_l5i_prn25_nav_bits.bin", random_symbols) do path
            _precompile_cnav_symbols(_precompile_cnav_messages(path))
        end
    galileo_e5a =
        _precompile_fixture_symbols("galileo_e5a_fnav_pages.bin", random_symbols) do path
            _precompile_galileo_stream(
                _precompile_e5a_pages(path),
                (1, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0),
                61,
            )
        end
    # E6-B C/NAV: 486 information bits per page. The HAS payload the pages carry
    # is not reproduced here — a structurally valid page (right sync, FEC and
    # interleaving) already compiles the sync search, deinterleaver, Viterbi and
    # CRC, which is the decoder's hot path.
    galileo_e6b = _precompile_galileo_stream(
        [Bool[isodd(k >> 2) for k = 1:486] for _ = 1:3],
        (1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0),
        123,
    )
    @compile_workload begin
        for (system, prn, symbols) in (
            (GPSL1CA(), 25, gps_l1ca),
            (GalileoE1B(), 2, galileo_e1b),
            (GalileoE1B_BOC11(), 2, galileo_e1b),
            (GPSL2CM(), 25, gps_l2c),
            (GPSL5I(), 25, gps_l5i),
            (GalileoE5aI(), 21, galileo_e5a),
            # E5b carries the same I/NAV pages as E1B and decodes them
            # identically, so the E1B stream exercises it for real.
            (GalileoE5bI(), 2, galileo_e1b),
            (GalileoE6B(), 2, galileo_e6b),
            (GPSL1C_D(), 7, random_symbols),
            (BeiDouB1I(), 7, random_symbols),
            (BeiDouB1C_D(), 7, random_symbols),
            (BeiDouB2aI(), 7, random_symbols),
            (BeiDouB2bI(), 7, random_symbols),
            (BeiDouB3I(), 7, random_symbols),
        )
            state = decode(GNSSDecoderState(system, prn), symbols, length(symbols))
            is_sat_healthy(state)
            is_decoding_completed_for_positioning(state)
        end
    end
end
