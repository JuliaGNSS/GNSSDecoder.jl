# GPS L1C-D (CNAV-2) decoder — TOI sync + subframe 2 (issue #38).
#
# The L1C-D civil data signal carries the CNAV-2 message as 18-second frames
# of 1800 symbols at 100 sps:
#
#   - Subframe 1: 52 symbols — a BCH(51,8)+parity encoding of a 9-bit TOI
#     count (0..399), used purely for frame sync (see `src/bch_toi.jl`).
#   - Subframes 2 + 3: 1748 symbols, block-interleaved (38×46, see
#     `src/deinterleave.jl`) over a 1200-symbol rate-½ LDPC codeword for
#     subframe 2 (600 info bits) followed by a 548-symbol rate-½ LDPC
#     codeword for subframe 3 (274 info bits).
#
# This slice decodes subframe 2 fully (clock + ephemeris, IS-GPS-800G
# Figure 3.5-1 / Table 3.5-1) and LDPC-decodes + CRC-checks subframe 3 but
# only records that a page was received; the per-page subframe-3 field
# parsing lands in issue #39.
#
# Frame sync follows the generic streaming framework in `src/gnss.jl`: the
# soft-symbol `CircularDeque` is sized to one full subframe (1800) plus the
# 52-symbol subframe-1 segment of the *next* frame (1852 total) so the BCH
# match can be confirmed at both ends of the window before locking.

# `import Aff3ct` for the qualified `Aff3ct.decode` call (it would clash with
# `GNSSDecoder.decode`); bring only the LDPC constructors into scope by name.
import Aff3ct
using Aff3ct: LDPCMatrix, LDPCBPDecoder

"Total symbols accumulated before a sync attempt: 1800 (one frame) + 52 (next subframe-1)."
const L1C_D_FRAME_LENGTH = 1800
const L1C_D_SUBFRAME1_LENGTH = 52
const L1C_D_WINDOW_LENGTH = L1C_D_FRAME_LENGTH + L1C_D_SUBFRAME1_LENGTH  # 1852

# Subframe 2/3 channel-symbol counts (IS-GPS-800G §3.2.3).
const L1C_D_SF2_SYMBOLS = 1200
const L1C_D_SF3_SYMBOLS = 548
const L1C_D_PAYLOAD_SYMBOLS = L1C_D_SF2_SYMBOLS + L1C_D_SF3_SYMBOLS  # 1748
const L1C_D_INTERLEAVER_ROWS = 38
const L1C_D_INTERLEAVER_COLS = 46  # 38 * 46 == 1748

# LDPC info-block sizes (rate ½): N channel symbols -> K info bits.
const L1C_D_SF2_INFO_BITS = 600
const L1C_D_SF3_INFO_BITS = 274
const L1C_D_CRC_BITS = 24

# Semi-major axis reference (IS-GPS-800G Table 3.5-1 footnote, meters).
const L1C_D_A_REF = 26_559_710.0
# Rate-of-right-ascension reference (IS-GPS-800G Table 3.5-1, semi-circles/sec).
const L1C_D_OMEGA_DOT_REF = -2.6e-9

"""
    GPSL1C_DConstants

WGS 84 constants and CNAV-2 message structure parameters for GPS L1C-D decoding.

The frame is modelled through the generic streaming framework: `preamble_length`
is the 52-symbol subframe-1 BCH segment of the *next* frame retained at the tail
of the sync window, and `syncro_sequence_length` is the 1800-symbol frame that is
drained once a subframe is decoded.

# Fields
$(TYPEDFIELDS)

# Reference
IS-GPS-800G, Sections 3.2 and 3.5, Table 3.5-1.
"""
Base.@kwdef struct GPSL1C_DConstants <: AbstractGNSSConstants
    "Frame length drained after each decoded subframe (1800 symbols)"
    syncro_sequence_length::Int = L1C_D_FRAME_LENGTH
    "Trailing next-frame subframe-1 BCH segment retained for sync (52 symbols)"
    preamble_length::Int = L1C_D_SUBFRAME1_LENGTH
    "Mathematical constant π (IS-GPS-800G)"
    PI::Float64 = 3.1415926535898
    "WGS 84 Earth rotation rate (rad/s)"
    Ω_dot_e::Float64 = 7.2921151467e-5
    "Speed of light (m/s)"
    c::Float64 = 2.99792458e8
    "WGS 84 Earth gravitational parameter (m³/s²)"
    μ::Float64 = 3.986005e14
    "Relativistic correction constant (s/√m)"
    F::Float64 = -4.442807633e-10
end

"""
    GPSL1C_DData

Decoded GPS L1C-D (CNAV-2) navigation message data.

Holds the subframe-2 clock, ephemeris, and accuracy parameters (IS-GPS-800G
Figure 3.5-1 / Table 3.5-1). Subframe-3 page contents are not parsed in this
slice (issue #39); only the count of CRC-valid subframe-3 pages received is
tracked. Field-naming follows [`GPSL1CAData`](@ref): semi-circle quantities are
converted to radians on decode (multiplied by π), all `Union{Nothing,…}` until
first decoded.

# Sync / timing
- `toi::Int`: Last validated Time-Of-Interval count (0..399), or `nothing`.
- `ITOW::Int64`: Interval time of week — number of two-hour epochs since the
  start of the week (subframe 2 bits 14-21).
- `WN::Int64`: Transmission week number, modulo-8192 (subframe 2 bits 1-13).
- `t_op::Int64`: Data predict time of week (seconds).

# Health / accuracy
- `l1c_health::Bool`: L1C signal health bit (false = OK, true = bad/unavailable).
- `ura_ed_index::Int64`: Ephemeris URA index (signed).
- `ura_ned0_index::Int64`, `ura_ned1_index::Int64`, `ura_ned2_index::Int64`:
  Clock URA indices.

# Ephemeris (Table 3.5-1)
- `t_0e::Int64`: Ephemeris/clock data reference time of week (seconds).
- `ΔA::Float64`: Semi-major axis difference at reference time (meters).
- `A_dot::Float64`: Change rate in semi-major axis (m/s).
- `Δn_0::Float64`: Mean motion difference from computed value (rad/s).
- `Δn_0_dot::Float64`: Rate of mean motion difference (rad/s²).
- `M_0::Float64`: Mean anomaly at reference time (rad).
- `e::Float64`: Eccentricity (dimensionless).
- `ω::Float64`: Argument of perigee (rad).
- `Ω_0::Float64`: Reference right ascension angle (rad).
- `i_0::Float64`: Inclination angle at reference time (rad).
- `ΔΩ_dot::Float64`: Rate of right ascension difference (rad/s).
- `i_0_dot::Float64`: Rate of inclination angle (rad/s).
- `C_is::Float64`, `C_ic::Float64`: Sine/cosine inclination harmonic corrections (rad).
- `C_rs::Float64`, `C_rc::Float64`: Sine/cosine orbit-radius harmonic corrections (m).
- `C_us::Float64`, `C_uc::Float64`: Sine/cosine argument-of-latitude harmonic corrections (rad).

# Clock (Table 3.5-1)
- `t_0c::Int64`: Clock data reference time of week (seconds); equals `t_0e` in CNAV-2.
- `a_f0::Float64`, `a_f1::Float64`, `a_f2::Float64`: Clock bias / drift / drift-rate.
- `T_GD::Float64`: L1/L2 P(Y) inter-signal correction (seconds).
- `ISC_L1CP::Float64`, `ISC_L1CD::Float64`: L1CP / L1CD inter-signal corrections (seconds).

# Subframe 3 (issue #39 deferral)
- `num_sf3_pages_received::Int`: Count of CRC-valid subframe-3 pages received.

# Reference
IS-GPS-800G, Figure 3.5-1 and Table 3.5-1.
"""
Base.@kwdef struct GPSL1C_DData <: AbstractGNSSData
    toi::Union{Nothing,Int} = nothing
    ITOW::Union{Nothing,Int64} = nothing
    WN::Union{Nothing,Int64} = nothing
    t_op::Union{Nothing,Int64} = nothing

    l1c_health::Union{Nothing,Bool} = nothing
    ura_ed_index::Union{Nothing,Int64} = nothing
    ura_ned0_index::Union{Nothing,Int64} = nothing
    ura_ned1_index::Union{Nothing,Int64} = nothing
    ura_ned2_index::Union{Nothing,Int64} = nothing

    t_0e::Union{Nothing,Int64} = nothing
    ΔA::Union{Nothing,Float64} = nothing
    A_dot::Union{Nothing,Float64} = nothing
    Δn_0::Union{Nothing,Float64} = nothing
    Δn_0_dot::Union{Nothing,Float64} = nothing
    M_0::Union{Nothing,Float64} = nothing
    e::Union{Nothing,Float64} = nothing
    ω::Union{Nothing,Float64} = nothing
    Ω_0::Union{Nothing,Float64} = nothing
    i_0::Union{Nothing,Float64} = nothing
    ΔΩ_dot::Union{Nothing,Float64} = nothing
    i_0_dot::Union{Nothing,Float64} = nothing
    C_is::Union{Nothing,Float64} = nothing
    C_ic::Union{Nothing,Float64} = nothing
    C_rs::Union{Nothing,Float64} = nothing
    C_rc::Union{Nothing,Float64} = nothing
    C_us::Union{Nothing,Float64} = nothing
    C_uc::Union{Nothing,Float64} = nothing

    t_0c::Union{Nothing,Int64} = nothing
    a_f0::Union{Nothing,Float64} = nothing
    a_f1::Union{Nothing,Float64} = nothing
    a_f2::Union{Nothing,Float64} = nothing
    T_GD::Union{Nothing,Float64} = nothing
    ISC_L1CP::Union{Nothing,Float64} = nothing
    ISC_L1CD::Union{Nothing,Float64} = nothing

    num_sf3_pages_received::Int = 0
end

function GPSL1C_DData(
    data::GPSL1C_DData;
    toi = data.toi,
    ITOW = data.ITOW,
    WN = data.WN,
    t_op = data.t_op,
    l1c_health = data.l1c_health,
    ura_ed_index = data.ura_ed_index,
    ura_ned0_index = data.ura_ned0_index,
    ura_ned1_index = data.ura_ned1_index,
    ura_ned2_index = data.ura_ned2_index,
    t_0e = data.t_0e,
    ΔA = data.ΔA,
    A_dot = data.A_dot,
    Δn_0 = data.Δn_0,
    Δn_0_dot = data.Δn_0_dot,
    M_0 = data.M_0,
    e = data.e,
    ω = data.ω,
    Ω_0 = data.Ω_0,
    i_0 = data.i_0,
    ΔΩ_dot = data.ΔΩ_dot,
    i_0_dot = data.i_0_dot,
    C_is = data.C_is,
    C_ic = data.C_ic,
    C_rs = data.C_rs,
    C_rc = data.C_rc,
    C_us = data.C_us,
    C_uc = data.C_uc,
    t_0c = data.t_0c,
    a_f0 = data.a_f0,
    a_f1 = data.a_f1,
    a_f2 = data.a_f2,
    T_GD = data.T_GD,
    ISC_L1CP = data.ISC_L1CP,
    ISC_L1CD = data.ISC_L1CD,
    num_sf3_pages_received = data.num_sf3_pages_received,
)
    GPSL1C_DData(
        toi,
        ITOW,
        WN,
        t_op,
        l1c_health,
        ura_ed_index,
        ura_ned0_index,
        ura_ned1_index,
        ura_ned2_index,
        t_0e,
        ΔA,
        A_dot,
        Δn_0,
        Δn_0_dot,
        M_0,
        e,
        ω,
        Ω_0,
        i_0,
        ΔΩ_dot,
        i_0_dot,
        C_is,
        C_ic,
        C_rs,
        C_rc,
        C_us,
        C_uc,
        t_0c,
        a_f0,
        a_f1,
        a_f2,
        T_GD,
        ISC_L1CP,
        ISC_L1CD,
        num_sf3_pages_received,
    )
end

"""
$(TYPEDEF)

Per-decoder cache for the GPS L1C-D signal.

Holds the soft-symbol `CircularDeque{Float32}` (capacity = 1852 = 1800 frame +
52 next-frame BCH segment) and the two Aff3ct LDPC belief-propagation decoders
(subframe 2: K=600, N=1200; subframe 3: K=274, N=548), built lazily from the
committed `.alist` parity matrices in `data/`. The decoder objects are mutable
Aff3ct handles; they are shared by reference through the otherwise-immutable
[`GNSSDecoderState`](@ref).

# Fields
$(TYPEDFIELDS)
"""
struct GPSL1C_DCache <: AbstractGNSSCache
    "Soft-symbol buffer (1852 = 1800 frame + 52 next subframe-1)"
    soft_buffer::CircularDeque{Float32}
    "Aff3ct LDPC BP decoder for subframe 2 (K=600, N=1200)"
    sf2_decoder::LDPCBPDecoder
    "Aff3ct LDPC BP decoder for subframe 3 (K=274, N=548)"
    sf3_decoder::LDPCBPDecoder
end

# Path to the committed LDPC `.alist` parity matrices. `pkgdir`-free: walk up
# from this file (src/gps/) to the package root, then into data/.
_l1c_d_data_path(name) = joinpath(@__DIR__, "..", "..", "data", name)

function GPSL1C_DCache()
    sf2_H = LDPCMatrix(_l1c_d_data_path("cnv2_sf2.alist"))
    sf3_H = LDPCMatrix(_l1c_d_data_path("cnv2_sf3.alist"))
    GPSL1C_DCache(
        CircularDeque{Float32}(L1C_D_WINDOW_LENGTH),
        LDPCBPDecoder(sf2_H; num_iterations = 50),
        LDPCBPDecoder(sf3_H; num_iterations = 50),
    )
end

# The LDPC decoder handles are stateless w.r.t. equality (they are runtime
# Aff3ct objects); two L1C-D caches are equal when their soft buffers match.
function Base.:(==)(a::GPSL1C_DCache, b::GPSL1C_DCache)
    deques_equal(a.soft_buffer, b.soft_buffer)
end

function is_subframe2_decoded(data::GPSL1C_DData)
    !isnothing(data.WN) &&
        !isnothing(data.ITOW) &&
        !isnothing(data.t_0e) &&
        !isnothing(data.ΔA) &&
        !isnothing(data.M_0) &&
        !isnothing(data.e) &&
        !isnothing(data.ω) &&
        !isnothing(data.Ω_0) &&
        !isnothing(data.i_0) &&
        !isnothing(data.a_f0) &&
        !isnothing(data.a_f1) &&
        !isnothing(data.l1c_health)
end

function is_decoding_completed_for_positioning(data::GPSL1C_DData)
    !isnothing(data.toi) && is_subframe2_decoded(data)
end

"""
$(TYPEDSIGNATURES)

Create a decoder state for GPS L1C-D (CNAV-2) navigation messages.

Wires up a [`GNSSDecoderState`](@ref) with a 1852-symbol soft-symbol buffer,
the 400-entry BCH(51,8) TOI codeword table (`src/bch_toi.jl`), and two Aff3ct
LDPC belief-propagation decoders loaded lazily from the committed `.alist`
parity matrices in `data/`.

# Arguments
- `prn::Int`: Pseudo-Random Noise code identifier (1-63 for L1C).

# Returns
- `GNSSDecoderState{GPSL1C_DData}`: Initialized decoder state for GPS L1C-D.

# Example
```julia
state = GPSL1C_DDecoderState(1)            # PRN 1
state = decode(state, soft_symbols, num_symbols)
if is_sat_healthy(state)
    # Use state.data for positioning
end
```

# See Also
- [`GNSSDecoderState`](@ref): The underlying state structure
- [`decode`](@ref): Decode soft symbols using this state
- [`reset_decoder_state`](@ref): Reset after signal loss
"""
function GPSL1C_DDecoderState(prn)
    GNSSDecoderState(
        prn,
        GPSL1C_DData(),
        GPSL1C_DData(),
        GPSL1C_DConstants(),
        GPSL1C_DCache(),
        nothing,
        false,
    )
end

function GNSSDecoderState(system::GPSL1C_D, prn)
    GPSL1C_DDecoderState(prn)
end

"""
$(TYPEDSIGNATURES)

Reset the GPS L1C-D decoder state after a signal loss or reacquisition.

Clears the in-flight sync state (soft-symbol buffer and TOI) and the validated
data, while preserving the long-lived CED/clock fields in `raw_data` so a
[`GNSSReceiver`] can re-use the satellite after reacquisition without re-decoding
all of subframe 2. Mirrors the semantics of the GPS L1 C/A and Galileo E1B
implementations.

# Arguments
- `state::GNSSDecoderState{<:GPSL1C_DData}`: Current GPS L1C-D decoder state.

# Returns
- `GNSSDecoderState{<:GPSL1C_DData}`: Reset decoder state with cleared buffers.

# See Also
- [`GPSL1C_DDecoderState`](@ref): Create a fresh decoder state
- [`decode`](@ref): Continue decoding after reset
"""
function reset_decoder_state(state::GNSSDecoderState{<:GPSL1C_DData})
    empty!(state.cache.soft_buffer)
    GNSSDecoderState(
        state;
        raw_data = GPSL1C_DData(state.raw_data; toi = nothing),
        data = GPSL1C_DData(),
        num_bits_after_valid_syncro_sequence = nothing,
        is_shifted_by_180_degrees = false,
    )
end

# ---- Sync ------------------------------------------------------------------
#
# Override the generic packed-buffer sync. L1C-D sync runs the BCH(51,8) TOI
# match directly on the soft-symbol deque: the first 52 symbols must equal the
# codeword for some `toi`, and the last 52 the codeword for `(toi+1) mod 400`,
# in either polarity. `sync_bch_toi` (src/bch_toi.jl) implements exactly this.

# `CircularDeque` indexing is O(1); collect a contiguous slice into a `Vector`
# for the (length-checked) `soft_to_hard_codeword`. Small (52 elements), runs
# only once per buffered window.
function _deque_slice(deque::CircularDeque{Float32}, start::Int, len::Int)
    v = Vector{Float32}(undef, len)
    @inbounds for i in 1:len
        v[i] = deque[start + i - 1]
    end
    return v
end

"""
    try_sync(state::GNSSDecoderState{<:GPSL1C_DData}) -> Union{Nothing,BCHToiSync}

Per-signal sync hook for GPS L1C-D. Hard-slices the leading 52 and trailing 52
soft symbols of the 1852-symbol window and runs the multi-subframe BCH(51,8)
match (`sync_bch_toi`). Returns the [`BCHToiSync`](@ref) on a match (carrying
the TOI and the detected polarity flip) or `nothing`.
"""
function try_sync(state::GNSSDecoderState{<:GPSL1C_DData})
    deque = soft_buffer(state)
    first52 = soft_to_hard_codeword(_deque_slice(deque, 1, L1C_D_SUBFRAME1_LENGTH))
    next52 = soft_to_hard_codeword(
        _deque_slice(deque, L1C_D_FRAME_LENGTH + 1, L1C_D_SUBFRAME1_LENGTH),
    )
    sync_bch_toi(first52, next52)
end

"""
    complement_buffer_if_necessary(state::GNSSDecoderState{<:GPSL1C_DData}, sync)

Resolve the 180° polarity ambiguity for L1C-D. `sync::BCHToiSync` already
carries the detected polarity from `sync_bch_toi`; record it on the state and
pass the sync object through unchanged for `decode_syncro_sequence`.
"""
function complement_buffer_if_necessary(
    state::GNSSDecoderState{<:GPSL1C_DData},
    sync::BCHToiSync,
)
    GNSSDecoderState(state; is_shifted_by_180_degrees = sync.polarity_flipped), sync
end

# ---- Subframe pipeline -----------------------------------------------------

"""
    decode_syncro_sequence(state::GNSSDecoderState{<:GPSL1C_DData}, sync::BCHToiSync)

Process one locked CNAV-2 frame. Enforces a monotonic TOI (each frame's TOI must
be `(previous + 1) mod 400`; a break resets the decoder to searching), then
deinterleaves and LDPC-decodes subframes 2 and 3 from the 1748-symbol payload,
CRC-checking each. Subframe 2 fields are parsed into `raw_data`; subframe 3 is
recorded as a received page (issue #39 will parse it).
"""
function decode_syncro_sequence(state::GNSSDecoderState{<:GPSL1C_DData}, sync::BCHToiSync)
    # Monotonic-TOI validation: a locked frame whose TOI does not follow the
    # previous one by exactly +1 (mod 400) indicates loss of frame lock.
    prev_toi = state.raw_data.toi
    if !isnothing(prev_toi) && (prev_toi + 1) % TOI_RANGE != sync.toi
        return reset_decoder_state(state)
    end
    state = GNSSDecoderState(state; raw_data = GPSL1C_DData(state.raw_data; toi = sync.toi))

    # Extract the 1748-symbol interleaved SF2+SF3 payload (symbols 53..1800),
    # applying the polarity flip by negating soft symbols up front.
    deque = soft_buffer(state)
    flip = state.is_shifted_by_180_degrees ? -1.0f0 : 1.0f0
    interleaved = Vector{Float32}(undef, L1C_D_PAYLOAD_SYMBOLS)
    @inbounds for i in 1:L1C_D_PAYLOAD_SYMBOLS
        interleaved[i] = flip * deque[L1C_D_SUBFRAME1_LENGTH + i]
    end

    deinterleaved = deinterleave(interleaved, L1C_D_INTERLEAVER_ROWS, L1C_D_INTERLEAVER_COLS)
    sf2_symbols = @view deinterleaved[1:L1C_D_SF2_SYMBOLS]
    sf3_symbols = @view deinterleaved[(L1C_D_SF2_SYMBOLS + 1):L1C_D_PAYLOAD_SYMBOLS]

    state = decode_subframe2(state, sf2_symbols)
    state = decode_subframe3(state, sf3_symbols)
    return state
end

"Run an Aff3ct LDPC BP decode on `symbols` and return the `info_bits`-long info block as a `Vector{Bool}` (true = 1)."
function ldpc_decode_bits(decoder::LDPCBPDecoder, symbols, info_bits::Int)
    # AFF3CT LLR convention matches ours: positive ⇒ bit 0, negative ⇒ bit 1.
    llr = collect(Float32, symbols)
    info = Aff3ct.decode(decoder, llr)
    bits = Vector{Bool}(undef, info_bits)
    @inbounds for i in 1:info_bits
        bits[i] = info[i] != 0
    end
    return bits
end

"Verify a CRC-24Q-protected info block (message bits followed by a 24-bit CRC). Returns true iff the CRC matches."
function crc24q_ok(info_bits::AbstractVector{Bool})
    crc24q(info_bits) == 0
end

# ---- Subframe 2 bit-field extraction (IS-GPS-800G Figure 3.5-1) ------------
#
# `bits` is the 600-bit subframe-2 info block, MSB first, 1-based. Helpers
# below read fields by their 1-based start bit and length.

"Read an unsigned integer field of `len` bits starting at 1-based bit `start` (MSB first)."
function _u(bits::AbstractVector{Bool}, start::Int, len::Int)
    v = 0
    @inbounds for i in 0:(len - 1)
        v = (v << 1) | (bits[start + i] ? 1 : 0)
    end
    return v
end

"Read a two's-complement signed integer field of `len` bits starting at 1-based bit `start`."
function _s(bits::AbstractVector{Bool}, start::Int, len::Int)
    v = _u(bits, start, len)
    (bits[start] ? v - (1 << len) : v)
end

function decode_subframe2(state::GNSSDecoderState{<:GPSL1C_DData}, sf2_symbols)
    bits = ldpc_decode_bits(state.cache.sf2_decoder, sf2_symbols, L1C_D_SF2_INFO_BITS)
    crc24q_ok(bits) || return state  # silently drop on CRC failure

    PI = state.constants.PI

    WN = _u(bits, 1, 13)
    ITOW = _u(bits, 14, 8)
    t_op = _u(bits, 22, 11) * 300
    l1c_health = bits[33]
    ura_ed_index = _s(bits, 34, 5)
    t_0e = _u(bits, 39, 11) * 300
    ΔA = _s(bits, 50, 26) * 2.0^-9
    A_dot = _s(bits, 76, 25) * 2.0^-21
    Δn_0 = _s(bits, 101, 17) * 2.0^-44 * PI
    Δn_0_dot = _s(bits, 118, 23) * 2.0^-57 * PI
    M_0 = _s(bits, 141, 33) * 2.0^-32 * PI
    e = _u(bits, 174, 33) * 2.0^-34
    ω = _s(bits, 207, 33) * 2.0^-32 * PI
    Ω_0 = _s(bits, 240, 33) * 2.0^-32 * PI
    i_0 = _s(bits, 273, 33) * 2.0^-32 * PI
    ΔΩ_dot = _s(bits, 306, 17) * 2.0^-44 * PI
    i_0_dot = _s(bits, 323, 15) * 2.0^-44 * PI
    C_is = _s(bits, 338, 16) * 2.0^-30
    C_ic = _s(bits, 354, 16) * 2.0^-30
    C_rs = _s(bits, 370, 24) * 2.0^-8
    C_rc = _s(bits, 394, 24) * 2.0^-8
    C_us = _s(bits, 418, 21) * 2.0^-30
    C_uc = _s(bits, 439, 21) * 2.0^-30
    ura_ned0_index = _s(bits, 460, 5)
    ura_ned1_index = _u(bits, 465, 3)
    ura_ned2_index = _u(bits, 468, 3)
    a_f0 = _s(bits, 471, 26) * 2.0^-35
    a_f1 = _s(bits, 497, 20) * 2.0^-48
    a_f2 = _s(bits, 517, 10) * 2.0^-60
    T_GD = _s(bits, 527, 13) * 2.0^-35
    ISC_L1CP = _s(bits, 540, 13) * 2.0^-35
    ISC_L1CD = _s(bits, 553, 13) * 2.0^-35

    raw = GPSL1C_DData(
        state.raw_data;
        WN,
        ITOW,
        t_op,
        l1c_health,
        ura_ed_index,
        t_0e,
        ΔA,
        A_dot,
        Δn_0,
        Δn_0_dot,
        M_0,
        e,
        ω,
        Ω_0,
        i_0,
        ΔΩ_dot,
        i_0_dot,
        C_is,
        C_ic,
        C_rs,
        C_rc,
        C_us,
        C_uc,
        ura_ned0_index,
        ura_ned1_index,
        ura_ned2_index,
        t_0c = t_0e,  # CNAV-2 shares one reference time for clock and ephemeris
        a_f0,
        a_f1,
        a_f2,
        T_GD,
        ISC_L1CP,
        ISC_L1CD,
    )
    GNSSDecoderState(state; raw_data = raw)
end

function decode_subframe3(state::GNSSDecoderState{<:GPSL1C_DData}, sf3_symbols)
    bits = ldpc_decode_bits(state.cache.sf3_decoder, sf3_symbols, L1C_D_SF3_INFO_BITS)
    # A subframe-3 page is "received" regardless of CRC outcome (issue #38
    # records the page; field parsing — and any CRC-gated trust — is #39).
    crc24q_ok(bits)  # currently informational only
    GNSSDecoderState(
        state;
        raw_data = GPSL1C_DData(
            state.raw_data;
            num_sf3_pages_received = state.raw_data.num_sf3_pages_received + 1,
        ),
    )
end

"""
    validate_data(state::GNSSDecoderState{<:GPSL1C_DData})

Promote `raw_data` to `data` once subframe 2 has been fully decoded. The
per-frame BCH re-check and monotonic-TOI enforcement are performed inline by
[`decode_syncro_sequence`](@ref); this hook only publishes validated data and
arms the streaming counter.
"""
function validate_data(state::GNSSDecoderState{<:GPSL1C_DData})
    if is_decoding_completed_for_positioning(state.raw_data)
        return GNSSDecoderState(
            state;
            data = state.raw_data,
            num_bits_after_valid_syncro_sequence = state.constants.preamble_length,
        )
    end
    return state
end

"""
$(TYPEDSIGNATURES)

Check if the GPS L1C-D satellite is healthy and usable for positioning.

Examines the 1-bit L1C signal health flag from subframe 2 (IS-GPS-800G
§3.5.3.4): a satellite is healthy iff the health bit is 0 (Signal OK).

!!! warning
    Requires subframe 2 to have been decoded; returns `false` until then.

# Arguments
- `state::GNSSDecoderState{<:GPSL1C_DData}`: GPS L1C-D decoder state.

# Returns
- `Bool`: `true` iff the L1C signal-health bit indicates OK.
"""
function is_sat_healthy(state::GNSSDecoderState{<:GPSL1C_DData})
    state.data.l1c_health === false
end
