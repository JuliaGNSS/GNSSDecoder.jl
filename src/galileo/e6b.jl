# Galileo E6-B (C/NAV) navigation message decoder — Galileo HAS SIS ICD,
# Issue 1.0 (in force; the E6-B/C spreading codes are specified separately in
# the Galileo E6-B/C Codes Technical Note, Issue 1).
#
# E6-B is the data component of the E6 signal at 1278.75 MHz; E6-C is a
# dataless pilot. E6-B carries the C/NAV message at 1000 sps, which is the
# broadcast channel of the Galileo High Accuracy Service (HAS): PPP-grade orbit,
# clock, code-bias and phase-bias corrections *to another signal's* broadcast
# ephemeris, rather than an ephemeris of its own.
#
# One C/NAV page per second, 1000 symbols (ICD Table 5):
#
#     Sync (16 sym, 1011011101110000) | 984 encoded symbols
#
# The 984 encoded symbols are the same Galileo K=7 rate-1/2 NSC convolutional
# code every other Galileo data channel uses (G1 = 0o171, G2 = 0o133, G2
# inverted, ICD Table 3 — explicitly "the same as the one used for other Galileo
# data channels"), over a 123 x 8 block interleaver (ICD Table 4). They decode
# to 492 trellis steps = 486 information bits + 6 tail bits:
#
#     Reserved (14) | HAS Page (448) | CRC (24)      [+ 6 tail bits]
#
# with CRC-24Q over the leading 462 bits (ICD §2.3.3) — the same polynomial as
# I/NAV, so `crc24q` is reused. The HAS Page is a 24-bit HAS Page Header plus a
# 424-bit HAS Encoded Page (ICD Table 6).
#
# THE OUTER LAYER is what makes E6-B unusual. A HAS message is *not* carried by
# one satellite: it is cut into k <= 32 non-encoded pages of 424 bits, encoded
# "vertically" with a systematic RS(255, 32) code over GF(256) into 255 encoded
# pages, and different encoded pages are handed to different satellites (ICD
# §6.3, "HPVRS"). A receiver collects any k pages carrying the same Message ID —
# from one satellite over k seconds, or from several satellites at once — and
# recovers the message by inverting the k x k submatrix of the generator matrix
# selected by the received Page IDs (ICD §6.4). That inversion, the field
# arithmetic, and the generator matrix live in `src/coding/reed_solomon.jl`.
#
# Consequently a *single-satellite* decoder state accumulates pages across time
# and only occasionally completes a message; a receiver tracking several E6-B
# satellites will complete messages much faster by pooling their pages, which
# this API supports by exposing the per-page store on the cache. Messages that
# do not complete within 150 s are discarded (ICD §6.4.1).
#
# Only Message Type 1 is defined (ICD Table 10); its body is a variable-length
# sequence of content blocks selected by six flags in the MT1 header (ICD
# Table 12). Block lengths depend on the *mask* — how many constellations,
# satellites and signals are corrected — which is itself broadcast in a mask
# block, possibly in a different message. A message whose Mask Flag is 0
# therefore cannot be parsed until a mask with its Mask ID has been received;
# such a message is held and re-parsed when a matching mask arrives.
#
# CROSS-CHECKS. The framing, header layout, CRC scope and interleaver
# orientation agree with GNSS-SDR (`Galileo_CNAV.h`, `galileo_cnav_message.cc`,
# `galileo_e6_has_msg_receiver.cc`) and with PocketSDR's C/NAV framing
# (`decode_E6B` / `decode_gal_syms(buff+16, 123, 8, bits)`); PocketSDR stops at
# the 486 raw bits and implements no HAS layer. Two deliberate differences from
# GNSS-SDR are noted at their sites: it recovers the message with a
# Berlekamp-Massey erasure decode of the full 255-symbol codeword instead of the
# ICD's generator-matrix inversion, and it maps the ICD's "data not available"
# and "shall not be used" sentinels to a 0 m correction downstream, where this
# decoder keeps them distinguishable (`nothing` and a flag).

# ---- Page framing (HAS SIS ICD, Issue 1.0, §2.3) ----------------------------

"""
One C/NAV page: 1000 symbols at 1000 sps = 1 second (ICD Table 5).
"""
const E6B_PAGE_SYMBOLS = 1000
"""
Leading synchronisation pattern `1011011101110000` = 0xB770 (ICD §2.3.1),
16 symbols, neither encoded nor interleaved.
"""
const E6B_SYNC_PATTERN = 0xB770
const E6B_SYNC_SYMBOLS = 16
"""
FEC-encoded symbols per page: 1000 - 16 (ICD Table 5).
"""
const E6B_ENCODED_SYMBOLS = 984
"""
Information bits recovered per page: 984 / 2 - 6 tail bits (ICD Table 5 counts
the 6 tail bits inside its 492).
"""
const E6B_PAGE_BITS = 486
"""
Columns of the C/NAV block interleaver (ICD Table 4: 8 rows by 123 columns); the
8 rows every Galileo channel shares are `GALILEO_INTERLEAVER_ROWS`.
"""
const E6B_INTERLEAVER_COLUMNS = 123
"""
Sync window: one page plus the next page's sync pattern.
"""
const E6B_WINDOW_SYMBOLS = E6B_PAGE_SYMBOLS + E6B_SYNC_SYMBOLS  # 1016
"""
Octets spanned by the 486 page bits, for the packed-word `crc24q` — 486 rounds up
to 61, whose two leading bits are zero and therefore neutral.
"""
const E6B_PAGE_OCTETS = cld(E6B_PAGE_BITS, 8)
"""
1-based bit position of the 24-bit HAS Page Header inside the 486-bit page
(after the 14 reserved bits).
"""
const E6B_HEADER_START = 15
"""
1-based bit position of the 424-bit HAS Encoded Page inside the 486-bit page.
"""
const E6B_ENCODED_PAGE_START = 39
"""
Octets per HAS Encoded Page: 424 bits / 8 (the RS code's J, ICD §6.3).
"""
const E6B_OCTETS_PER_PAGE = 53
"""
HAS dummy-page header `hex[AF3BC3]` (ICD §2.4.1). Broadcast whenever a
satellite has no valid HAS data; such pages are discarded.
"""
const E6B_DUMMY_PAGE_HEADER = UInt32(0xAF3BC3)
"""
RS(255, 32) code parameters of the HPVRS outer layer (ICD §6.2).
"""
const E6B_RS_CODE_LENGTH = 255
const E6B_RS_CODE_DIMENSION = 32
"""
The systematic RS(255, 32) generator matrix of the HPVRS outer layer (ICD Annex
B), built once at load time from the narrow-sense generator polynomial.

Both its arguments are compile-time constants and `rs_erasure_decode` only reads
it, so it is shared by every decoder state rather than rebuilt and stored per
satellite — 8160 read-only octets and ~70 µs each otherwise.
"""
const E6B_GENERATOR_MATRIX = rs_systematic_generator_matrix(
    GALILEO_HAS_GF256,
    E6B_RS_CODE_LENGTH,
    E6B_RS_CODE_DIMENSION,
)
"""
A Message ID must be completed within 150 s or its pages are discarded
(ICD §6.4.1), counted here in *accepted* C/NAV pages.

One C/NAV page is one second, so on an unbroken stream the two are the same
count — but the clock is `GalileoE6BCache.page_counter`, which only advances on a
page that passed sync and CRC-24Q. Pages lost to a marginal signal therefore
stretch the window beyond 150 s of wall time (a group can outlive its deadline by
however long the decoder was not decoding). This is a deliberate approximation:
a symbol-domain decoder has no notion of "now", and erring towards keeping pages
costs a stale reassembly at worst, where erring the other way would drop
recoverable messages on every fade. `reset_decoder_state` clears the store
outright for the same reason, since an outage is the one case where the
approximation is badly wrong.
"""
const E6B_MESSAGE_TIMEOUT_PAGES = 150

"""
Validity Interval Index → validity interval in seconds (ICD Table 23). Index 15
is reserved; [`e6b_validity_interval`](@ref) reports it as `nothing`.
"""
const E6B_VALIDITY_INTERVALS =
    (5, 10, 15, 20, 30, 60, 90, 120, 180, 240, 300, 600, 900, 1800, 3600)

"""
    e6b_validity_interval(index) -> Union{Nothing,Int}

Validity interval in seconds for a 4-bit Validity Interval Index (ICD Table 23),
or `nothing` for the reserved index 15.
"""
e6b_validity_interval(index::Integer) = get(E6B_VALIDITY_INTERVALS, index + 1, nothing)

# GNSS ID values (ICD Table 18) and the width of the Reference IOD field each
# implies (ICD Table 26). GNSS IDs 1 and 3-15 are reserved; a mask naming one
# makes every downstream block length unknowable, so the body parse stops.
const E6B_GNSS_ID_GPS = 0
const E6B_GNSS_ID_GALILEO = 2

"""
    e6b_iod_ref_length(GNSS_ID) -> Union{Nothing,Int}

Width in bits of the Reference IOD field for a GNSS (ICD Table 26): 8 for GPS
(IODE/IODC), 10 for Galileo (IODnav). `nothing` for a reserved GNSS ID.
"""
e6b_iod_ref_length(GNSS_ID::Integer) =
    GNSS_ID == E6B_GNSS_ID_GPS ? 8 : GNSS_ID == E6B_GNSS_ID_GALILEO ? 10 : nothing

# ---- HAS status --------------------------------------------------------------

"""
    HASStatus

Galileo High Accuracy Service status, broadcast in every HAS Page Header.

# Values

  - `has_test_mode`: HAS service testing activities ongoing; nominal performance may not be met (value 0)
  - `has_operational_mode`: HAS is expected to provide nominal performance (value 1)
  - `has_status_reserved`: Reserved (value 2)
  - `has_do_not_use`: Users shall stop using HAS from all satellites and discard previously received messages (value 3)

# Reference

Galileo HAS SIS ICD, Issue 1.0, Table 9
"""
@enum HASStatus begin
    has_test_mode
    has_operational_mode
    has_status_reserved
    has_do_not_use
end

# ---- Constants ---------------------------------------------------------------

"""
$(TYPEDEF)

Constants for the Galileo E6-B C/NAV (HAS) decoder — Galileo HAS SIS ICD,
Issue 1.0.

Unlike every other decoder in this package these carry *only* the sync geometry
the shared `decode` loop needs. C/NAV broadcasts corrections to another signal's
ephemeris, never an ephemeris of its own, so there is no satellite position to
compute here and no π, μ, F, speed of light or Earth rotation rate to compute it
with. Carrying them anyway would imply orbit math this decoder does not do —
uniformity of *interface* is already what `AbstractGNSSConstants` provides, and
the type lattice states the same distinction one level up (`GalileoE6BData` is an
`AbstractGalileoData` but not an `AbstractGalileoEphemerisData`).

# Fields

$(TYPEDFIELDS)
"""
Base.@kwdef struct GalileoE6BConstants <: AbstractGNSSConstants
    """
    Page length drained after each decoded page (1000 symbols)
    """
    syncro_sequence_length::Int = E6B_PAGE_SYMBOLS
    """
    Synchronisation pattern 0xB770 = `1011011101110000` (ICD §2.3.1), MSB first
    """
    preamble::UInt16 = E6B_SYNC_PATTERN
    """
    Trailing next-page sync segment retained for sync (16 symbols)
    """
    preamble_length::Int = E6B_SYNC_SYMBOLS
end

# ---- Storage bounds (HAS SIS ICD, Issue 1.0) ----------------------------------
#
# `decode!` allocates nothing, so every variable-length HAS structure is backed
# by a buffer sized once, at construction, to the largest value the ICD's field
# widths allow, and overwritten in place afterwards. These are those maxima.

"""
Most non-encoded pages in one HAS message — `MS` is a 5-bit "size - 1" field
(ICD Table 8), and equals the RS code dimension.
"""
const E6B_MAX_MESSAGE_PAGES = 32
"""
Most octets in one reassembled HAS message: 32 pages × 53 octets.
"""
const E6B_MAX_MESSAGE_OCTETS = E6B_MAX_MESSAGE_PAGES * 53
"""
Most bits in one reassembled HAS message (32 × 424) — the budget every
variable-length block is bounded by.
"""
const E6B_MAX_MESSAGE_BITS = 8 * E6B_MAX_MESSAGE_OCTETS
"""
Number of Message IDs (5-bit MID, ICD Table 8), i.e. of concurrent page groups.
"""
const E6B_NUM_MESSAGE_IDS = 32
"""
Number of Mask IDs (5-bit field, ICD Table 12).
"""
const E6B_NUM_MASK_IDS = 32
"""
Most constellations in one mask — `Nsys` is a 4-bit field (ICD Table 15).
"""
const E6B_MAX_SYSTEMS = 15
"""
Width of the Satellite Mask, i.e. most satellites per constellation (ICD Table 19).
"""
const E6B_MAX_SATELLITES = 40
"""
Width of the Signal Mask, i.e. most signals per constellation (ICD Table 20).
"""
const E6B_MAX_SIGNALS = 16

# ---- Decoded HAS records (ICD §5) -------------------------------------------
#
# A mask is small and bounded (at most 15 constellations of at most 40
# satellites and 16 signals), so it is an immutable, pointer-free *value*: its
# index lists and Cell Mask are stored inline rather than in heap vectors. That
# is what lets masks be copied into the preallocated Mask ID store, into a
# message and into the validated `data` without allocating, and without any two
# of them sharing a container.

"""
$(TYPEDEF)

The indices named by the set bits of a HAS bit mask (MSB = index 0), counting
from `first_index`: an immutable `AbstractVector{Int}` computed from the mask
itself, so it stores no list and never allocates.

It is what [`GalileoHASSatelliteMask`](@ref)'s `SVIDs` (40-bit Satellite Mask,
`first_index = 1`: Galileo SVID / GPS PRN) and `signal_indices` (16-bit Signal
Mask, `first_index = 0`) are. It compares equal to a `Vector` of the same
indices.

# Fields

$(TYPEDFIELDS)
"""
struct GalileoHASMaskIndices <: AbstractVector{Int}
    """
    The raw mask, right-aligned in `width` bits
    """
    mask::UInt64
    """
    Width of the mask in bits (40 or 16)
    """
    width::Int
    """
    Index named by the mask's most significant bit
    """
    first_index::Int
end

Base.size(indices::GalileoHASMaskIndices) = (count_ones(indices.mask),)
Base.IndexStyle(::Type{GalileoHASMaskIndices}) = IndexLinear()

# The bit position (0 = MSB) of the next set bit at or after `index`, or `width`.
@inline function _next_set_index(indices::GalileoHASMaskIndices, index::Int)
    while index < indices.width
        (indices.mask >> (indices.width - 1 - index)) & 0x1 == 0x1 && return index
        index += 1
    end
    return index
end

function Base.getindex(indices::GalileoHASMaskIndices, i::Int)
    @boundscheck checkbounds(indices, i)
    index = _next_set_index(indices, 0)
    for _ = 2:i
        index = _next_set_index(indices, index + 1)
    end
    return index + indices.first_index
end

function Base.iterate(indices::GalileoHASMaskIndices, index::Int = 0)
    index = _next_set_index(indices, index)
    index >= indices.width && return nothing
    return (index + indices.first_index, index + 1)
end

"""
$(TYPEDEF)

An `Nsat × Nsig` HAS Cell Mask (ICD §5.2.1.5) as an immutable, inline
`AbstractMatrix{Bool}`: one row of up to 16 bits per satellite, stored in a
fixed 40-row tuple so it never allocates. It compares equal to a
`Matrix{Bool}` with the same cells, and any `AbstractMatrix{Bool}` of at most
40 × 16 converts to it.

# Fields

$(TYPEDFIELDS)
"""
struct GalileoHASCellMask <: AbstractMatrix{Bool}
    """
    Row `r` holds satellite `r`'s cells as broadcast, `num_signals` bits
    right-aligned with signal 1 in the most significant of them; rows past
    `num_satellites` are zero
    """
    rows::NTuple{E6B_MAX_SATELLITES,UInt16}
    """
    Number of rows (`Nsat`)
    """
    num_satellites::Int
    """
    Number of columns (`Nsig`)
    """
    num_signals::Int
end

Base.size(cells::GalileoHASCellMask) = (cells.num_satellites, cells.num_signals)

function Base.getindex(cells::GalileoHASCellMask, row::Int, column::Int)
    @boundscheck checkbounds(cells, row, column)
    return (cells.rows[row] >> (cells.num_signals - column)) & 0x1 == 0x1
end

function GalileoHASCellMask(cells::AbstractMatrix{Bool})
    num_satellites, num_signals = size(cells)
    (num_satellites <= E6B_MAX_SATELLITES && num_signals <= E6B_MAX_SIGNALS) ||
        throw(ArgumentError("a HAS Cell Mask is at most 40 × 16"))
    rows = ntuple(Val(E6B_MAX_SATELLITES)) do row
        value = 0x0000
        if row <= num_satellites
            for column = 1:num_signals
                value = (value << 1) | UInt16(cells[row, column])
            end
        end
        value
    end
    GalileoHASCellMask(rows, num_satellites, num_signals)
end

Base.convert(::Type{GalileoHASCellMask}, cells::GalileoHASCellMask) = cells
Base.convert(::Type{GalileoHASCellMask}, cells::AbstractMatrix{Bool}) =
    GalileoHASCellMask(cells)

"""
    GalileoHASSatelliteMask

The HAS mask for one constellation: which satellites are corrected, which
signals carry biases, and which broadcast navigation message the orbit and clock
corrections refer to.

The raw ICD bit fields are kept alongside the expanded lists, since the raw
masks are what later blocks' lengths are computed from. An immutable,
pointer-free value: the expanded lists are views computed from the raw masks
(see [`GalileoHASMaskIndices`](@ref)), and the Cell Mask is stored inline.

# Fields

  - `GNSS_ID::Int`: GNSS index — 0 = GPS, 2 = Galileo (Table 18)
  - `satellite_mask::UInt64`: 40-bit Satellite Mask, MSB = satellite index 0 (Table 19)
  - `signal_mask::UInt16`: 16-bit Signal Mask, MSB = signal index 0 (Table 20)
  - `cell_mask::Union{Nothing,GalileoHASCellMask}`: `Nsat × Nsig` Cell Mask (an `AbstractMatrix{Bool}`), or `nothing` when the Cell Mask Availability Flag is 0 (biases then cover every masked satellite/signal pair)
  - `nav_message_index::Int`: Navigation Message Index — 0 = I/NAV (Galileo) / LNAV (GPS) (Table 21)
  - `SVIDs::GalileoHASMaskIndices`: Satellite IDs of the masked satellites (satellite index + 1, i.e. Galileo SVID / GPS PRN), an `AbstractVector{Int}` derived from `satellite_mask`
  - `signal_indices::GalileoHASMaskIndices`: Signal indices of the masked signals (0-based, per Table 20), an `AbstractVector{Int}` derived from `signal_mask`

# Reference

Galileo HAS SIS ICD, Issue 1.0, Tables 16-21
"""
Base.@kwdef struct GalileoHASSatelliteMask
    GNSS_ID::Int
    satellite_mask::UInt64
    signal_mask::UInt16
    cell_mask::Union{Nothing,GalileoHASCellMask} = nothing
    nav_message_index::Int
    SVIDs::GalileoHASMaskIndices = GalileoHASMaskIndices(satellite_mask, 40, 1)
    signal_indices::GalileoHASMaskIndices = GalileoHASMaskIndices(signal_mask, 16, 0)
end

Base.:(==)(a::GalileoHASSatelliteMask, b::GalileoHASSatelliteMask) = fields_equal(a, b)

"""
Placeholder filling the unused slots of a [`GalileoHASSatelliteMaskList`](@ref).
"""
const E6B_EMPTY_SATELLITE_MASK = GalileoHASSatelliteMask(;
    GNSS_ID = 0,
    satellite_mask = UInt64(0),
    signal_mask = UInt16(0),
    nav_message_index = 0,
)

"""
$(TYPEDEF)

The per-constellation masks of one HAS Mask block, as an immutable, inline
`AbstractVector{GalileoHASSatelliteMask}` of at most 15 entries (the ICD's
4-bit `Nsys`). Compares equal to a `Vector` of the same masks.

# Fields

$(TYPEDFIELDS)
"""
struct GalileoHASSatelliteMaskList <: AbstractVector{GalileoHASSatelliteMask}
    """
    The masks in broadcast order; slots past `length` hold a placeholder
    """
    items::NTuple{E6B_MAX_SYSTEMS,GalileoHASSatelliteMask}
    """
    Number of masks (`Nsys`)
    """
    length::Int
end

function GalileoHASSatelliteMaskList(masks::AbstractVector{GalileoHASSatelliteMask})
    length(masks) <= E6B_MAX_SYSTEMS ||
        throw(ArgumentError("a HAS mask covers at most $E6B_MAX_SYSTEMS constellations"))
    items = ntuple(
        i -> i <= length(masks) ? masks[i] : E6B_EMPTY_SATELLITE_MASK,
        Val(E6B_MAX_SYSTEMS),
    )
    GalileoHASSatelliteMaskList(items, length(masks))
end

Base.size(list::GalileoHASSatelliteMaskList) = (list.length,)
Base.IndexStyle(::Type{GalileoHASSatelliteMaskList}) = IndexLinear()
function Base.getindex(list::GalileoHASSatelliteMaskList, i::Int)
    @boundscheck checkbounds(list, i)
    return @inbounds list.items[i]
end

Base.convert(::Type{GalileoHASSatelliteMaskList}, list::GalileoHASSatelliteMaskList) = list
Base.convert(
    ::Type{GalileoHASSatelliteMaskList},
    masks::AbstractVector{GalileoHASSatelliteMask},
) = GalileoHASSatelliteMaskList(masks)

"""
    GalileoHASMask

One complete HAS Mask block: the per-constellation masks it defines, under the
Mask ID that later messages reference. An immutable, pointer-free value, so
storing it (in the Mask ID store, a message, the validated `data`) copies it.

# Fields

  - `mask_id::Int`: Mask ID this mask defines (0-31)
  - `satellite_masks::GalileoHASSatelliteMaskList`: one entry per corrected constellation, in broadcast order (an `AbstractVector{GalileoHASSatelliteMask}`)

# Reference

Galileo HAS SIS ICD, Issue 1.0, Table 15
"""
Base.@kwdef struct GalileoHASMask
    mask_id::Int
    satellite_masks::GalileoHASSatelliteMaskList
end

Base.:(==)(a::GalileoHASMask, b::GalileoHASMask) = fields_equal(a, b)

"""
    GalileoHASOrbitCorrection

Orbit correction for one satellite, in the satellite-centred NTW frame
(radial / in-track / cross-track, ICD §7.2).

A field is `nothing` where the ICD's "data not available" sentinel was
broadcast (the most negative two's-complement value: -10.24 m radial,
-16.384 m in- and cross-track). This is deliberately *not* folded to zero — a
zero correction and an absent correction are different facts, and GNSS-SDR's
choice to map both to 0 m downstream loses that.

# Fields

  - `GNSS_ID::Int`: GNSS index the satellite belongs to (Table 18)
  - `SVID::Int`: Satellite ID (Galileo SVID / GPS PRN)
  - `IOD_ref::Int`: Reference IOD of the corrected broadcast navigation data — IODnav for Galileo, IODE/IODC for GPS (Table 26)
  - `δ_radial::Union{Nothing,Float64}`: Delta Radial correction (meters, LSB 0.0025)
  - `δ_in_track::Union{Nothing,Float64}`: Delta In-Track correction (meters, LSB 0.008)
  - `δ_cross_track::Union{Nothing,Float64}`: Delta Cross-Track correction (meters, LSB 0.008)

# Reference

Galileo HAS SIS ICD, Issue 1.0, Tables 24-25
"""
Base.@kwdef struct GalileoHASOrbitCorrection
    GNSS_ID::Int
    SVID::Int
    IOD_ref::Int
    δ_radial::Union{Nothing,Float64} = nothing
    δ_in_track::Union{Nothing,Float64} = nothing
    δ_cross_track::Union{Nothing,Float64} = nothing
end

"""
    GalileoHASClockCorrection

Clock correction for one satellite (ICD §7.3).

`δ_clock` already has the constellation's Delta Clock Multiplier applied, so it
is the correction in meters ready to use; `multiplier` is reported alongside for
traceability. `δ_clock` is `nothing` for the "data not available" sentinel
(raw -4096); `do_not_use` marks the distinct "satellite shall not be used"
sentinel (raw +4095), which carries no correction either but means something
stronger.

# Fields

  - `GNSS_ID::Int`: GNSS index the satellite belongs to (Table 18)
  - `SVID::Int`: Satellite ID (Galileo SVID / GPS PRN)
  - `multiplier::Int`: Delta Clock Multiplier applied, 1-4 (Table 29)
  - `δ_clock::Union{Nothing,Float64}`: Delta clock correction (meters, LSB 0.0025 × `multiplier`)
  - `do_not_use::Bool`: The satellite shall not be used (Table 31)

# Reference

Galileo HAS SIS ICD, Issue 1.0, Tables 28-34
"""
Base.@kwdef struct GalileoHASClockCorrection
    GNSS_ID::Int
    SVID::Int
    multiplier::Int
    δ_clock::Union{Nothing,Float64} = nothing
    do_not_use::Bool = false
end

"""
    GalileoHASCodeBias

Code bias for one satellite/signal cell (ICD §7.4). `bias` is `nothing` for the
"data not available" sentinel (raw -1024).

# Fields

  - `GNSS_ID::Int`: GNSS index (Table 18)
  - `SVID::Int`: Satellite ID (Galileo SVID / GPS PRN)
  - `signal_index::Int`: Signal index within the constellation's Signal Mask (0-based, Table 20)
  - `bias::Union{Nothing,Float64}`: Code bias (meters, LSB 0.02)

# Reference

Galileo HAS SIS ICD, Issue 1.0, Tables 36-37
"""
Base.@kwdef struct GalileoHASCodeBias
    GNSS_ID::Int
    SVID::Int
    signal_index::Int
    bias::Union{Nothing,Float64} = nothing
end

"""
    GalileoHASPhaseBias

Phase bias for one satellite/signal cell (ICD §7.5). `bias` is `nothing` for the
"data not available" sentinel (raw -1024). `phase_discontinuity_indicator`
increments whenever the fixed ambiguity for this satellite and signal must be
re-initialised (ICD §5.2.6.1).

# Fields

  - `GNSS_ID::Int`: GNSS index (Table 18)
  - `SVID::Int`: Satellite ID (Galileo SVID / GPS PRN)
  - `signal_index::Int`: Signal index within the constellation's Signal Mask (0-based, Table 20)
  - `bias::Union{Nothing,Float64}`: Phase bias (cycles, LSB 0.01)
  - `phase_discontinuity_indicator::Int`: Phase Discontinuity Indicator (0-3)

# Reference

Galileo HAS SIS ICD, Issue 1.0, Tables 39-40
"""
Base.@kwdef struct GalileoHASPhaseBias
    GNSS_ID::Int
    SVID::Int
    signal_index::Int
    bias::Union{Nothing,Float64} = nothing
    phase_discontinuity_indicator::Int
end

"""
    e6b_block_capacity(T) -> Int

Most entries of type `T` one content block can carry. Every entry consumes at
least a fixed number of bits of the message — 45 for an orbit correction (an
8-bit GPS IOD plus 37 correction bits), 13 for a clock correction and for a
phase bias (11 + 2), 11 for a code bias — so the 13568-bit message budget bounds
the count, whatever the mask says.
"""
e6b_block_capacity(::Type{GalileoHASOrbitCorrection}) = cld(E6B_MAX_MESSAGE_BITS, 45)
e6b_block_capacity(::Type{GalileoHASClockCorrection}) = cld(E6B_MAX_MESSAGE_BITS, 13)
e6b_block_capacity(::Type{GalileoHASCodeBias}) = cld(E6B_MAX_MESSAGE_BITS, 11)
e6b_block_capacity(::Type{GalileoHASPhaseBias}) = cld(E6B_MAX_MESSAGE_BITS, 13)

"""
    GalileoHASCorrectionBlock{T}

One HAS MT1 content block: a validity interval, the message header context it
was broadcast under, and the per-satellite (or per-cell) corrections themselves.

Every block carries its own Validity Interval Index (ICD §5.2.2.1) starting at
the message's Time Of Hour, so blocks of one message can — and routinely do —
expire at different times. `mask_id` and `IOD_set_id` identify the satellite set
and the broadcast-ephemeris issue the corrections apply to (ICD §7.6).

!!! warning "Overwritten in place"

    A block is a preallocated, mutable buffer owned by the decoder state:
    [`decode!`](@ref) **overwrites** its fields and the contents of its
    `corrections` vector when a later message carries a block of the same kind.
    Copy it (`copy(state)`, or [`decode`](@ref)) to keep a snapshot.

# Fields

  - `TOH::Int`: Time Of Hour of the message that carried this block (seconds into the GST hour, 0-3599)
  - `mask_id::Int`: Mask ID the corrections are keyed to
  - `IOD_set_id::Int`: IOD Set ID the corrections are keyed to
  - `validity_interval::Union{Nothing,Int}`: Validity interval in seconds from `TOH` (`nothing` for the reserved index 15)
  - `corrections::Vector{T}`: the block's entries, in broadcast order — a buffer preallocated to the most entries a message can carry, resized within that capacity

# Reference

Galileo HAS SIS ICD, Issue 1.0, Tables 22, 27, 32, 35, 38
"""
mutable struct GalileoHASCorrectionBlock{T}
    TOH::Int
    mask_id::Int
    IOD_set_id::Int
    validity_interval::Union{Nothing,Int}
    const corrections::Vector{T}
end

function GalileoHASCorrectionBlock{T}(;
    TOH::Integer,
    mask_id::Integer,
    IOD_set_id::Integer,
    validity_interval::Union{Nothing,Integer},
    corrections::AbstractVector,
) where {T}
    buffer = sizehint!(T[], max(e6b_block_capacity(T), length(corrections)))
    append!(buffer, corrections)
    GalileoHASCorrectionBlock{T}(TOH, mask_id, IOD_set_id, validity_interval, buffer)
end

"""
An empty correction block of entry type `T` whose `corrections` buffer holds the
most entries one message can carry — the preallocated storage the decoder
overwrites.
"""
preallocated_block(::Type{T}) where {T} = GalileoHASCorrectionBlock{T}(;
    TOH = 0,
    mask_id = 0,
    IOD_set_id = 0,
    validity_interval = nothing,
    corrections = T[],
)

# The `corrections` vector (and mutability) make the default `===` too strict;
# see `fields_equal`.
Base.:(==)(a::GalileoHASCorrectionBlock, b::GalileoHASCorrectionBlock) = fields_equal(a, b)

"""
    overwrite!(dst::GalileoHASCorrectionBlock, src::GalileoHASCorrectionBlock) -> dst

Make `dst` an exact copy of `src` by **overwriting** its fields and the
contents of its `corrections` buffer (resized within its capacity).
"""
function overwrite!(
    dst::GalileoHASCorrectionBlock{T},
    src::GalileoHASCorrectionBlock{T},
) where {T}
    dst === src && return dst
    dst.TOH = src.TOH
    dst.mask_id = src.mask_id
    dst.IOD_set_id = src.IOD_set_id
    dst.validity_interval = src.validity_interval
    resize!(dst.corrections, length(src.corrections))
    copyto!(dst.corrections, src.corrections)
    return dst
end

"""
    GalileoHASMessage

One completely received and decoded HAS message — the atomic unit the ICD
defines, reassembled from `MS` HAS Encoded Pages by the RS erasure decoder.

Only Message Type 1 is specified (ICD Table 10), so in practice
`message_type == 1` and the six block fields are populated according to the
flags of the MT1 header; a flag that was 0 leaves its field `nothing`.

!!! warning "Overwritten in place"

    A message is a preallocated, mutable record owned by the decoder state:
    [`decode!`](@ref) **overwrites** it with the next completed message. In
    `state.data` its block fields are the very blocks `state.data` publishes as
    the latest of each kind (`data.message.clock_corrections === data.clock_corrections` when the message carried one). Copy the state to
    keep a snapshot.

# Fields

  - `message_id::Int`: Message ID (MID) the pages carried (0-31)
  - `message_type::Int`: Message Type (1 = satellite corrections)
  - `message_size::Int`: Message size in non-encoded pages (`MS`, 1-32)
  - `TOH::Int`: Time Of Hour (seconds into the GST hour, 0-3599)
  - `mask_id::Int`: Mask ID (0-31)
  - `IOD_set_id::Int`: IOD Set ID (0-31)
  - `mask::Union{Nothing,GalileoHASMask}`: Mask block, when the Mask Flag was set
  - `orbit_corrections::Union{Nothing,GalileoHASCorrectionBlock{GalileoHASOrbitCorrection}}`: Orbit Corrections block
  - `clock_corrections::Union{Nothing,GalileoHASCorrectionBlock{GalileoHASClockCorrection}}`: Clock Full-Set Corrections block
  - `clock_subset_corrections::Union{Nothing,GalileoHASCorrectionBlock{GalileoHASClockCorrection}}`: Clock Subset Corrections block
  - `code_biases::Union{Nothing,GalileoHASCorrectionBlock{GalileoHASCodeBias}}`: Code Biases block
  - `phase_biases::Union{Nothing,GalileoHASCorrectionBlock{GalileoHASPhaseBias}}`: Phase Biases block

# Reference

Galileo HAS SIS ICD, Issue 1.0, Tables 11-14
"""
Base.@kwdef mutable struct GalileoHASMessage
    message_id::Int
    message_type::Int
    message_size::Int
    TOH::Int
    mask_id::Int
    IOD_set_id::Int
    mask::Union{Nothing,GalileoHASMask} = nothing
    orbit_corrections::Union{Nothing,GalileoHASCorrectionBlock{GalileoHASOrbitCorrection}} =
        nothing
    clock_corrections::Union{Nothing,GalileoHASCorrectionBlock{GalileoHASClockCorrection}} =
        nothing
    clock_subset_corrections::Union{
        Nothing,
        GalileoHASCorrectionBlock{GalileoHASClockCorrection},
    } = nothing
    code_biases::Union{Nothing,GalileoHASCorrectionBlock{GalileoHASCodeBias}} = nothing
    phase_biases::Union{Nothing,GalileoHASCorrectionBlock{GalileoHASPhaseBias}} = nothing
end

"""
An empty message record, the preallocated storage the decoder overwrites.
"""
preallocated_message() = GalileoHASMessage(;
    message_id = 0,
    message_type = 0,
    message_size = 0,
    TOH = 0,
    mask_id = 0,
    IOD_set_id = 0,
)

# Mutable, and the nested blocks carry vectors: compare by value.
Base.:(==)(a::GalileoHASMessage, b::GalileoHASMessage) = fields_equal(a, b)

"""
    e6b_overwrite_message!(dst, src, orbit, clock, clock_subset, code, phase) -> dst

**Overwrite** message record `dst` with `src`'s header and mask, pointing each
block field at the given block where `src` carried that kind of block and at
`nothing` where it did not. The blocks are the caller's (already overwritten)
copies of `src`'s, so `dst` never references a block `src` owns.
"""
function e6b_overwrite_message!(
    dst::GalileoHASMessage,
    src::GalileoHASMessage,
    orbit,
    clock,
    clock_subset,
    code,
    phase,
)
    dst.message_id = src.message_id
    dst.message_type = src.message_type
    dst.message_size = src.message_size
    dst.TOH = src.TOH
    dst.mask_id = src.mask_id
    dst.IOD_set_id = src.IOD_set_id
    dst.mask = src.mask
    dst.orbit_corrections = isnothing(src.orbit_corrections) ? nothing : orbit
    dst.clock_corrections = isnothing(src.clock_corrections) ? nothing : clock
    dst.clock_subset_corrections =
        isnothing(src.clock_subset_corrections) ? nothing : clock_subset
    dst.code_biases = isnothing(src.code_biases) ? nothing : code
    dst.phase_biases = isnothing(src.phase_biases) ? nothing : phase
    return dst
end

"""
$(TYPEDEF)

One preallocated block of each MT1 content kind — the buffers a message's
blocks are parsed into before they are copied out.

# Fields

$(TYPEDFIELDS)
"""
struct GalileoHASBlockBuffers
    orbit_corrections::GalileoHASCorrectionBlock{GalileoHASOrbitCorrection}
    clock_corrections::GalileoHASCorrectionBlock{GalileoHASClockCorrection}
    clock_subset_corrections::GalileoHASCorrectionBlock{GalileoHASClockCorrection}
    code_biases::GalileoHASCorrectionBlock{GalileoHASCodeBias}
    phase_biases::GalileoHASCorrectionBlock{GalileoHASPhaseBias}
end

GalileoHASBlockBuffers() = GalileoHASBlockBuffers(
    preallocated_block(GalileoHASOrbitCorrection),
    preallocated_block(GalileoHASClockCorrection),
    preallocated_block(GalileoHASClockCorrection),
    preallocated_block(GalileoHASCodeBias),
    preallocated_block(GalileoHASPhaseBias),
)

# ---- Decoded data container --------------------------------------------------

"""
    GalileoE6BData <: AbstractGalileoData

Decoded Galileo E6-B C/NAV (High Accuracy Service) data.

Two views of the same stream are kept, because both are genuinely useful:

  - `message` is the most recently *completed* HAS message, exactly as the ICD
    defines it — one atomic unit, with whichever content blocks its header
    flagged.
  - `masks` and the five correction-block fields accumulate the latest of each
    kind across messages. This is what a correction consumer wants: HAS
    routinely splits a mask + orbit + bias message from a clock-only message
    (ICD §5.1), so no single message holds a usable set.

Every field here has passed the per-page CRC-24Q and the RS erasure decode, so
`raw_data` and `data` track each other (there is no cross-message issue-of-data
vote to run — see `validate_data`).

# Storage

Every container here — the message record, the Mask ID store and the five
correction blocks — is preallocated by the decoder state (to the ICD's
maxima) and **overwritten in place** by [`decode!`](@ref); `data` has its own
set, which `validate_data` overwrites with a copy of `raw_data`'s, so the two
never share a container.

# Service status

  - `HAS_status::HASStatus`: HAS status from the most recent valid page (Table 9)

# Decoded content

  - `message::GalileoHASMessage`: the most recently completed message in full,
    including its header — `message.TOH`, `.mask_id`, `.IOD_set_id` and
    `.message_id`. Those are deliberately *not* mirrored as flat fields here: the
    blocks below can come from *different* messages, so each carries its own
    `TOH` / `mask_id` / `IOD_set_id`, and those are the ones to age a correction
    against. A flat `data.TOH` alongside `data.orbit_corrections.TOH` would look
    like the same fact and is not.
  - `masks::SlotDictionary{GalileoHASMask,32}`: every mask received so far, keyed by Mask ID (0-31)
  - `orbit_corrections::GalileoHASCorrectionBlock{GalileoHASOrbitCorrection}`: latest Orbit Corrections block
  - `clock_corrections::GalileoHASCorrectionBlock{GalileoHASClockCorrection}`: latest Clock Full-Set Corrections block
  - `clock_subset_corrections::GalileoHASCorrectionBlock{GalileoHASClockCorrection}`: latest Clock Subset Corrections block
  - `code_biases::GalileoHASCorrectionBlock{GalileoHASCodeBias}`: latest Code Biases block
  - `phase_biases::GalileoHASCorrectionBlock{GalileoHASPhaseBias}`: latest Phase Biases block

# Reference

Galileo HAS SIS ICD, Issue 1.0
"""
Base.@kwdef struct GalileoE6BData <: AbstractGalileoData
    HAS_status::Union{Nothing,HASStatus} = nothing
    message::Union{Nothing,GalileoHASMessage} = nothing
    masks::Union{Nothing,SlotDictionary{GalileoHASMask,E6B_NUM_MASK_IDS}} = nothing
    orbit_corrections::Union{Nothing,GalileoHASCorrectionBlock{GalileoHASOrbitCorrection}} =
        nothing
    clock_corrections::Union{Nothing,GalileoHASCorrectionBlock{GalileoHASClockCorrection}} =
        nothing
    clock_subset_corrections::Union{
        Nothing,
        GalileoHASCorrectionBlock{GalileoHASClockCorrection},
    } = nothing
    code_biases::Union{Nothing,GalileoHASCorrectionBlock{GalileoHASCodeBias}} = nothing
    phase_biases::Union{Nothing,GalileoHASCorrectionBlock{GalileoHASPhaseBias}} = nothing
end

@inline function GalileoE6BData(
    data::GalileoE6BData;
    HAS_status = data.HAS_status,
    message = data.message,
    masks = data.masks,
    orbit_corrections = data.orbit_corrections,
    clock_corrections = data.clock_corrections,
    clock_subset_corrections = data.clock_subset_corrections,
    code_biases = data.code_biases,
    phase_biases = data.phase_biases,
)
    GalileoE6BData(
        HAS_status,
        message,
        masks,
        orbit_corrections,
        clock_corrections,
        clock_subset_corrections,
        code_biases,
        phase_biases,
    )
end

# The mutable `masks` store, message record and correction blocks make the
# default struct `==` (which falls back to `===`) too strict.
Base.:(==)(a::GalileoE6BData, b::GalileoE6BData) = fields_equal(a, b)

# Every container field at its ICD size: the message record, one mask slot per
# Mask ID (0-31), and each correction block's buffer at the most entries one
# message can carry (`e6b_block_capacity`).
function preallocated_data(::Type{GalileoE6BData})
    blocks = GalileoHASBlockBuffers()
    GalileoE6BData(;
        message = preallocated_message(),
        masks = SlotDictionary{GalileoHASMask,E6B_NUM_MASK_IDS}(),
        blocks.orbit_corrections,
        blocks.clock_corrections,
        blocks.clock_subset_corrections,
        blocks.code_biases,
        blocks.phase_biases,
    )
end

"""
$(TYPEDSIGNATURES)

Always `false` for Galileo E6-B.

C/NAV carries no ephemeris, clock polynomial or week number of its own: the HAS
message is a set of *corrections* to the broadcast navigation data of another
signal (Galileo I/NAV or GPS LNAV, selected per constellation by the mask's
Navigation Message Index). A positioning engine therefore pairs an E6-B decoder
with an I/NAV or LNAV decoder rather than using it alone, so this readiness gate
— which asks whether *this* satellite's own positioning set is complete — can
never be satisfied here.

Use `state.data.orbit_corrections`, `.clock_corrections`, `.code_biases` and
`.phase_biases` (each with its own validity interval and `IOD_set_id`) together
with the corresponding ephemeris decoder instead.
"""
is_decoding_completed_for_positioning(data::GalileoE6BData) = false

# C/NAV stamps each HAS correction block with a time *of hour*, not a time of
# week, and broadcasts no week number to place it in — see the type docstring on
# why those stay on the blocks that own them rather than being mirrored flat.
get_time_of_week(::GalileoE6BData) = nothing

# HAS carries orbit, clock and bias corrections; inter-system time offsets are
# not among them. The `GNSS_ID` fields in the satellite masks are constellation
# indices, not time-offset identifiers.
get_time_offset(::GNSSDecoderState{<:GalileoE6BData}, ::TimeSystem) = nothing

# ---- Page store (ICD §6.4) ---------------------------------------------------

"""
$(TYPEDEF)

Encoded pages collected so far for one Message ID, awaiting the `message_size`
distinct pages the Reed-Solomon erasure decoder needs (ICD §6.4).

Mutable and mutated in place inside the decoder cache: it is exactly the kind of
"still partial" state `CONTEXT.md` says belongs there. The page store holds one
preallocated group per Message ID, sized for the largest message (32 pages),
and a new message under that ID **overwrites** it (see `e6b_reopen_group!`).

# Fields

$(TYPEDFIELDS)
"""
mutable struct GalileoHASPageGroup
    """
    Message Type of the collected pages (pages of a different type never mix)
    """
    message_type::Int
    """
    Message size `MS` in non-encoded pages — the `k` of the RS decode
    """
    message_size::Int
    """
    Value of the cache's page counter when this group was opened, for the 150 s timeout
    """
    opened_at::Int
    """
    HAS Page IDs of the collected pages, in arrival order — its length is how many
    are held, and reaching `message_size` is what completes the group (capacity 32)
    """
    const page_ids::Vector{Int}
    """
    Collected pages, `32 × 53` octets; row `i` belongs to `page_ids[i]`, so rows
    beyond `length(page_ids)` are not filled (or hold an earlier message's pages)
    """
    const octets::Matrix{UInt8}
end

GalileoHASPageGroup(message_type::Int, message_size::Int, opened_at::Int) =
    GalileoHASPageGroup(
        message_type,
        message_size,
        opened_at,
        sizehint!(Int[], E6B_MAX_MESSAGE_PAGES),
        zeros(UInt8, E6B_MAX_MESSAGE_PAGES, E6B_OCTETS_PER_PAGE),
    )

"""
**Overwrite** `group` so it is an empty group for a new message: the given
header fields, no pages collected.
"""
function e6b_reopen_group!(
    group::GalileoHASPageGroup,
    message_type::Int,
    message_size::Int,
    opened_at::Int,
)
    group.message_type = message_type
    group.message_size = message_size
    group.opened_at = opened_at
    empty!(group.page_ids)
    return group
end

# Being *mutable* makes the default `==` reference equality outright — not just
# for the `Vector`/`Matrix` fields — so without this two decoders fed the same
# stream would compare unequal for as long as either holds a partial message,
# which is most of the time. `GalileoE6BCache` compares its page store by value,
# so this is what makes that comparison mean anything. Only the filled rows of
# `octets` are compared: the rest is leftover from earlier messages.
function Base.:(==)(a::GalileoHASPageGroup, b::GalileoHASPageGroup)
    a.message_type == b.message_type &&
    a.message_size == b.message_size &&
    a.opened_at == b.opened_at &&
    a.page_ids == b.page_ids || return false
    for j = 1:E6B_OCTETS_PER_PAGE, i = 1:length(a.page_ids)
        a.octets[i, j] == b.octets[i, j] || return false
    end
    return true
end

"""
The page store: one slot per Message ID (0-31), each holding its own
preallocated [`GalileoHASPageGroup`](@ref) whether or not the slot is occupied.
"""
const GalileoHASPageStore = SlotDictionary{GalileoHASPageGroup,E6B_NUM_MESSAGE_IDS}

function GalileoHASPageStore()
    store = SlotDictionary{GalileoHASPageGroup,E6B_NUM_MESSAGE_IDS}(
        [GalileoHASPageGroup(0, 0, 0) for _ = 1:E6B_NUM_MESSAGE_IDS],
        SlotIndices(zeros(Bool, E6B_NUM_MESSAGE_IDS), 0),
    )
    return store
end

# The groups are mutable and overwritten in place, so a copy of the store (the
# one `copy(state)` makes via `duplicate`) must get groups of its own, where the
# generic `SlotDictionary` copy would share them.
Base.copy(store::GalileoHASPageStore) = GalileoHASPageStore(
    [e6b_copy_group(group) for group in store.values],
    SlotIndices(copy(store.indices.occupied), length(store)),
)

function e6b_copy_group(group::GalileoHASPageGroup)
    copied = GalileoHASPageGroup(group.message_type, group.message_size, group.opened_at)
    append!(copied.page_ids, group.page_ids)
    copyto!(copied.octets, group.octets)
    return copied
end

"""
$(TYPEDEF)

A reassembled HAS message held back because its body references a Mask ID that
has not been received yet (ICD §5.1.1.1) — every block's length derives from the
mask, so nothing past the header can be parsed without it.

One slot is enough: HAS broadcasts a defining mask every few messages, so the
orphan worth keeping is the newest. It carries the same `opened_at` page stamp as
[`GalileoHASPageGroup`](@ref) and expires on the same clock, so a mask arriving
long afterwards cannot resurrect corrections whose validity intervals have run
out. The cache preallocates one and **overwrites** it with each newly held
message.

# Fields

$(TYPEDFIELDS)
"""
mutable struct GalileoHASPendingMessage
    """
    Message ID the message's pages carried
    """
    message_id::Int
    """
    Message Type of the reassembled message
    """
    message_type::Int
    """
    Message size `MS` in non-encoded pages
    """
    message_size::Int
    """
    Value of the cache's page counter when the message was held, for the timeout
    """
    opened_at::Int
    """
    The reassembled message octets, awaiting a mask (capacity 32 × 53)
    """
    const octets::Vector{UInt8}
end

GalileoHASPendingMessage() =
    GalileoHASPendingMessage(0, 0, 0, 0, sizehint!(UInt8[], E6B_MAX_MESSAGE_OCTETS))

# `octets` is a Vector, so the default struct `==` would be reference equality.
Base.:(==)(a::GalileoHASPendingMessage, b::GalileoHASPendingMessage) = fields_equal(a, b)

# ---- HAS message bit reader --------------------------------------------------
#
# The MT1 body is a variable-length, unaligned bit stream up to 32 × 424 =
# 13568 bits long, and every block's length depends on values read earlier in
# it, so it cannot be indexed by precomputed offsets the way a fixed-layout word
# can. A sequential MSB-first reader over the reassembled octets is the honest
# representation; it also makes "ran off the end of the message" — which happens
# whenever a mask is stale or a future block type is appended (ICD §5.1 warns
# forward-compatible receivers to expect exactly that) — a single check.

"""
$(TYPEDEF)

Sequential MSB-first bit reader over the octets of a reassembled HAS message.
The decoder keeps one in its cache and rewinds it onto each message
(`e6b_rewind!`), overwriting its fields.

# Fields

$(TYPEDFIELDS)
"""
mutable struct HASBitReader
    """
    Buffer holding the reassembled message octets (possibly longer than the message)
    """
    octets::Vector{UInt8}
    """
    Number of leading octets of `octets` that make up the message
    """
    num_octets::Int
    """
    Number of bits consumed so far
    """
    position::Int
end

HASBitReader(octets::Vector{UInt8}) = HASBitReader(octets, length(octets), 0)

"""
Point `reader` at the first `num_octets` octets of `octets`, from the start —
**overwrites** the reader's fields.
"""
function e6b_rewind!(reader::HASBitReader, octets::Vector{UInt8}, num_octets::Int)
    num_octets <= length(octets) ||
        throw(DimensionMismatch("message longer than its octet buffer"))
    reader.octets = octets
    reader.num_octets = num_octets
    reader.position = 0
    return reader
end

"""
Bits left unread in the message.
"""
bits_remaining(reader::HASBitReader) = 8 * reader.num_octets - reader.position

"""
    peek_bits(reader, position, num_bits) -> UInt64

The `num_bits` (at most 64) bits starting `position` bits into the message,
MSB-first, without consuming anything.
"""
function peek_bits(reader::HASBitReader, position::Int, num_bits::Int)
    value = UInt64(0)
    @inbounds for p = position:(position+num_bits-1)
        byte = reader.octets[(p>>3)+1]
        value = (value << 1) | UInt64((byte >> (7 - (p & 7))) & 0x01)
    end
    return value
end

"""
    read_bits!(reader, num_bits) -> UInt64

Consume the next `num_bits` (at most 64) MSB-first. Callers must have checked
[`bits_remaining`](@ref) first.
"""
function read_bits!(reader::HASBitReader, num_bits::Int)
    value = peek_bits(reader, reader.position, num_bits)
    reader.position += num_bits
    return value
end

"""
    read_signed_bits!(reader, num_bits) -> Int64

Consume `num_bits` as a two's-complement integer, sign bit in the MSB (the ICD's
convention for every correction field).
"""
function read_signed_bits!(reader::HASBitReader, num_bits::Int)
    # `get_twos_complement_num` is the package's single sign-extension rule
    # (`src/bit_fiddling.jl`), including its 32-bit-safe routing through UInt64;
    # the whole field is the value, so it is read at offset 1 of its own width.
    get_twos_complement_num(read_bits!(reader, num_bits), num_bits, 1, num_bits)
end

"""
$(TYPEDEF)

Everything [`parse_has_message!`](@ref) writes into, preallocated once: the bit
reader, a scratch list for the constellations of a Mask block, one block of
each content kind, and the message record the parse returns (whose block fields
point into `blocks`). **Overwritten** by every parse, so the returned message is
valid only until the next one — the decoder copies it into `raw_data` straight
away.

# Fields

$(TYPEDFIELDS)
"""
struct GalileoHASParser
    """
    Bit reader, rewound onto each message
    """
    reader::HASBitReader
    """
    Per-constellation masks of the Mask block being parsed (capacity 15)
    """
    satellite_masks::Vector{GalileoHASSatelliteMask}
    """
    The content blocks being parsed
    """
    blocks::GalileoHASBlockBuffers
    """
    The parsed message
    """
    message::GalileoHASMessage
end

GalileoHASParser() = GalileoHASParser(
    HASBitReader(UInt8[]),
    sizehint!(GalileoHASSatelliteMask[], E6B_MAX_SYSTEMS),
    GalileoHASBlockBuffers(),
    preallocated_message(),
)

# ---- Cache -------------------------------------------------------------------

"""
$(TYPEDEF)

Per-decoder cache for Galileo E6-B.

Beyond the shared soft-symbol deque this holds the page-level FEC scratch, the
long-lived Viterbi decoder, the in-flight HAS page store — the partial-message
state that is the whole point of the HPVRS outer layer — and the preallocated
buffers the reassembly, the parse and the decoded data are written into (the RS
generator matrix is a shared constant, `E6B_GENERATOR_MATRIX`).

Everything here is **overwritten in place** by [`decode!`](@ref): the page
store, the page counter, the held-back message, and every scratch buffer (the
`CONTEXT.md` rule: decoded fields in `data`/`raw_data`, genuinely-in-flight
state in the cache).

# Fields

$(TYPEDFIELDS)
"""
struct GalileoE6BCache <: AbstractGNSSCache
    """
    Soft-symbol buffer (1016 = 1000 page + 16 next-page sync pattern)
    """
    soft_buffer::CircularDeque{Float32}
    """
    Polarity-resolved 984-symbol FEC window copied out per sync attempt
    """
    fec_window::Vector{Float32}
    """
    Viterbi decoder and its scratch buffers (K = 486, N = 984), built once and reused across pages
    """
    viterbi::GalileoViterbiScratch
    """
    Encoded pages collected per Message ID (0-31), awaiting completion; one preallocated group per slot
    """
    page_groups::GalileoHASPageStore
    """
    Count of accepted C/NAV pages — one per second, the clock for the ICD's 150 s message timeout
    """
    page_counter::Base.RefValue{Int}
    """
    A completed message whose body needs a Mask ID not yet received, held for re-parsing
    (`nothing`, or the preallocated `held_message`)
    """
    pending_message::Base.RefValue{Union{Nothing,GalileoHASPendingMessage}}
    """
    Preallocated record `pending_message` points at while a message is held
    """
    held_message::GalileoHASPendingMessage
    """
    Working matrices of the RS erasure decode (32 × 32)
    """
    rs_scratch::RSErasureScratch
    """
    Reassembled message octets (32 × 53), overwritten by each reassembly
    """
    message_octets::Vector{UInt8}
    """
    Reader, scratch blocks and message record the MT1 parse writes into
    """
    parser::GalileoHASParser
    """
    Preallocated containers `raw_data` and `data` are decoded into
    """
    storage::DataStorage{GalileoE6BData}
end

GalileoE6BCache() = GalileoE6BCache(
    CircularDeque{Float32}(E6B_WINDOW_SYMBOLS),
    Vector{Float32}(undef, E6B_ENCODED_SYMBOLS),
    GalileoViterbiScratch(E6B_PAGE_BITS, E6B_ENCODED_SYMBOLS),
    GalileoHASPageStore(),
    Ref(0),
    Ref{Union{Nothing,GalileoHASPendingMessage}}(nothing),
    GalileoHASPendingMessage(),
    RSErasureScratch(E6B_RS_CODE_DIMENSION),
    zeros(UInt8, E6B_MAX_MESSAGE_OCTETS),
    GalileoHASParser(),
    DataStorage{GalileoE6BData}(),
)

# The Viterbi handle, the FEC, RS and parse scratch and the preallocated data
# storage are derived, not state, so they are excluded — but the page store, the page counter and the held orphan are
# genuine in-flight state and are compared, the way `GalileoINAVCache` compares
# its stitched even page and almanac chain and `BeiDouDNAVCache` its pending
# pages. (`BeiDouB2bCache` compares only its deque because it holds nothing
# else.) Without this, a cache holding 31 of 32 pages plus an orphan would
# compare equal to an empty one.
function Base.:(==)(a::GalileoE6BCache, b::GalileoE6BCache)
    deques_equal(a.soft_buffer, b.soft_buffer) &&
        a.page_groups == b.page_groups &&
        a.page_counter[] == b.page_counter[] &&
        a.pending_message[] == b.pending_message[]
end

# ---- Decoder state -----------------------------------------------------------

"""
$(TYPEDSIGNATURES)

Create a decoder state for Galileo E6-B C/NAV (High Accuracy Service) messages.

Initializes a [`GNSSDecoderState`](@ref) configured for decoding the Galileo HAS
message from the FEC-encoded 1000 sps soft symbols of the E6-B component. Each
sync attempt matches the 16-symbol sync pattern `1011011101110000` at both ends
of the 1016-symbol window, 123×8 deinterleaves and Viterbi-decodes the 984
encoded symbols to a 486-bit page, and gates it on CRC-24Q. Valid, non-dummy
pages are accumulated per Message ID until `MS` distinct HAS Page IDs are held,
at which point the Reed-Solomon erasure decoder recovers the HAS message and its
Message Type 1 content blocks are parsed into a [`GalileoE6BData`](@ref).

!!! note "One satellite is slow; several are fast"

    A HAS message needs `MS` (up to 32) *distinct* encoded pages. One satellite
    broadcasts one page per second, so a single-satellite decoder needs up to 32
    seconds per message. HAS is designed for pages to be pooled across
    satellites — a receiver tracking several E6-B satellites completes messages
    far sooner. Each decoder state here accumulates only its own satellite's
    pages; combining them across satellites is a receiver-level concern.

# Arguments

  - `prn::Int`: Pseudo-Random Noise code identifier (1-50 for Galileo E6)

# Returns

  - `GNSSDecoderState{GalileoE6BData}`: Initialized decoder state for Galileo E6-B

# Example

```julia
state = GalileoE6BDecoderState(1)  # Create decoder for PRN 1
state = decode(state, soft_symbols, num_symbols)
if !isnothing(state.data.orbit_corrections)
    # Apply HAS corrections to the I/NAV ephemeris of the masked satellites
end
```

# See Also

  - [`GNSSDecoderState`](@ref): The underlying state structure
  - [`GalileoE1BDecoderState`](@ref): The I/NAV decoder whose ephemeris HAS corrects
  - [`decode`](@ref): Decode soft symbols using this state
  - [`reset_decoder_state`](@ref): Reset after signal loss
"""
function GalileoE6BDecoderState(prn)
    GNSSDecoderState(
        prn,
        GalileoE6BData(),
        GalileoE6BData(),
        GalileoE6BConstants(),
        GalileoE6BCache(),
        nothing,
        false,
    )
end

# Dispatch from a GNSSSignals system type. C/NAV rides on the E6-B (data)
# component — `GalileoE6B` — while `GalileoE6C` is the dataless pilot, so only
# `GalileoE6B` maps to a decoder.
function GNSSDecoderState(system::GalileoE6B, prn)
    GalileoE6BDecoderState(prn)
end

# The signal this decoder demodulates. Signal metadata is forwarded through it
# (see `src/gps/l1ca.jl`), so the state reports the E6 band and 1000 Hz symbol
# rate.
get_signal_type(::GalileoE6BConstants) = GalileoE6B

"""
$(TYPEDSIGNATURES)

Reset the Galileo E6-B decoder state after a signal loss or reacquisition.

Clears the soft-symbol buffer and the in-flight HAS page store, and drops the
published correction blocks, while keeping the received masks in `raw_data`:
masks change only when the corrected satellite/signal set changes, so a mask
survives an outage and lets a clock-only message be parsed immediately after
reacquisition.

The page store *is* cleared, deliberately. The ICD's 150 s message timeout is
counted in received pages (one page = one second of signal), so a decoder that
went dark for ten minutes would otherwise resume with pages that look fresh but
are not.

In place: the soft-symbol buffer and the page store are emptied (overwritten)
and the held message is dropped; the preallocated containers stay in the cache,
to be overwritten by the next decode. Allocates nothing.

# Arguments

  - `state::GNSSDecoderState{<:GalileoE6BData}`: Current Galileo E6-B decoder state

# Returns

  - `GNSSDecoderState{<:GalileoE6BData}`: Reset decoder state with cleared buffers

# See Also

  - [`GalileoE6BDecoderState`](@ref): Create a fresh decoder state
  - [`decode`](@ref): Continue decoding after reset
"""
function reset_decoder_state!(state::GNSSDecoderState{<:GalileoE6BData})
    empty!(state.cache.soft_buffer)
    empty!(state.cache.page_groups)
    state.cache.pending_message[] = nothing
    GNSSDecoderState(
        state;
        raw_data = GalileoE6BData(;
            HAS_status = state.raw_data.HAS_status,
            masks = state.raw_data.masks,
        ),
        data = GalileoE6BData(),
        num_bits_after_valid_syncro_sequence = nothing,
    )
end

# No `packed_buffer_type` method: E6-B overrides `try_sync` and reads the 32
# sync bits straight from the soft buffer, so the 1016-bit packed window the
# default would build (once per symbol, at 1000 sps) is never needed.

# ---- Page FEC and sync -------------------------------------------------------

"""
$(TYPEDEF)

Result of a successful C/NAV page sync: the CRC-validated 486-bit page and the
resolved polarity. Produced by `try_sync`, threaded through
`complement_buffer_if_necessary` to `decode_syncro_sequence`.

# Fields

$(TYPEDFIELDS)
"""
struct GalileoE6BSync
    """
    CRC-validated page, packed MSB-first (the first reserved bit at bit 1)
    """
    page::UInt512
    """
    Whether the symbol stream is 180-degrees phase shifted
    """
    polarity_flipped::Bool
end

"""
    try_sync(state::GNSSDecoderState{<:GalileoE6BData}) -> Union{Nothing,GalileoE6BSync}

C/NAV page sync (HAS SIS ICD, Issue 1.0, §2.3). Two gates, cheapest first:

 1. The 16-bit sync pattern `1011011101110000` must appear at both ends of the
    1016-symbol window (start of this page and start of the next), both upright
    or both inverted (the shared `find_preamble` rule, resolving the 180-degree
    carrier ambiguity).
 2. The 984 encoded symbols must Viterbi-decode to a 486-bit page whose
    CRC-24Q over the leading 462 bits matches the broadcast checksum.

Only then is the page handed to `decode_syncro_sequence`, so a corrupted page
can never enter the HAS page store — which matters more here than elsewhere,
because the RS erasure decoder trusts every collected page absolutely: one bad
octet corrupts the whole reassembled message, with no syndrome left to catch it
(ICD §6.4 models the channel as a *binary erasure* channel precisely because
the page CRC removes errors first).
"""
function try_sync(state::GNSSDecoderState{<:GalileoE6BData})
    polarity_flipped = find_preamble_in_deque(
        soft_buffer(state),
        state.constants.preamble,
        E6B_SYNC_SYMBOLS,
        E6B_PAGE_SYMBOLS,
    )
    isnothing(polarity_flipped) && return nothing
    window = copy_soft_window!(
        state.cache.fec_window,
        soft_buffer(state),
        E6B_SYNC_SYMBOLS,
        E6B_ENCODED_SYMBOLS,
        polarity_flipped,
    )
    page = galileo_viterbi(state.cache.viterbi, window, E6B_INTERLEAVER_COLUMNS, UInt512)
    # CRC-24Q over the whole 486-bit page (message ++ checksum) must be zero.
    # 486 bits is not a whole number of octets; the leading two zero bits of the
    # 61-octet representation are neutral because CRC-24Q initialises its
    # register to zero. GNSS-SDR and PocketSDR pad the same way.
    crc24q(page, E6B_PAGE_OCTETS) == 0 || return nothing
    return GalileoE6BSync(page, polarity_flipped)
end

"""
    complement_buffer_if_necessary(state::GNSSDecoderState{<:GalileoE6BData}, sync)

Record the polarity resolved by `try_sync` on the state and pass the
[`GalileoE6BSync`](@ref) through unchanged (its page bits are already
polarity-resolved and CRC-validated).
"""
function complement_buffer_if_necessary(
    state::GNSSDecoderState{<:GalileoE6BData},
    sync::GalileoE6BSync,
)
    GNSSDecoderState(state; is_shifted_by_180_degrees = sync.polarity_flipped), sync
end

# ---- MT1 parsing (ICD §5) ----------------------------------------------------
#
# Every parser below writes into preallocated storage (a `GalileoHASParser`):
# a block parser **overwrites** the block it is handed and returns it, or
# returns `nothing` when the message runs out of bits — in which case the
# block's contents are unspecified and it is not published.

"""
$(TYPEDEF)

Everything a content-block parser needs beyond the bit reader: the mask its
lengths derive from, and the three MT1 header fields every block it produces is
stamped with.

All five content blocks of an MT1 body take exactly this, so it travels as one
value rather than four repeated parameters, and
[`e6b_start_block!`](@ref) stamps a block with it. The mask belongs here for the
same reason the header fields do: *every* block's length is derived from it
(HAS SIS ICD, Issue 1.0, §5.2), so no parser can run without one.

# Fields

$(TYPEDFIELDS)
"""
struct GalileoHASBlockContext
    """
    The mask the corrected satellite/signal set and every block length come from
    """
    mask::GalileoHASMask
    """
    Time Of Hour of the carrying message (seconds into the GST hour, 0-3599)
    """
    TOH::Int
    """
    Mask ID the corrections are keyed to
    """
    mask_id::Int
    """
    IOD Set ID the corrections are keyed to
    """
    IOD_set_id::Int
end

"""
    e6b_start_block!(block, context, validity_interval) -> block

**Overwrite** `block` with the header context it is broadcast under and an
empty `corrections` list (its capacity kept), ready for a parser to push the
block's entries.
"""
function e6b_start_block!(
    block::GalileoHASCorrectionBlock,
    context::GalileoHASBlockContext,
    validity_interval::Union{Nothing,Int},
)
    block.TOH = context.TOH
    block.mask_id = context.mask_id
    block.IOD_set_id = context.IOD_set_id
    block.validity_interval = validity_interval
    empty!(block.corrections)
    return block
end

"""
    e6b_expand_mask(mask, width; first_index = 0) -> GalileoHASMaskIndices

The indices a `width`-bit HAS mask's set bits name, MSB = index 0, counting from
`first_index` — as a [`GalileoHASMaskIndices`](@ref), which computes them from
the mask instead of storing them.

Both masks the ICD defines are this walk and differ only in those two numbers:
the 40-bit Satellite Mask names Galileo SVID / GPS PRN `index + 1` (Table 19),
and the 16-bit Signal Mask names 0-based signal indices (Table 20).
"""
e6b_expand_mask(mask::Unsigned, width::Int; first_index::Int = 0) =
    GalileoHASMaskIndices(UInt64(mask), width, first_index)

"""
    parse_has_mask_block!(parser, mask_id) -> Union{Nothing,GalileoHASMask}

Parse the Mask block (ICD §5.2.1) from `parser.reader`, using
`parser.satellite_masks` as scratch (overwritten). Returns `nothing` if the
message runs out of bits, declares the reserved `Nsys` value 0, or names a
reserved GNSS ID — in every case every following block's length is unknown, so
the body parse must stop.

Note the 6 reserved bits after the per-constellation masks are consumed
unconditionally, whereas GNSS-SDR consumes them (and the `Nsys` field itself)
only when `Nsys != 0`. That is not a framing bug on its part — it clears its
`have_mask` flag and gates every later block on it, so it bails out rather than
desynchronising — but this parser rejects `Nsys == 0` outright instead, which is
what the ICD asks for.
"""
function parse_has_mask_block!(parser::GalileoHASParser, mask_id::Int)
    reader = parser.reader
    bits_remaining(reader) >= 4 || return nothing
    num_systems = Int(read_bits!(reader, 4))
    # `Nsys` takes "values from 1 to 15 (value "0" is Reserved)" (ICD §5.2.1). A
    # zero would otherwise yield a mask correcting nothing, which is worse than
    # useless: it would be cached under a real Mask ID and then *satisfy* later
    # `Mask Flag = 0` messages, so they would parse to empty correction blocks
    # instead of being held until the genuine mask arrives.
    num_systems == 0 && return nothing
    masks = empty!(parser.satellite_masks)
    for _ = 1:num_systems
        bits_remaining(reader) >= 4 + 40 + 16 + 1 || return nothing
        GNSS_ID = Int(read_bits!(reader, 4))
        # A reserved GNSS ID makes the Reference IOD width, and therefore every
        # later block length, unknowable (ICD Table 26).
        isnothing(e6b_iod_ref_length(GNSS_ID)) && return nothing
        satellite_mask = read_bits!(reader, 40)
        signal_mask = UInt16(read_bits!(reader, 16))
        cell_mask_available = read_bits!(reader, 1) == 1
        SVIDs = e6b_expand_mask(satellite_mask, 40; first_index = 1)
        signal_indices = e6b_expand_mask(signal_mask, 16)
        cell_mask = nothing
        if cell_mask_available
            # L_CM = Nsig · Nsat, read satellite-major: the table is "read from
            # left to right and from top to bottom" of Nsat rows by Nsig columns
            # (ICD §5.2.1.5, Eq. 3). Each row is one Nsig-bit field.
            num_satellites = length(SVIDs)
            num_signals = length(signal_indices)
            bits_remaining(reader) >= num_satellites * num_signals + 3 || return nothing
            start = reader.position
            rows = ntuple(Val(E6B_MAX_SATELLITES)) do row
                row <= num_satellites ?
                UInt16(peek_bits(reader, start + (row - 1) * num_signals, num_signals)) : 0x0000
            end
            reader.position += num_satellites * num_signals
            cell_mask = GalileoHASCellMask(rows, num_satellites, num_signals)
        end
        bits_remaining(reader) >= 3 || return nothing
        nav_message_index = Int(read_bits!(reader, 3))
        push!(
            masks,
            GalileoHASSatelliteMask(
                GNSS_ID,
                satellite_mask,
                signal_mask,
                cell_mask,
                nav_message_index,
                SVIDs,
                signal_indices,
            ),
        )
    end
    # 6 reserved bits close the Mask block (ICD Table 15).
    bits_remaining(reader) >= 6 || return nothing
    read_bits!(reader, 6)
    return GalileoHASMask(mask_id, GalileoHASSatelliteMaskList(masks))
end

"""
    parse_has_orbit_block!(block, reader, context)

Parse the Orbit Corrections block (ICD §5.2.2) into `block` (overwritten).
Returns `nothing` on a truncated message.
"""
function parse_has_orbit_block!(
    block::GalileoHASCorrectionBlock{GalileoHASOrbitCorrection},
    reader::HASBitReader,
    context::GalileoHASBlockContext,
)
    bits_remaining(reader) >= 4 || return nothing
    validity_interval = e6b_validity_interval(Int(read_bits!(reader, 4)))
    corrections = e6b_start_block!(block, context, validity_interval).corrections
    for satellite_mask in context.mask.satellite_masks
        iod_length = e6b_iod_ref_length(satellite_mask.GNSS_ID)
        isnothing(iod_length) && return nothing
        for SVID in satellite_mask.SVIDs
            bits_remaining(reader) >= iod_length + 13 + 12 + 12 || return nothing
            IOD_ref = Int(read_bits!(reader, iod_length))
            # "Data not available" is the most negative value of each field
            # (ICD Table 25): -2^12 radial, -2^11 in-/cross-track.
            radial_raw = read_signed_bits!(reader, 13)
            in_track_raw = read_signed_bits!(reader, 12)
            cross_track_raw = read_signed_bits!(reader, 12)
            push!(
                corrections,
                GalileoHASOrbitCorrection(;
                    GNSS_ID = satellite_mask.GNSS_ID,
                    SVID,
                    IOD_ref,
                    δ_radial = radial_raw == -4096 ? nothing : radial_raw * 0.0025,
                    δ_in_track = in_track_raw == -2048 ? nothing : in_track_raw * 0.008,
                    δ_cross_track = cross_track_raw == -2048 ? nothing :
                                    cross_track_raw * 0.008,
                ),
            )
        end
    end
    return block
end

"""
Decode one 13-bit Delta Clock Correction field into (value in meters,
`do_not_use`) given the constellation's multiplier (ICD Table 31).
"""
function e6b_delta_clock(raw::Int64, multiplier::Int)
    raw == -4096 && return (nothing, false)   # data not available
    raw == 4095 && return (nothing, true)     # satellite shall not be used
    return (raw * 0.0025 * multiplier, false)
end

"""
    parse_has_clock_full_set_block!(block, reader, context)

Parse the Clock Full-Set Corrections block (ICD §5.2.3) into `block`
(overwritten): one 2-bit Delta Clock Multiplier per constellation of the mask,
then one 13-bit correction per corrected satellite. Returns `nothing` on a
truncated message.

The 2-bit multiplier field maps `0…3` to multipliers `1…4` (ICD Table 29); the
raw field is *not* the multiplier, and reading it as one scales every correction
in the block by the wrong integer. Table 29 is the only authority for that
mapping: the ICD's own worked example (the Annex D attachment embedded in the
PDF, transcribed in `test/has_test_vectors.jl`) prints the *raw* Delta Clock
Multiplier fields and the *un-multiplied* corrections, so it pins the field
positions and the 0.0025 m LSB but cannot confirm or refute the multiplier
table either way.
"""
function parse_has_clock_full_set_block!(
    block::GalileoHASCorrectionBlock{GalileoHASClockCorrection},
    reader::HASBitReader,
    context::GalileoHASBlockContext,
)
    num_systems = length(context.mask.satellite_masks)
    bits_remaining(reader) >= 4 + 2 * num_systems || return nothing
    validity_interval = e6b_validity_interval(Int(read_bits!(reader, 4)))
    # The multipliers precede every correction; they are read in place from
    # their positions here rather than collected into a list first.
    multipliers_start = reader.position
    reader.position += 2 * num_systems
    corrections = e6b_start_block!(block, context, validity_interval).corrections
    for (system_index, satellite_mask) in enumerate(context.mask.satellite_masks)
        multiplier =
            Int(peek_bits(reader, multipliers_start + 2 * (system_index - 1), 2)) + 1
        for SVID in satellite_mask.SVIDs
            bits_remaining(reader) >= 13 || return nothing
            value, do_not_use = e6b_delta_clock(read_signed_bits!(reader, 13), multiplier)
            push!(
                corrections,
                GalileoHASClockCorrection(;
                    GNSS_ID = satellite_mask.GNSS_ID,
                    SVID,
                    multiplier,
                    δ_clock = value,
                    do_not_use,
                ),
            )
        end
    end
    return block
end

"""
    parse_has_clock_subset_block!(block, reader, context)

Parse the Clock Subset Corrections block (ICD §5.2.4) into `block`
(overwritten): corrections for a subset of the mask's satellites, each
constellation carrying its own satellite submask whose length is that
constellation's masked-satellite count. Returns `nothing` on a truncated
message, an unknown GNSS ID, or `Nsys_sub == 0`.

The submask is `Nsat` bits — the number of ones in that constellation's Satellite
Mask, not one fewer. GNSS-SDR's loop stops a bit short here and then writes the
corrections into an unsized vector, so its clock-subset support is unfinished and
disabled downstream; this block is decoded on the same footing as the others.

Unlike `Nsys` in the Mask block, the ICD does not declare `Nsys_sub == 0`
reserved (§5.2.4, Table 32) — but a block flagged present that then corrects no
constellation at all is vacuous, and accepting it would publish an empty
`clock_subset_corrections` over a good one. It is treated as malformed, as
GNSS-SDR does.
"""
function parse_has_clock_subset_block!(
    block::GalileoHASCorrectionBlock{GalileoHASClockCorrection},
    reader::HASBitReader,
    context::GalileoHASBlockContext,
)
    bits_remaining(reader) >= 8 || return nothing
    validity_interval = e6b_validity_interval(Int(read_bits!(reader, 4)))
    num_subset_systems = Int(read_bits!(reader, 4))
    num_subset_systems == 0 && return nothing
    corrections = e6b_start_block!(block, context, validity_interval).corrections
    satellite_masks = context.mask.satellite_masks
    for _ = 1:num_subset_systems
        bits_remaining(reader) >= 6 || return nothing
        GNSS_ID = Int(read_bits!(reader, 4))
        multiplier = Int(read_bits!(reader, 2)) + 1
        index = 0
        for (i, satellite_mask) in enumerate(satellite_masks)
            if satellite_mask.GNSS_ID == GNSS_ID
                index = i
                break
            end
        end
        # A subset naming a constellation the mask does not cover leaves the
        # submask length unknown.
        index == 0 && return nothing
        SVIDs = satellite_masks[index].SVIDs
        num_satellites = length(SVIDs)
        bits_remaining(reader) >= num_satellites || return nothing
        # At most 40 bits: the whole submask fits one word, MSB = first SVID.
        submask = read_bits!(reader, num_satellites)
        for (i, SVID) in enumerate(SVIDs)
            (submask >> (num_satellites - i)) & 0x1 == 0x1 || continue
            bits_remaining(reader) >= 13 || return nothing
            value, do_not_use = e6b_delta_clock(read_signed_bits!(reader, 13), multiplier)
            push!(
                corrections,
                GalileoHASClockCorrection(;
                    GNSS_ID,
                    SVID,
                    multiplier,
                    δ_clock = value,
                    do_not_use,
                ),
            )
        end
    end
    return block
end

"""
    e6b_foreach_cell(f, mask) -> Bool

Walk the mask's satellite/signal cells in broadcast order — constellation, then
satellite, then signal — skipping the pairs a Cell Mask excludes, and call
`f(GNSS_ID, SVID, signal_index)` on each. `f` returns `false` to abort (a
truncated message); `e6b_foreach_cell` then returns `false` too, and `true` when
every cell was visited.

The code-bias and phase-bias blocks are indexed by exactly this walk (ICD §5.2.5,
§5.2.6), so it is defined once. Stating it twice would make "the phase biases are
iterated exactly as the code biases" a comment to be trusted rather than a fact
of the code — and a divergence would silently attach every bias to the wrong
cell, which no CRC downstream can catch.
"""
function e6b_foreach_cell(f, mask::GalileoHASMask)
    for satellite_mask in mask.satellite_masks
        cell_mask = satellite_mask.cell_mask
        for (row, SVID) in enumerate(satellite_mask.SVIDs)
            for (column, signal_index) in enumerate(satellite_mask.signal_indices)
                if !isnothing(cell_mask) && !cell_mask[row, column]
                    continue
                end
                f(satellite_mask.GNSS_ID, SVID, signal_index) || return false
            end
        end
    end
    return true
end

"""
    parse_has_code_bias_block!(block, reader, context)

Parse the Code Biases block (ICD §5.2.5) into `block` (overwritten): one 11-bit
bias per cell of the mask, over [`e6b_foreach_cell`](@ref). Returns `nothing`
on a truncated message.
"""
function parse_has_code_bias_block!(
    block::GalileoHASCorrectionBlock{GalileoHASCodeBias},
    reader::HASBitReader,
    context::GalileoHASBlockContext,
)
    bits_remaining(reader) >= 4 || return nothing
    validity_interval = e6b_validity_interval(Int(read_bits!(reader, 4)))
    biases = e6b_start_block!(block, context, validity_interval).corrections
    complete = e6b_foreach_cell(context.mask) do GNSS_ID, SVID, signal_index
        bits_remaining(reader) >= 11 || return false
        raw = read_signed_bits!(reader, 11)
        push!(
            biases,
            GalileoHASCodeBias(;
                GNSS_ID,
                SVID,
                signal_index,
                bias = raw == -1024 ? nothing : raw * 0.02,
            ),
        )
        return true
    end
    complete || return nothing
    return block
end

"""
    parse_has_phase_bias_block!(block, reader, context)

Parse the Phase Biases block (ICD §5.2.6) into `block` (overwritten): per cell
an 11-bit bias immediately followed by its 2-bit Phase Discontinuity Indicator,
over the same [`e6b_foreach_cell`](@ref) walk the code biases use. Returns
`nothing` on a truncated message.
"""
function parse_has_phase_bias_block!(
    block::GalileoHASCorrectionBlock{GalileoHASPhaseBias},
    reader::HASBitReader,
    context::GalileoHASBlockContext,
)
    bits_remaining(reader) >= 4 || return nothing
    validity_interval = e6b_validity_interval(Int(read_bits!(reader, 4)))
    biases = e6b_start_block!(block, context, validity_interval).corrections
    complete = e6b_foreach_cell(context.mask) do GNSS_ID, SVID, signal_index
        bits_remaining(reader) >= 13 || return false
        raw = read_signed_bits!(reader, 11)
        phase_discontinuity_indicator = Int(read_bits!(reader, 2))
        push!(
            biases,
            GalileoHASPhaseBias(;
                GNSS_ID,
                SVID,
                signal_index,
                bias = raw == -1024 ? nothing : raw * 0.01,
                phase_discontinuity_indicator,
            ),
        )
        return true
    end
    complete || return nothing
    return block
end

"""
Run one content-block parser if its header flag is set and the reader is still
aligned (every block before it parsed), returning `(block or nothing, aligned)`.
"""
@inline function e6b_parse_block(
    parse!,
    aligned::Bool,
    flag::Bool,
    block,
    reader::HASBitReader,
    context::GalileoHASBlockContext,
)
    (aligned && flag) || return (nothing, aligned)
    parsed = parse!(block, reader, context)
    return (parsed, !isnothing(parsed))
end

"""
    parse_has_message!(parser, octets, num_octets, message_id, message_type, message_size, masks)
        -> Union{Nothing,Symbol,GalileoHASMessage}

Parse the reassembled HAS message in the first `num_octets` of `octets`,
**overwriting** `parser`'s reader, scratch, blocks and message record. Returns

  - `parser.message` on success — valid until the next parse, its block fields
    pointing into `parser.blocks`,
  - `:mask_unavailable` when the body needs a Mask ID that has not been received
    (the caller holds the message and retries once a matching mask arrives), or
  - `nothing` when the message is not parseable at all (reserved Message Type,
    out-of-range Time Of Hour, or a body that runs out of bits).

`masks` is the store of masks received so far (anything `haskey`/`getindex`-able
by Mask ID, or `nothing`).

Blocks are read strictly in header-flag order (ICD §5.1). A block that fails to
parse ends the body — everything after it in the stream is unaligned — but the
blocks already parsed are kept, since each is self-contained.

The MT1 header's 4 reserved bits are skipped, per the ICD's instruction that
"forward-compatible receivers shall account for the possibility of HAS MT1
messages containing additional non-decodable information at the end of the
decodable one": trailing bits we cannot interpret are simply left unread.
"""
function parse_has_message!(
    parser::GalileoHASParser,
    octets::Vector{UInt8},
    num_octets::Int,
    message_id::Int,
    message_type::Int,
    message_size::Int,
    masks,
)
    # Only Message Type 1 is defined (ICD Table 10).
    message_type == 1 || return nothing
    reader = e6b_rewind!(parser.reader, octets, num_octets)
    bits_remaining(reader) >= 32 || return nothing
    TOH = Int(read_bits!(reader, 12))
    # TOH is seconds into the hour, so 3600-4095 cannot occur (ICD Table 13).
    TOH <= 3599 || return nothing
    mask_flag = read_bits!(reader, 1) == 1
    orbit_flag = read_bits!(reader, 1) == 1
    clock_full_set_flag = read_bits!(reader, 1) == 1
    clock_subset_flag = read_bits!(reader, 1) == 1
    code_bias_flag = read_bits!(reader, 1) == 1
    phase_bias_flag = read_bits!(reader, 1) == 1
    read_bits!(reader, 4)   # reserved (ICD Table 12)
    mask_id = Int(read_bits!(reader, 5))
    IOD_set_id = Int(read_bits!(reader, 5))

    # The message record is overwritten from here on.
    message = parser.message
    message.message_id = message_id
    message.message_type = message_type
    message.message_size = message_size
    message.TOH = TOH
    message.mask_id = mask_id
    message.IOD_set_id = IOD_set_id
    message.mask = nothing
    message.orbit_corrections = nothing
    message.clock_corrections = nothing
    message.clock_subset_corrections = nothing
    message.code_biases = nothing
    message.phase_biases = nothing

    mask = nothing
    if mask_flag
        mask = parse_has_mask_block!(parser, mask_id)
        isnothing(mask) && return nothing
    elseif !isnothing(masks) && haskey(masks, mask_id)
        # Mask Flag = 0 relates the body to a mask already defined by another
        # message with the same Mask ID (ICD §5.1.1.1).
        mask = masks[mask_id]
    end
    if isnothing(mask)
        # Every body block's length is derived from the mask, so without one only
        # a body that has no blocks at all can be parsed.
        (
            orbit_flag ||
            clock_full_set_flag ||
            clock_subset_flag ||
            code_bias_flag ||
            phase_bias_flag
        ) && return :mask_unavailable
        return message
    end
    context = GalileoHASBlockContext(mask, TOH, mask_id, IOD_set_id)

    # Content blocks appear in flag order (ICD Table 14). A block that fails to
    # parse leaves the reader unaligned, so everything after it is unreadable —
    # but the blocks already parsed are self-contained and are kept. Stating that
    # rule once, in `e6b_parse_block`, keeps it from having to be re-established
    # at each of the five call sites (and from being quietly omitted at the last
    # one).
    blocks = parser.blocks
    aligned = true
    orbit, aligned = e6b_parse_block(
        parse_has_orbit_block!,
        aligned,
        orbit_flag,
        blocks.orbit_corrections,
        reader,
        context,
    )
    clock, aligned = e6b_parse_block(
        parse_has_clock_full_set_block!,
        aligned,
        clock_full_set_flag,
        blocks.clock_corrections,
        reader,
        context,
    )
    clock_subset, aligned = e6b_parse_block(
        parse_has_clock_subset_block!,
        aligned,
        clock_subset_flag,
        blocks.clock_subset_corrections,
        reader,
        context,
    )
    code, aligned = e6b_parse_block(
        parse_has_code_bias_block!,
        aligned,
        code_bias_flag,
        blocks.code_biases,
        reader,
        context,
    )
    phase, _ = e6b_parse_block(
        parse_has_phase_bias_block!,
        aligned,
        phase_bias_flag,
        blocks.phase_biases,
        reader,
        context,
    )
    message.mask = mask_flag ? mask : nothing
    message.orbit_corrections = orbit
    message.clock_corrections = clock
    message.clock_subset_corrections = clock_subset
    message.code_biases = code
    message.phase_biases = phase
    return message
end

"""
    parse_has_message(octets, message_id, message_type, message_size, masks)
        -> Union{Nothing,Symbol,GalileoHASMessage}

[`parse_has_message!`](@ref) of the whole of `octets` into freshly allocated
storage — a convenience for inspecting a single message outside a decoder. The
returned message owns its blocks.
"""
parse_has_message(
    octets::Vector{UInt8},
    message_id::Int,
    message_type::Int,
    message_size::Int,
    masks,
) = parse_has_message!(
    GalileoHASParser(),
    octets,
    length(octets),
    message_id,
    message_type,
    message_size,
    masks,
)

# ---- Page accumulation and message assembly (ICD §6.4) ----------------------

"""
    e6b_reassemble_message!(octets, scratch, group) -> Union{Nothing,Int}

Run the HPVRS erasure decode on a complete page group, **overwriting** the
leading `message_size × 53` entries of `octets` with the message octets in page
order (and the working matrices in `scratch`). Returns the number of message
octets, or `nothing` if the decode is impossible.

The ICD's own recipe (§6.4): take the `k` rows of the generator matrix named by
the received Page IDs and its first `k` columns — the trailing `32 - k`
information symbols are known zeros, since "pages Ck+1 … CK contain only
zeroes" (§6.3) — invert that `k × k` matrix over GF(256) once, and apply it to
each of the 53 octet columns. GNSS-SDR instead runs a Berlekamp-Massey erasure
decode over the full 255-symbol codeword; both recover the same message, but the
matrix form is what the ICD specifies and needs no error-locator machinery.
"""
function e6b_reassemble_message!(
    octets::Vector{UInt8},
    scratch::RSErasureScratch,
    group::GalileoHASPageGroup,
)
    k = group.message_size
    # Row `i` of the decode is non-encoded page `M_i`; the message is those
    # pages concatenated, which is the row-by-row order `rs_erasure_decode!`
    # writes in.
    decoded = rs_erasure_decode!(
        octets,
        scratch,
        GALILEO_HAS_GF256,
        E6B_GENERATOR_MATRIX,
        group.page_ids,
        group.octets,
        k,
    )
    isnothing(decoded) && return nothing
    return k * E6B_OCTETS_PER_PAGE
end

"""
    e6b_reassemble_message(group) -> Union{Nothing,Vector{UInt8}}

[`e6b_reassemble_message!`](@ref) into a freshly allocated vector of exactly the
message's octets.
"""
function e6b_reassemble_message(group::GalileoHASPageGroup)
    octets = Vector{UInt8}(undef, group.message_size * E6B_OCTETS_PER_PAGE)
    num_octets =
        e6b_reassemble_message!(octets, RSErasureScratch(E6B_RS_CODE_DIMENSION), group)
    return isnothing(num_octets) ? nothing : octets
end

"""
    e6b_merge_block!(current, spare, new) -> Union{Nothing,GalileoHASCorrectionBlock}

The latest block of one kind after a message: `current` unchanged when the
message did not carry one (`new === nothing`), otherwise `new` copied into —
**overwriting** — `current`, or the preallocated `spare` when there is no
current block yet.
"""
e6b_merge_block!(current, spare, new) =
    isnothing(new) ? current : overwrite!(something(current, spare), new)

"""
    e6b_merge_message!(data, storage, message) -> GalileoE6BData

Fold a decoded [`GalileoHASMessage`](@ref) into the accumulated data container
by **overwriting** its containers in place: publish it as the latest message,
remember any mask it defined, and replace the latest instance of each content
block it carried. Containers `data` does not have yet are taken from the
preallocated `storage` (the cache's `storage.raw`). `message` is copied, so the
parser may overwrite it afterwards.
"""
function e6b_merge_message!(
    data::GalileoE6BData,
    storage::GalileoE6BData,
    message::GalileoHASMessage,
)
    masks = data.masks
    mask = message.mask
    if !isnothing(mask)
        masks = writable_container(data.masks, storage.masks)
        set!(masks, mask.mask_id, mask)
    end
    orbit = e6b_merge_block!(
        data.orbit_corrections,
        storage.orbit_corrections,
        message.orbit_corrections,
    )
    clock = e6b_merge_block!(
        data.clock_corrections,
        storage.clock_corrections,
        message.clock_corrections,
    )
    clock_subset = e6b_merge_block!(
        data.clock_subset_corrections,
        storage.clock_subset_corrections,
        message.clock_subset_corrections,
    )
    code = e6b_merge_block!(data.code_biases, storage.code_biases, message.code_biases)
    phase = e6b_merge_block!(data.phase_biases, storage.phase_biases, message.phase_biases)
    latest = e6b_overwrite_message!(
        something(data.message, storage.message),
        message,
        orbit,
        clock,
        clock_subset,
        code,
        phase,
    )
    GalileoE6BData(data.HAS_status, latest, masks, orbit, clock, clock_subset, code, phase)
end

"""
Drop the in-flight state that has aged past the ICD's 150 s limit (§6.4.1): page
groups that never completed, and a held orphan message whose mask never arrived.

One C/NAV page is one second, so the count of accepted pages is the clock — see
`E6B_MESSAGE_TIMEOUT_PAGES` for why lost pages make that an approximation rather
than an equality.

The orphan expires on the same clock as the page groups for the same reason it is
held at all: it is worth re-parsing because a mask is due within a few messages.
Past the timeout the corrections it carries have outlived their validity
intervals, so publishing them on a late mask would be worse than dropping them.
"""
function e6b_expire_stale!(cache::GalileoE6BCache)
    now = cache.page_counter[]
    held = cache.pending_message[]
    if !isnothing(held) && now - held.opened_at >= E6B_MESSAGE_TIMEOUT_PAGES
        cache.pending_message[] = nothing
    end
    groups = cache.page_groups
    for message_id = 0:(E6B_NUM_MESSAGE_IDS-1)
        if haskey(groups, message_id) &&
           now - groups[message_id].opened_at >= E6B_MESSAGE_TIMEOUT_PAGES
            unset!(groups, message_id)
        end
    end
    return cache
end

"""
Add one HAS Encoded Page to its Message ID's group, returning the group if it is
now complete (ready for the RS decode) and `nothing` otherwise. `unpack!` is
called with the group's newly claimed octet row to fill in.

A Message ID is reused over time, so a page whose Message Type or Message Size
disagrees with the group's starts a fresh group — the previous message is gone,
and the slot's preallocated group is **overwritten** (`e6b_reopen_group!`).
Duplicate Page IDs are ignored: the RS decode needs `k` *distinct* rows.
"""
function e6b_collect_page!(
    unpack!,
    cache::GalileoE6BCache,
    message_id::Int,
    message_type::Int,
    message_size::Int,
    page_id::Int,
)
    groups = cache.page_groups
    if haskey(groups, message_id) &&
       groups[message_id].message_type == message_type &&
       groups[message_id].message_size == message_size
        group = groups[message_id]
    else
        # Every slot owns a preallocated group, occupied or not.
        group = e6b_reopen_group!(
            groups.values[message_id+1],
            message_type,
            message_size,
            cache.page_counter[],
        )
        set!(groups, message_id, group)
    end
    page_id in group.page_ids && return nothing
    # A group is complete at `message_size` ≤ 32 distinct pages, so it never
    # grows past its 32 preallocated rows.
    push!(group.page_ids, page_id)
    # `unpack!` writes the page's octets into the row that has just been claimed;
    # taking a callback rather than a vector lets the caller unpack in place.
    unpack!(@view group.octets[length(group.page_ids), :])
    return length(group.page_ids) == group.message_size ? group : nothing
end

"""
    decode_syncro_sequence(state::GNSSDecoderState{<:GalileoE6BData}, sync::GalileoE6BSync)

Consume one CRC-validated C/NAV page: read its HAS Page Header, discard dummy
and unusable pages, add its HAS Encoded Page to the store for its Message ID,
and — once `MS` distinct pages are held — recover and parse the HAS message.
Overwrites the cache's page store and scratch buffers and `raw_data`'s
containers in place.

Pages are dropped, before reaching the store, when

  - the header is the dummy-page pattern `hex[AF3BC3]` (ICD §2.4.1) — expected
    during nominal operation from satellites with no HAS data to send;
  - the HAS Status is `has_status_reserved` or `has_do_not_use` (ICD Table 9).
    On `has_do_not_use` the ICD requires more than dropping the page: users
    "shall stop using HAS from all satellites and discard previously received
    messages", so the page store and every published correction are cleared;
  - the Page ID is 0, which is reserved (ICD Table 8);
  - the Message Type is not 1, the only type defined (ICD Table 10).
"""
function decode_syncro_sequence(
    state::GNSSDecoderState{<:GalileoE6BData},
    sync::GalileoE6BSync,
)
    cache = state.cache
    page = sync.page
    cache.page_counter[] += 1
    e6b_expire_stale!(cache)
    # Every page reaching here framed correctly and passed CRC-24Q, so it re-arms
    # the framework's symbol counter whether or not its contents turn out usable:
    # the field counts symbols since the last valid sync, which is how a consumer
    # tells a locked decoder from an unsynchronised one. Without this it would
    # stay `nothing` forever on E6-B, since `decode` only ever increments a
    # counter that is already armed. `preamble_length` is the trailing next-page
    # sync `drain_after_sync!` leaves buffered (cf. `gps/cnav.jl`).
    state = GNSSDecoderState(
        state;
        num_bits_after_valid_syncro_sequence = state.constants.preamble_length,
    )

    header = UInt32(get_bits(page, E6B_PAGE_BITS, E6B_HEADER_START, 24))
    # Dummy pages carry no HAS data at all (ICD §2.4.1); their header fields are
    # not meaningful, so they are dropped before the status is even read.
    header == E6B_DUMMY_PAGE_HEADER && return state

    HAS_status = HASStatus(get_bits(page, E6B_PAGE_BITS, E6B_HEADER_START, 2))
    if HAS_status == has_do_not_use
        # Table 9 requires more than dropping the page: everything received so
        # far goes, so the status is carried into an otherwise empty container
        # rather than merged into the existing one. (The preallocated containers
        # stay in the cache's storage, to be overwritten by the next message.)
        empty!(cache.page_groups)
        cache.pending_message[] = nothing
        return GNSSDecoderState(
            state;
            raw_data = GalileoE6BData(; HAS_status),
            data = GalileoE6BData(; HAS_status),
        )
    end
    state = GNSSDecoderState(state; raw_data = GalileoE6BData(state.raw_data; HAS_status))
    HAS_status == has_status_reserved && return state

    message_type = Int(get_bits(page, E6B_PAGE_BITS, E6B_HEADER_START + 4, 2))
    message_id = Int(get_bits(page, E6B_PAGE_BITS, E6B_HEADER_START + 6, 5))
    # MS is broadcast as "size - 1": "0" = 1 page … "31" = 32 pages (ICD Table 8).
    message_size = Int(get_bits(page, E6B_PAGE_BITS, E6B_HEADER_START + 11, 5)) + 1
    page_id = Int(get_bits(page, E6B_PAGE_BITS, E6B_HEADER_START + 16, 8))
    # Page ID 0 is reserved (ICD Table 8), and Message Type 1 is the only type
    # defined (ICD Table 10).
    (page_id == 0 || message_type != 1) && return state
    # Page IDs message_size+1 … 32 are the zero padding of a short message:
    # "Pages Ck+1, …, CK contain only zeroes and are excluded from transmission"
    # (ICD §6.3). They carry no information and would make the decoding matrix
    # singular, so they are never collected.
    (message_size < page_id <= E6B_RS_CODE_DIMENSION) && return state

    group = e6b_collect_page!(cache, message_id, message_type, message_size, page_id) do row
        # Unpack the 53 HAS Encoded Page octets straight into the group's
        # matrix row rather than through a temporary vector.
        @inbounds for j = 1:E6B_OCTETS_PER_PAGE
            row[j] = UInt8(
                get_bits(page, E6B_PAGE_BITS, E6B_ENCODED_PAGE_START + 8 * (j - 1), 8),
            )
        end
    end
    isnothing(group) && return state

    num_octets = e6b_reassemble_message!(cache.message_octets, cache.rs_scratch, group)
    unset!(cache.page_groups, message_id)
    isnothing(num_octets) && return state

    return e6b_apply_message(state, message_id, message_type, message_size, num_octets)
end

"""
Hold a message whose mask has not arrived: **overwrite** the cache's
preallocated held-message record with it (its octets copied out of the
reassembly buffer, which the next message reuses) and point `pending_message`
at it. One slot is enough: only the most recent orphan is worth keeping.
"""
function e6b_hold_message!(
    cache::GalileoE6BCache,
    message_id::Int,
    message_type::Int,
    message_size::Int,
    num_octets::Int,
)
    held = something(cache.pending_message[], cache.held_message)
    held.message_id = message_id
    held.message_type = message_type
    held.message_size = message_size
    held.opened_at = cache.page_counter[]
    resize!(held.octets, num_octets)
    copyto!(held.octets, 1, cache.message_octets, 1, num_octets)
    cache.pending_message[] = held
    return cache
end

"""
Parse the reassembled message in the cache's `message_octets` and fold it into
`raw_data` (overwriting its containers), handling the mask-not-yet-received
case by holding the message for a later retry, and retrying any previously held
message once a new mask has been learnt.
"""
function e6b_apply_message(
    state::GNSSDecoderState{<:GalileoE6BData},
    message_id::Int,
    message_type::Int,
    message_size::Int,
    num_octets::Int,
)
    cache = state.cache
    parser = cache.parser
    storage = cache.storage.raw
    octets = cache.message_octets
    parsed = parse_has_message!(
        parser,
        octets,
        num_octets,
        message_id,
        message_type,
        message_size,
        state.raw_data.masks,
    )
    if parsed === :mask_unavailable
        # Hold the message; HAS broadcasts a defining mask every few messages,
        # and re-parsing then recovers corrections that would otherwise be lost.
        e6b_hold_message!(cache, message_id, message_type, message_size, num_octets)
        return state
    end
    parsed isa GalileoHASMessage || return state
    raw = state.raw_data
    # A message that carried a *new* mask may unlock the held-back one; a message
    # that merely referenced an existing one cannot, so only check when
    # `parsed.mask` is present.
    #
    # ORDER MATTERS. The held message completed *earlier* than `parsed`, so it is
    # merged *first* and `parsed` folded over the top: `e6b_merge_message!` lets
    # the later of two messages win every field it carries, so merging the orphan
    # last would publish it as `message` (documented as the most recently
    # completed one) and would replace a fresher correction block with a staler
    # one of the same kind. Retrying needs the mask `parsed` carries but not the
    # rest of it, so the mask is stored on its own for the retry — `parsed`
    # stores it there anyway when it is merged below.
    held = cache.pending_message[]
    new_mask = parsed.mask
    if !isnothing(new_mask) && !isnothing(held)
        masks = writable_container(raw.masks, storage.masks)
        set!(masks, new_mask.mask_id, new_mask)
        raw = GalileoE6BData(raw; masks)
        retried = parse_has_message!(
            parser,
            held.octets,
            length(held.octets),
            held.message_id,
            held.message_type,
            held.message_size,
            masks,
        )
        if retried isa GalileoHASMessage
            cache.pending_message[] = nothing
            raw = e6b_merge_message!(raw, storage, retried)
        end
        # The retry overwrote the parser, so `parsed` is parsed again — the same
        # octets give the same message (it carries its own mask, so the store
        # just updated does not change it).
        parsed = parse_has_message!(
            parser,
            octets,
            num_octets,
            message_id,
            message_type,
            message_size,
            raw.masks,
        )
        parsed isa GalileoHASMessage || return GNSSDecoderState(state; raw_data = raw)
    end
    return GNSSDecoderState(state; raw_data = e6b_merge_message!(raw, storage, parsed))
end

"""
    validate_data(state::GNSSDecoderState{<:GalileoE6BData})

Publish `raw_data` as `data` unconditionally — as a copy: every container of
`raw_data` is copied into (**overwriting**) the matching preallocated container
of the cache's `storage.validated`, so `data` shares nothing with `raw_data` and
holds still while later pages are decoded into it. The published message's
block fields point at the published blocks, as they do in `raw_data`.

Every other decoder here withholds `data` until a cross-message consistency
check passes — matching issues of data, a repeated broadcast, a plausible time
of week. C/NAV has no such check to make and needs none: a HAS message is
reassembled only from pages that each passed CRC-24Q, under one Message ID with
one Message Size, and the ICD's whole erasure-channel model rests on that gate
(§6.4). There is nothing left to corroborate, and withholding corrections that
carry their own validity intervals would only make them stale.
"""
function validate_data(state::GNSSDecoderState{<:GalileoE6BData})
    raw = state.raw_data
    validated = state.cache.storage.validated
    orbit = publish!(validated.orbit_corrections, raw.orbit_corrections)
    clock = publish!(validated.clock_corrections, raw.clock_corrections)
    clock_subset =
        publish!(validated.clock_subset_corrections, raw.clock_subset_corrections)
    code = publish!(validated.code_biases, raw.code_biases)
    phase = publish!(validated.phase_biases, raw.phase_biases)
    message =
        isnothing(raw.message) ? nothing :
        e6b_overwrite_message!(
            validated.message,
            raw.message,
            orbit,
            clock,
            clock_subset,
            code,
            phase,
        )
    data = GalileoE6BData(
        raw.HAS_status,
        message,
        publish!(validated.masks, raw.masks),
        orbit,
        clock,
        clock_subset,
        code,
        phase,
    )
    return GNSSDecoderState(state; data)
end

"""
$(TYPEDSIGNATURES)

Report whether this satellite's High Accuracy Service data is usable.

E6-B C/NAV broadcasts no health flag for the transmitting satellite — it carries
corrections, not that satellite's own navigation data — so what this reports is
the *service* status from the HAS Page Header (ICD Table 9): `true` only while
HAS is in operational mode, when nominal performance is expected.

Test mode (`has_test_mode`) returns `false`: the data decodes normally and is
published in `state.data`, but the ICD warns that "nominal performance may not
be met", so it should not be trusted silently. `has_do_not_use` also returns
`false`, and additionally causes every previously received message to be
discarded (see `decode_syncro_sequence`).

# Arguments

  - `state::GNSSDecoderState{<:GalileoE6BData}`: Galileo E6-B decoder state

# Returns

  - `Bool`: `true` iff the last valid page reported operational HAS status

# See Also

  - [`GalileoE6BDecoderState`](@ref): Create decoder state
  - [`HASStatus`](@ref): The broadcast status values
"""
function is_sat_healthy(state::GNSSDecoderState{<:GalileoE6BData})
    state.data.HAS_status == has_operational_mode
end
