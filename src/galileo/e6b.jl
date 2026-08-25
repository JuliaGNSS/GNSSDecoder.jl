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
# arithmetic, and the generator matrix live in `src/reed_solomon.jl`.
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
    e6b_iod_ref_length(gnss_id) -> Union{Nothing,Int}

Width in bits of the Reference IOD field for a GNSS (ICD Table 26): 8 for GPS
(IODE/IODC), 10 for Galileo (IODnav). `nothing` for a reserved GNSS ID.
"""
e6b_iod_ref_length(gnss_id::Integer) =
    gnss_id == E6B_GNSS_ID_GPS ? 8 : gnss_id == E6B_GNSS_ID_GALILEO ? 10 : nothing

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

# ---- Decoded HAS records (ICD §5) -------------------------------------------

"""
    GalileoHASSatelliteMask

The HAS mask for one constellation: which satellites are corrected, which
signals carry biases, and which broadcast navigation message the orbit and clock
corrections refer to.

The raw ICD bit fields are kept alongside the expanded lists, since the raw
masks are what later blocks' lengths are computed from.

# Fields

  - `gnss_id::Int`: GNSS index — 0 = GPS, 2 = Galileo (Table 18)
  - `satellite_mask::UInt64`: 40-bit Satellite Mask, MSB = satellite index 0 (Table 19)
  - `signal_mask::UInt16`: 16-bit Signal Mask, MSB = signal index 0 (Table 20)
  - `cell_mask::Union{Nothing,Matrix{Bool}}`: `Nsat × Nsig` Cell Mask, or `nothing` when the Cell Mask Availability Flag is 0 (biases then cover every masked satellite/signal pair)
  - `nav_message_index::Int`: Navigation Message Index — 0 = I/NAV (Galileo) / LNAV (GPS) (Table 21)
  - `svids::Vector{Int}`: Satellite IDs of the masked satellites (satellite index + 1, i.e. Galileo SVID / GPS PRN)
  - `signal_indices::Vector{Int}`: Signal indices of the masked signals (0-based, per Table 20)

# Reference

Galileo HAS SIS ICD, Issue 1.0, Tables 16-21
"""
Base.@kwdef struct GalileoHASSatelliteMask
    gnss_id::Int
    satellite_mask::UInt64
    signal_mask::UInt16
    cell_mask::Union{Nothing,Matrix{Bool}} = nothing
    nav_message_index::Int
    svids::Vector{Int}
    signal_indices::Vector{Int}
end

Base.:(==)(a::GalileoHASSatelliteMask, b::GalileoHASSatelliteMask) = fields_equal(a, b)

"""
    GalileoHASMask

One complete HAS Mask block: the per-constellation masks it defines, under the
Mask ID that later messages reference.

# Fields

  - `mask_id::Int`: Mask ID this mask defines (0-31)
  - `satellite_masks::Vector{GalileoHASSatelliteMask}`: one entry per corrected constellation, in broadcast order

# Reference

Galileo HAS SIS ICD, Issue 1.0, Table 15
"""
Base.@kwdef struct GalileoHASMask
    mask_id::Int
    satellite_masks::Vector{GalileoHASSatelliteMask}
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

  - `gnss_id::Int`: GNSS index the satellite belongs to (Table 18)
  - `svid::Int`: Satellite ID (Galileo SVID / GPS PRN)
  - `IOD_ref::Int`: Reference IOD of the corrected broadcast navigation data — IODnav for Galileo, IODE/IODC for GPS (Table 26)
  - `δ_radial::Union{Nothing,Float64}`: Delta Radial correction (meters, LSB 0.0025)
  - `δ_in_track::Union{Nothing,Float64}`: Delta In-Track correction (meters, LSB 0.008)
  - `δ_cross_track::Union{Nothing,Float64}`: Delta Cross-Track correction (meters, LSB 0.008)

# Reference

Galileo HAS SIS ICD, Issue 1.0, Tables 24-25
"""
Base.@kwdef struct GalileoHASOrbitCorrection
    gnss_id::Int
    svid::Int
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

  - `gnss_id::Int`: GNSS index the satellite belongs to (Table 18)
  - `svid::Int`: Satellite ID (Galileo SVID / GPS PRN)
  - `multiplier::Int`: Delta Clock Multiplier applied, 1-4 (Table 29)
  - `δ_clock::Union{Nothing,Float64}`: Delta clock correction (meters, LSB 0.0025 × `multiplier`)
  - `do_not_use::Bool`: The satellite shall not be used (Table 31)

# Reference

Galileo HAS SIS ICD, Issue 1.0, Tables 28-34
"""
Base.@kwdef struct GalileoHASClockCorrection
    gnss_id::Int
    svid::Int
    multiplier::Int
    δ_clock::Union{Nothing,Float64} = nothing
    do_not_use::Bool = false
end

"""
    GalileoHASCodeBias

Code bias for one satellite/signal cell (ICD §7.4). `bias` is `nothing` for the
"data not available" sentinel (raw -1024).

# Fields

  - `gnss_id::Int`: GNSS index (Table 18)
  - `svid::Int`: Satellite ID (Galileo SVID / GPS PRN)
  - `signal_index::Int`: Signal index within the constellation's Signal Mask (0-based, Table 20)
  - `bias::Union{Nothing,Float64}`: Code bias (meters, LSB 0.02)

# Reference

Galileo HAS SIS ICD, Issue 1.0, Tables 36-37
"""
Base.@kwdef struct GalileoHASCodeBias
    gnss_id::Int
    svid::Int
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

  - `gnss_id::Int`: GNSS index (Table 18)
  - `svid::Int`: Satellite ID (Galileo SVID / GPS PRN)
  - `signal_index::Int`: Signal index within the constellation's Signal Mask (0-based, Table 20)
  - `bias::Union{Nothing,Float64}`: Phase bias (cycles, LSB 0.01)
  - `phase_discontinuity_indicator::Int`: Phase Discontinuity Indicator (0-3)

# Reference

Galileo HAS SIS ICD, Issue 1.0, Tables 39-40
"""
Base.@kwdef struct GalileoHASPhaseBias
    gnss_id::Int
    svid::Int
    signal_index::Int
    bias::Union{Nothing,Float64} = nothing
    phase_discontinuity_indicator::Int
end

"""
    GalileoHASCorrectionBlock{T}

One HAS MT1 content block: a validity interval, the message header context it
was broadcast under, and the per-satellite (or per-cell) corrections themselves.

Every block carries its own Validity Interval Index (ICD §5.2.2.1) starting at
the message's Time Of Hour, so blocks of one message can — and routinely do —
expire at different times. `mask_id` and `IOD_set_id` identify the satellite set
and the broadcast-ephemeris issue the corrections apply to (ICD §7.6).

# Fields

  - `TOH::Int`: Time Of Hour of the message that carried this block (seconds into the GST hour, 0-3599)
  - `mask_id::Int`: Mask ID the corrections are keyed to
  - `IOD_set_id::Int`: IOD Set ID the corrections are keyed to
  - `validity_interval::Union{Nothing,Int}`: Validity interval in seconds from `TOH` (`nothing` for the reserved index 15)
  - `corrections::Vector{T}`: the block's entries, in broadcast order

# Reference

Galileo HAS SIS ICD, Issue 1.0, Tables 22, 27, 32, 35, 38
"""
Base.@kwdef struct GalileoHASCorrectionBlock{T}
    TOH::Int
    mask_id::Int
    IOD_set_id::Int
    validity_interval::Union{Nothing,Int}
    corrections::Vector{T}
end

# The `corrections` vector makes the default `===` too strict; see `fields_equal`.
Base.:(==)(a::GalileoHASCorrectionBlock, b::GalileoHASCorrectionBlock) = fields_equal(a, b)

"""
    GalileoHASMessage

One completely received and decoded HAS message — the atomic unit the ICD
defines, reassembled from `MS` HAS Encoded Pages by the RS erasure decoder.

Only Message Type 1 is specified (ICD Table 10), so in practice
`message_type == 1` and the six block fields are populated according to the
flags of the MT1 header; a flag that was 0 leaves its field `nothing`.

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
Base.@kwdef struct GalileoHASMessage
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

# Same reason as `GalileoHASCorrectionBlock` above: the nested blocks and the
# mask carry vectors.
Base.:(==)(a::GalileoHASMessage, b::GalileoHASMessage) = fields_equal(a, b)

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
  - `masks::Dictionary{Int,GalileoHASMask}`: every mask received so far, keyed by Mask ID
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
    masks::Union{Nothing,Dictionary{Int,GalileoHASMask}} = nothing
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

function GalileoE6BData(
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

# The mutable `masks::Dictionary` and the `Vector` fields inside the correction
# blocks make the default struct `==` (which falls back to `===`) too strict.
Base.:(==)(a::GalileoE6BData, b::GalileoE6BData) = fields_equal(a, b)

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

# ---- Page store (ICD §6.4) ---------------------------------------------------

"""
$(TYPEDEF)

Encoded pages collected so far for one Message ID, awaiting the `message_size`
distinct pages the Reed-Solomon erasure decoder needs (ICD §6.4).

Mutable and mutated in place inside the decoder cache: it is exactly the kind of
"still partial" state `CONTEXT.md` says belongs there.

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
    are held, and reaching `message_size` is what completes the group
    """
    page_ids::Vector{Int}
    """
    Collected pages, `message_size × 53` octets; row `i` belongs to `page_ids[i]`,
    so rows beyond `length(page_ids)` are not yet filled
    """
    octets::Matrix{UInt8}
end

GalileoHASPageGroup(message_type::Int, message_size::Int, opened_at::Int) =
    GalileoHASPageGroup(
        message_type,
        message_size,
        opened_at,
        Int[],
        zeros(UInt8, message_size, E6B_OCTETS_PER_PAGE),
    )

# Being *mutable* makes the default `==` reference equality outright — not just
# for the `Vector`/`Matrix` fields — so without this two decoders fed the same
# stream would compare unequal for as long as either holds a partial message,
# which is most of the time. `GalileoE6BCache` compares its page store by value,
# so this is what makes that comparison mean anything.
Base.:(==)(a::GalileoHASPageGroup, b::GalileoHASPageGroup) = fields_equal(a, b)

"""
$(TYPEDEF)

A reassembled HAS message held back because its body references a Mask ID that
has not been received yet (ICD §5.1.1.1) — every block's length derives from the
mask, so nothing past the header can be parsed without it.

One slot is enough: HAS broadcasts a defining mask every few messages, so the
orphan worth keeping is the newest. It carries the same `opened_at` page stamp as
[`GalileoHASPageGroup`](@ref) and expires on the same clock, so a mask arriving
long afterwards cannot resurrect corrections whose validity intervals have run
out.

# Fields

$(TYPEDFIELDS)
"""
struct GalileoHASPendingMessage
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
    The reassembled message octets, awaiting a mask
    """
    octets::Vector{UInt8}
end

# `octets` is a Vector, so the default struct `==` would be reference equality.
Base.:(==)(a::GalileoHASPendingMessage, b::GalileoHASPendingMessage) = fields_equal(a, b)

# ---- Cache -------------------------------------------------------------------

"""
$(TYPEDEF)

Per-decoder cache for Galileo E6-B.

Beyond the shared soft-symbol deque this holds the page-level FEC scratch, the
long-lived Viterbi decoder, the RS generator matrix (built once — 255 × 32
octets), and the in-flight HAS page store: the partial-message state that is the
whole point of the HPVRS outer layer.

The page store, the page counter, and the held-back message are mutated in
place (the `CONTEXT.md` rule: immutable decoded fields in `data`/`raw_data`,
genuinely-in-flight state in the cache).

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
    AFF3CT K=7 NSC Viterbi decoder (K = 486, N = 984), built once and reused across pages
    """
    viterbi_decoder::Aff3ct.ConvViterbiDecoder
    """
    Encoded pages collected per Message ID, awaiting completion
    """
    page_groups::Dictionary{Int,GalileoHASPageGroup}
    """
    Count of accepted C/NAV pages — one per second, the clock for the ICD's 150 s message timeout
    """
    page_counter::Base.RefValue{Int}
    """
    A completed message whose body needs a Mask ID not yet received, held for re-parsing
    """
    pending_message::Base.RefValue{Union{Nothing,GalileoHASPendingMessage}}
end

GalileoE6BCache() = GalileoE6BCache(
    CircularDeque{Float32}(E6B_WINDOW_SYMBOLS),
    Vector{Float32}(undef, E6B_ENCODED_SYMBOLS),
    Aff3ct.ConvViterbiDecoder(E6B_PAGE_BITS, E6B_ENCODED_SYMBOLS, GALILEO_VITERBI_POLY),
    Dictionary{Int,GalileoHASPageGroup}(),
    Ref(0),
    Ref{Union{Nothing,GalileoHASPendingMessage}}(nothing),
)

# The Viterbi handle and the FEC scratch window are derived, not state, so they
# are excluded — but the page store, the page counter and the held orphan are
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

# Arguments

  - `state::GNSSDecoderState{<:GalileoE6BData}`: Current Galileo E6-B decoder state

# Returns

  - `GNSSDecoderState{<:GalileoE6BData}`: Reset decoder state with cleared buffers

# See Also

  - [`GalileoE6BDecoderState`](@ref): Create a fresh decoder state
  - [`decode`](@ref): Continue decoding after reset
"""
function reset_decoder_state(state::GNSSDecoderState{<:GalileoE6BData})
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
    galileo_e6b_viterbi(decoder, soft_page) -> UInt512

Recover one C/NAV page's 486 information bits from `soft_page` — the 984
polarity-corrected `Float32` LLR soft symbols following the page's 16-symbol
sync pattern.

Thin wrapper over the shared [`galileo_viterbi`](@ref) with C/NAV's 8×123
interleaver shape (HAS SIS ICD Table 4) and `UInt512` payload type; the FEC
itself is the Galileo-wide K=7 NSC code (ICD Table 3), G2 inverted.
"""
galileo_e6b_viterbi(
    decoder::Aff3ct.ConvViterbiDecoder,
    soft_page::AbstractVector{Float32},
) = galileo_viterbi(decoder, soft_page, E6B_INTERLEAVER_COLUMNS, UInt512)

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
    page = galileo_e6b_viterbi(state.cache.viterbi_decoder, window)
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

# Fields

$(TYPEDFIELDS)
"""
mutable struct HASBitReader
    """
    Reassembled message octets
    """
    octets::Vector{UInt8}
    """
    Number of bits consumed so far
    """
    position::Int
end

HASBitReader(octets::Vector{UInt8}) = HASBitReader(octets, 0)

"""
Bits left unread in the message.
"""
bits_remaining(reader::HASBitReader) = 8 * length(reader.octets) - reader.position

"""
    read_bits!(reader, num_bits) -> UInt64

Consume the next `num_bits` (at most 64) MSB-first. Callers must have checked
[`bits_remaining`](@ref) first.
"""
function read_bits!(reader::HASBitReader, num_bits::Int)
    value = UInt64(0)
    @inbounds for _ = 1:num_bits
        byte = reader.octets[(reader.position>>3)+1]
        bit = (byte >> (7 - (reader.position & 7))) & 0x01
        value = (value << 1) | UInt64(bit)
        reader.position += 1
    end
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

# ---- MT1 parsing (ICD §5) ----------------------------------------------------

"""
$(TYPEDEF)

Everything a content-block parser needs beyond the bit reader: the mask its
lengths derive from, and the three MT1 header fields every block it produces is
stamped with.

All five content blocks of an MT1 body take exactly this, so it travels as one
value rather than four repeated parameters, and
[`has_correction_block`](@ref) turns it plus a block's own two fields into the
block itself. The mask belongs here for the same reason the header fields do:
*every* block's length is derived from it (HAS SIS ICD, Issue 1.0, §5.2), so no
parser can run without one.

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
    has_correction_block(context, validity_interval, corrections) -> GalileoHASCorrectionBlock

Stamp a parsed block's entries with the header context they were broadcast
under. The element type is taken from `corrections`, so each parser names its
correction type once, where it builds them.
"""
has_correction_block(
    context::GalileoHASBlockContext,
    validity_interval::Union{Nothing,Int},
    corrections::Vector{T},
) where {T} = GalileoHASCorrectionBlock{T}(;
    context.TOH,
    context.mask_id,
    context.IOD_set_id,
    validity_interval,
    corrections,
)

"""
    e6b_expand_mask(mask, width; first_index = 0) -> Vector{Int}

Expand a `width`-bit HAS mask into the indices its set bits name, MSB = index 0,
counting from `first_index`.

Both masks the ICD defines are this walk and differ only in those two numbers:
the 40-bit Satellite Mask names Galileo SVID / GPS PRN `index + 1` (Table 19),
and the 16-bit Signal Mask names 0-based signal indices (Table 20).
"""
function e6b_expand_mask(mask::Unsigned, width::Int; first_index::Int = 0)
    indices = Int[]
    for index = 0:(width-1)
        if (mask >> (width - 1 - index)) & one(mask) == 1
            push!(indices, index + first_index)
        end
    end
    return indices
end

"""
    parse_has_mask_block!(reader, mask_id) -> Union{Nothing,GalileoHASMask}

Parse the Mask block (ICD §5.2.1). Returns `nothing` if the message runs out of
bits, declares the reserved `Nsys` value 0, or names a reserved GNSS ID — in
every case every following block's length is unknown, so the body parse must
stop.

Note the 6 reserved bits after the per-constellation masks are consumed
unconditionally, whereas GNSS-SDR consumes them (and the `Nsys` field itself)
only when `Nsys != 0`. That is not a framing bug on its part — it clears its
`have_mask` flag and gates every later block on it, so it bails out rather than
desynchronising — but this parser rejects `Nsys == 0` outright instead, which is
what the ICD asks for.
"""
function parse_has_mask_block!(reader::HASBitReader, mask_id::Int)
    bits_remaining(reader) >= 4 || return nothing
    num_systems = Int(read_bits!(reader, 4))
    # `Nsys` takes "values from 1 to 15 (value "0" is Reserved)" (ICD §5.2.1). A
    # zero would otherwise yield a mask correcting nothing, which is worse than
    # useless: it would be cached under a real Mask ID and then *satisfy* later
    # `Mask Flag = 0` messages, so they would parse to empty correction blocks
    # instead of being held until the genuine mask arrives.
    num_systems == 0 && return nothing
    masks = GalileoHASSatelliteMask[]
    for _ = 1:num_systems
        bits_remaining(reader) >= 4 + 40 + 16 + 1 || return nothing
        gnss_id = Int(read_bits!(reader, 4))
        # A reserved GNSS ID makes the Reference IOD width, and therefore every
        # later block length, unknowable (ICD Table 26).
        isnothing(e6b_iod_ref_length(gnss_id)) && return nothing
        satellite_mask = read_bits!(reader, 40)
        signal_mask = UInt16(read_bits!(reader, 16))
        cell_mask_available = read_bits!(reader, 1) == 1
        svids = e6b_expand_mask(satellite_mask, 40; first_index = 1)
        signal_indices = e6b_expand_mask(signal_mask, 16)
        cell_mask = nothing
        if cell_mask_available
            # L_CM = Nsig · Nsat, read satellite-major: the table is "read from
            # left to right and from top to bottom" of Nsat rows by Nsig columns
            # (ICD §5.2.1.5, Eq. 3).
            num_cells = length(svids) * length(signal_indices)
            bits_remaining(reader) >= num_cells + 3 || return nothing
            cells = Matrix{Bool}(undef, length(svids), length(signal_indices))
            for row = 1:length(svids), column = 1:length(signal_indices)
                cells[row, column] = read_bits!(reader, 1) == 1
            end
            cell_mask = cells
        end
        bits_remaining(reader) >= 3 || return nothing
        nav_message_index = Int(read_bits!(reader, 3))
        push!(
            masks,
            GalileoHASSatelliteMask(;
                gnss_id,
                satellite_mask,
                signal_mask,
                cell_mask,
                nav_message_index,
                svids,
                signal_indices,
            ),
        )
    end
    # 6 reserved bits close the Mask block (ICD Table 15).
    bits_remaining(reader) >= 6 || return nothing
    read_bits!(reader, 6)
    return GalileoHASMask(; mask_id, satellite_masks = masks)
end

"""
    parse_has_orbit_block!(reader, context)

Parse the Orbit Corrections block (ICD §5.2.2). Returns `nothing` on a
truncated message.
"""
function parse_has_orbit_block!(reader::HASBitReader, context::GalileoHASBlockContext)
    bits_remaining(reader) >= 4 || return nothing
    validity_interval = e6b_validity_interval(Int(read_bits!(reader, 4)))
    corrections = GalileoHASOrbitCorrection[]
    for satellite_mask in context.mask.satellite_masks
        iod_length = e6b_iod_ref_length(satellite_mask.gnss_id)
        isnothing(iod_length) && return nothing
        for svid in satellite_mask.svids
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
                    gnss_id = satellite_mask.gnss_id,
                    svid,
                    IOD_ref,
                    δ_radial = radial_raw == -4096 ? nothing : radial_raw * 0.0025,
                    δ_in_track = in_track_raw == -2048 ? nothing : in_track_raw * 0.008,
                    δ_cross_track = cross_track_raw == -2048 ? nothing :
                                    cross_track_raw * 0.008,
                ),
            )
        end
    end
    return has_correction_block(context, validity_interval, corrections)
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
    parse_has_clock_full_set_block!(reader, context)

Parse the Clock Full-Set Corrections block (ICD §5.2.3): one 2-bit Delta Clock
Multiplier per constellation of the mask, then one 13-bit correction per
corrected satellite. Returns `nothing` on a truncated message.

The 2-bit multiplier field maps `0…3` to multipliers `1…4` (ICD Table 29); the
raw field is *not* the multiplier. HAS SIS ICD Annex D's worked example prints
the raw field (`2` for a constellation whose corrections need ×3), which is easy
to misread as the multiplier itself.
"""
function parse_has_clock_full_set_block!(
    reader::HASBitReader,
    context::GalileoHASBlockContext,
)
    num_systems = length(context.mask.satellite_masks)
    bits_remaining(reader) >= 4 + 2 * num_systems || return nothing
    validity_interval = e6b_validity_interval(Int(read_bits!(reader, 4)))
    multipliers = [Int(read_bits!(reader, 2)) + 1 for _ = 1:num_systems]
    corrections = GalileoHASClockCorrection[]
    for (system_index, satellite_mask) in enumerate(context.mask.satellite_masks)
        multiplier = multipliers[system_index]
        for svid in satellite_mask.svids
            bits_remaining(reader) >= 13 || return nothing
            value, do_not_use = e6b_delta_clock(read_signed_bits!(reader, 13), multiplier)
            push!(
                corrections,
                GalileoHASClockCorrection(;
                    gnss_id = satellite_mask.gnss_id,
                    svid,
                    multiplier,
                    δ_clock = value,
                    do_not_use,
                ),
            )
        end
    end
    return has_correction_block(context, validity_interval, corrections)
end

"""
    parse_has_clock_subset_block!(reader, context)

Parse the Clock Subset Corrections block (ICD §5.2.4): corrections for a subset
of the mask's satellites, each constellation carrying its own satellite submask
whose length is that constellation's masked-satellite count. Returns `nothing`
on a truncated message, an unknown GNSS ID, or `Nsys_sub == 0`.

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
    reader::HASBitReader,
    context::GalileoHASBlockContext,
)
    bits_remaining(reader) >= 8 || return nothing
    validity_interval = e6b_validity_interval(Int(read_bits!(reader, 4)))
    num_subset_systems = Int(read_bits!(reader, 4))
    num_subset_systems == 0 && return nothing
    corrections = GalileoHASClockCorrection[]
    for _ = 1:num_subset_systems
        bits_remaining(reader) >= 6 || return nothing
        gnss_id = Int(read_bits!(reader, 4))
        multiplier = Int(read_bits!(reader, 2)) + 1
        index = findfirst(m -> m.gnss_id == gnss_id, context.mask.satellite_masks)
        # A subset naming a constellation the mask does not cover leaves the
        # submask length unknown.
        isnothing(index) && return nothing
        svids = context.mask.satellite_masks[index].svids
        bits_remaining(reader) >= length(svids) || return nothing
        subset_svids = Int[]
        for svid in svids
            read_bits!(reader, 1) == 1 && push!(subset_svids, svid)
        end
        for svid in subset_svids
            bits_remaining(reader) >= 13 || return nothing
            value, do_not_use = e6b_delta_clock(read_signed_bits!(reader, 13), multiplier)
            push!(
                corrections,
                GalileoHASClockCorrection(;
                    gnss_id,
                    svid,
                    multiplier,
                    δ_clock = value,
                    do_not_use,
                ),
            )
        end
    end
    return has_correction_block(context, validity_interval, corrections)
end

"""
    e6b_foreach_cell(f, mask) -> Bool

Walk the mask's satellite/signal cells in broadcast order — constellation, then
satellite, then signal — skipping the pairs a Cell Mask excludes, and call
`f(gnss_id, svid, signal_index)` on each. `f` returns `false` to abort (a
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
        for (row, svid) in enumerate(satellite_mask.svids)
            for (column, signal_index) in enumerate(satellite_mask.signal_indices)
                if !isnothing(satellite_mask.cell_mask) &&
                   !satellite_mask.cell_mask[row, column]
                    continue
                end
                f(satellite_mask.gnss_id, svid, signal_index) || return false
            end
        end
    end
    return true
end

"""
    parse_has_code_bias_block!(reader, context)

Parse the Code Biases block (ICD §5.2.5): one 11-bit bias per cell of the mask,
over [`e6b_foreach_cell`](@ref). Returns `nothing` on a truncated message.
"""
function parse_has_code_bias_block!(reader::HASBitReader, context::GalileoHASBlockContext)
    bits_remaining(reader) >= 4 || return nothing
    validity_interval = e6b_validity_interval(Int(read_bits!(reader, 4)))
    biases = GalileoHASCodeBias[]
    complete = e6b_foreach_cell(context.mask) do gnss_id, svid, signal_index
        bits_remaining(reader) >= 11 || return false
        raw = read_signed_bits!(reader, 11)
        push!(
            biases,
            GalileoHASCodeBias(;
                gnss_id,
                svid,
                signal_index,
                bias = raw == -1024 ? nothing : raw * 0.02,
            ),
        )
        return true
    end
    complete || return nothing
    return has_correction_block(context, validity_interval, biases)
end

"""
    parse_has_phase_bias_block!(reader, context)

Parse the Phase Biases block (ICD §5.2.6): per cell an 11-bit bias immediately
followed by its 2-bit Phase Discontinuity Indicator, over the same
[`e6b_foreach_cell`](@ref) walk the code biases use. Returns `nothing` on a
truncated message.
"""
function parse_has_phase_bias_block!(reader::HASBitReader, context::GalileoHASBlockContext)
    bits_remaining(reader) >= 4 || return nothing
    validity_interval = e6b_validity_interval(Int(read_bits!(reader, 4)))
    biases = GalileoHASPhaseBias[]
    complete = e6b_foreach_cell(context.mask) do gnss_id, svid, signal_index
        bits_remaining(reader) >= 13 || return false
        raw = read_signed_bits!(reader, 11)
        phase_discontinuity_indicator = Int(read_bits!(reader, 2))
        push!(
            biases,
            GalileoHASPhaseBias(;
                gnss_id,
                svid,
                signal_index,
                bias = raw == -1024 ? nothing : raw * 0.01,
                phase_discontinuity_indicator,
            ),
        )
        return true
    end
    complete || return nothing
    return has_correction_block(context, validity_interval, biases)
end

"""
    parse_has_message(octets, message_id, message_type, message_size, masks)
        -> Union{Nothing,Symbol,GalileoHASMessage}

Parse a reassembled HAS message. Returns

  - a [`GalileoHASMessage`](@ref) on success,
  - `:mask_unavailable` when the body needs a Mask ID that has not been received
    (the caller holds the message and retries once a matching mask arrives), or
  - `nothing` when the message is not parseable at all (reserved Message Type,
    out-of-range Time Of Hour, or a body that runs out of bits).

Blocks are read strictly in header-flag order (ICD §5.1). A block that fails to
parse ends the body — everything after it in the stream is unaligned — but the
blocks already parsed are kept, since each is self-contained.

The MT1 header's 4 reserved bits are skipped, per the ICD's instruction that
"forward-compatible receivers shall account for the possibility of HAS MT1
messages containing additional non-decodable information at the end of the
decodable one": trailing bits we cannot interpret are simply left unread.
"""
function parse_has_message(
    octets::Vector{UInt8},
    message_id::Int,
    message_type::Int,
    message_size::Int,
    masks::Union{Nothing,Dictionary{Int,GalileoHASMask}},
)
    # Only Message Type 1 is defined (ICD Table 10).
    message_type == 1 || return nothing
    reader = HASBitReader(octets)
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

    mask = nothing
    if mask_flag
        mask = parse_has_mask_block!(reader, mask_id)
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
        return GalileoHASMessage(;
            message_id,
            message_type,
            message_size,
            TOH,
            mask_id,
            IOD_set_id,
        )
    end
    context = GalileoHASBlockContext(mask, TOH, mask_id, IOD_set_id)

    # Content blocks appear in flag order (ICD Table 14). A block that fails to
    # parse leaves the reader unaligned, so everything after it is unreadable —
    # but the blocks already parsed are self-contained and are kept. Stating that
    # rule once, here, keeps it from having to be re-established at each of the
    # five call sites (and from being quietly omitted at the last one).
    aligned = true
    function parse_block(flag::Bool, parser)
        (aligned && flag) || return nothing
        block = parser(reader, context)
        aligned = !isnothing(block)
        return block
    end
    orbit_corrections = parse_block(orbit_flag, parse_has_orbit_block!)
    clock_corrections = parse_block(clock_full_set_flag, parse_has_clock_full_set_block!)
    clock_subset_corrections = parse_block(clock_subset_flag, parse_has_clock_subset_block!)
    code_biases = parse_block(code_bias_flag, parse_has_code_bias_block!)
    phase_biases = parse_block(phase_bias_flag, parse_has_phase_bias_block!)
    return GalileoHASMessage(;
        message_id,
        message_type,
        message_size,
        TOH,
        mask_id,
        IOD_set_id,
        mask = mask_flag ? mask : nothing,
        orbit_corrections,
        clock_corrections,
        clock_subset_corrections,
        code_biases,
        phase_biases,
    )
end

# ---- Page accumulation and message assembly (ICD §6.4) ----------------------

"""
    e6b_reassemble_message(group) -> Union{Nothing,Vector{UInt8}}

Run the HPVRS erasure decode on a complete page group and return the
`message_size × 53` message octets in page order, or `nothing` if the decode is
impossible.

The ICD's own recipe (§6.4): take the `k` rows of the generator matrix named by
the received Page IDs and its first `k` columns — the trailing `32 - k`
information symbols are known zeros, since "pages Ck+1 … CK contain only
zeroes" (§6.3) — invert that `k × k` matrix over GF(256) once, and apply it to
each of the 53 octet columns. GNSS-SDR instead runs a Berlekamp-Massey erasure
decode over the full 255-symbol codeword; both recover the same message, but the
matrix form is what the ICD specifies and needs no error-locator machinery.
"""
function e6b_reassemble_message(group::GalileoHASPageGroup)
    k = group.message_size
    received = view(group.octets, 1:k, :)
    decoded = rs_erasure_decode(
        GALILEO_HAS_GF256,
        E6B_GENERATOR_MATRIX,
        group.page_ids,
        received,
        k,
    )
    isnothing(decoded) && return nothing
    # Row `i` is non-encoded page `M_i`; the message is those pages concatenated,
    # so the k × 53 result is read out row by row.
    return vec(permutedims(decoded))
end

"""
Keep a block the new message carried, else the one already published. Unlike
`something` this tolerates both being `nothing` — the common case, since a
message carries only two or three of the six blocks.
"""
_e6b_latest(new, old) = isnothing(new) ? old : new

"""
Fold a decoded [`GalileoHASMessage`](@ref) into the accumulated data container:
publish it as the latest message, remember any mask it defined, and replace the
latest instance of each content block it carried.
"""
function e6b_merge_message(data::GalileoE6BData, message::GalileoHASMessage)
    masks = data.masks
    if !isnothing(message.mask)
        masks = _merge_keyed(masks, message.mask.mask_id, message.mask)
    end
    GalileoE6BData(
        data;
        message,
        masks,
        orbit_corrections = _e6b_latest(message.orbit_corrections, data.orbit_corrections),
        clock_corrections = _e6b_latest(message.clock_corrections, data.clock_corrections),
        clock_subset_corrections = _e6b_latest(
            message.clock_subset_corrections,
            data.clock_subset_corrections,
        ),
        code_biases = _e6b_latest(message.code_biases, data.code_biases),
        phase_biases = _e6b_latest(message.phase_biases, data.phase_biases),
    )
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
    stale = Int[]
    for (message_id, group) in pairs(cache.page_groups)
        now - group.opened_at >= E6B_MESSAGE_TIMEOUT_PAGES && push!(stale, message_id)
    end
    for message_id in stale
        unset!(cache.page_groups, message_id)
    end
    return cache
end

"""
Add one HAS Encoded Page to its Message ID's group, returning the group if it is
now complete (ready for the RS decode) and `nothing` otherwise.

A Message ID is reused over time, so a page whose Message Type or Message Size
disagrees with the group's starts a fresh group — the previous message is gone.
Duplicate Page IDs are ignored: the RS decode needs `k` *distinct* rows.
"""
function e6b_collect_page!(
    cache::GalileoE6BCache,
    message_id::Int,
    message_type::Int,
    message_size::Int,
    page_id::Int,
    octets::AbstractVector{UInt8},
)
    group = get(cache.page_groups, message_id, nothing)
    if isnothing(group) ||
       group.message_type != message_type ||
       group.message_size != message_size
        group = GalileoHASPageGroup(message_type, message_size, cache.page_counter[])
        set!(cache.page_groups, message_id, group)
    end
    page_id in group.page_ids && return nothing
    push!(group.page_ids, page_id)
    group.octets[length(group.page_ids), :] .= octets
    return length(group.page_ids) == group.message_size ? group : nothing
end

"""
    decode_syncro_sequence(state::GNSSDecoderState{<:GalileoE6BData}, sync::GalileoE6BSync)

Consume one CRC-validated C/NAV page: read its HAS Page Header, discard dummy
and unusable pages, add its HAS Encoded Page to the store for its Message ID,
and — once `MS` distinct pages are held — recover and parse the HAS message.

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
        # rather than merged into the existing one.
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

    encoded_page = Vector{UInt8}(undef, E6B_OCTETS_PER_PAGE)
    @inbounds for j = 1:E6B_OCTETS_PER_PAGE
        encoded_page[j] =
            UInt8(get_bits(page, E6B_PAGE_BITS, E6B_ENCODED_PAGE_START + 8 * (j - 1), 8))
    end

    group = e6b_collect_page!(
        cache,
        message_id,
        message_type,
        message_size,
        page_id,
        encoded_page,
    )
    isnothing(group) && return state

    octets = e6b_reassemble_message(group)
    unset!(cache.page_groups, message_id)
    isnothing(octets) && return state

    return e6b_apply_message(state, message_id, message_type, message_size, octets)
end

"""
Parse a reassembled message and fold it into `raw_data`, handling the
mask-not-yet-received case by holding the message for a later retry, and
retrying any previously held message once a new mask has been learnt.
"""
function e6b_apply_message(
    state::GNSSDecoderState{<:GalileoE6BData},
    message_id::Int,
    message_type::Int,
    message_size::Int,
    octets::Vector{UInt8},
)
    cache = state.cache
    parsed = parse_has_message(
        octets,
        message_id,
        message_type,
        message_size,
        state.raw_data.masks,
    )
    if parsed === :mask_unavailable
        # Hold the message; HAS broadcasts a defining mask every few messages,
        # and re-parsing then recovers corrections that would otherwise be lost.
        # One slot is enough: only the most recent orphan is worth keeping.
        cache.pending_message[] = GalileoHASPendingMessage(
            message_id,
            message_type,
            message_size,
            cache.page_counter[],
            octets,
        )
        return state
    end
    isnothing(parsed) && return state
    raw = state.raw_data
    # A message that carried a *new* mask may unlock the held-back one; a message
    # that merely referenced an existing one cannot, so only check when
    # `parsed.mask` is present.
    #
    # ORDER MATTERS. The held message completed *earlier* than `parsed`, so it is
    # merged *first* and `parsed` folded over the top: `e6b_merge_message` lets
    # the later of two messages win every field it carries, so merging the orphan
    # last would publish it as `message` (documented as the most recently
    # completed one) and would replace a fresher correction block with a staler
    # one of the same kind. Retrying needs the mask `parsed` carries but not the
    # rest of it, so the mask is merged into the dictionary on its own for the
    # retry.
    held = cache.pending_message[]
    if !isnothing(parsed.mask) && !isnothing(held)
        retried = parse_has_message(
            held.octets,
            held.message_id,
            held.message_type,
            held.message_size,
            _merge_keyed(raw.masks, parsed.mask.mask_id, parsed.mask),
        )
        if retried isa GalileoHASMessage
            cache.pending_message[] = nothing
            raw = e6b_merge_message(raw, retried)
        end
    end
    return GNSSDecoderState(state; raw_data = e6b_merge_message(raw, parsed))
end

"""
    validate_data(state::GNSSDecoderState{<:GalileoE6BData})

Publish `raw_data` as `data` unconditionally.

Every other decoder here withholds `data` until a cross-message consistency
check passes — matching issues of data, a repeated broadcast, a plausible time
of week. C/NAV has no such check to make and needs none: a HAS message is
reassembled only from pages that each passed CRC-24Q, under one Message ID with
one Message Size, and the ICD's whole erasure-channel model rests on that gate
(§6.4). There is nothing left to corroborate, and withholding corrections that
carry their own validity intervals would only make them stale.
"""
function validate_data(state::GNSSDecoderState{<:GalileoE6BData})
    GNSSDecoderState(state; data = state.raw_data)
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
