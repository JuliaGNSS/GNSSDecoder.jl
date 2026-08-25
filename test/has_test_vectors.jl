# Galileo HAS (E6-B C/NAV) reference test vectors, transcribed from the official
# Galileo High Accuracy Service Signal-In-Space Interface Control Document
# (HAS SIS ICD), Issue 1.0 — the version in force at the European GNSS Service
# Centre (https://www.gsc-europa.eu/electronic-library/programme-reference-documents,
# "Galileo HAS in force Reference Documents"). © European Union 2022.
#
# Reproduced here, unmodified, solely so this package's decoder can be verified
# against the specification's own worked examples — the interoperability use the
# ICD is published for, and the same thing GNSS-SDR does with this example in
# `tests/unit-tests/system-parameters/has_decoding_test.cc`. The document's Terms
# of Use and Disclaimers (ICD pages i-ii) apply to this data.
#
# Everything is stored as MSB-first hexadecimal strings, one per page or matrix
# row; the ICD prints the encoded pages in decimal and the decoded messages in
# hex, so only the former were converted. Four independent things are captured:
#
#  1. `HAS_ANNEX_C_PAGE_HEX` — Annex C's sample C/NAV page: the 486 information
#     bits as broadcast (zero-padded to 512), with the header fields and the 53
#     HAS Encoded Page octets the ICD decodes from it. Pins down the page bit
#     layout (Reserved 14 | Header 24 | Encoded Page 424 | CRC 24), the CRC-24Q
#     scope, and the HAS Page Header field order — the one part of the format
#     where the ICD's rotated-text table is easy to misread.
#  2. `HAS_ANNEX_B_ROWS_HEX` — four rows of the Annex B RS(255, 32) systematic
#     generator matrix, as spot checks that this package *generates* the matrix
#     the ICD tabulates. Rows 33 and 34 are the first two parity rows, 254 and
#     255 the last two.
#  3. `HAS_GENERATOR_POLYNOMIAL_HEAD` / `_TAIL` — coefficients of the RS
#     generator polynomial from ICD Table 42.
#  4. `HAS_EXAMPLE_1` / `HAS_EXAMPLE_2` — Annex D's two complete HPVRS decoding
#     examples: the received encoded pages with their Page IDs, and the message
#     octets the ICD states they reassemble to. The expected *parsed* field
#     values are asserted in `test/galileo_e6b.jl` against the same annex.
#
# The Annex B and Annex D data ship as attachments inside the ICD PDF
# (`Galileo-HAS-SIS-ICD_1.0_Annex_B_Reed_Solomon_Generator_Matrix.txt` and
# `Galileo-HAS-SIS-ICD-1.0_Annex_D_HAS_Message_Decoding_Example.txt`); extract
# them with any tool that reads PDF embedded files (`pdfdetach -saveall`) to
# re-derive this file.

"""
Concatenate MSB-first hex strings into one octet vector.
"""
has_octets(hex_strings...) = mapreduce(hex2bytes, vcat, hex_strings)

# ---- Annex C: one sample C/NAV page -----------------------------------------

"""
486-bit C/NAV page (Annex C), MSB first, zero-padded to 512 bits.
"""
const HAS_ANNEX_C_PAGE_HEX = "fffc17b8de11ef1d27adf5c5d0911e23ed151a4630009cabaf05524b31bad56962038986eb8c5c688f742f958bf235bf623988a70a79f632677d0c4690000000"

"""
Header fields the ICD decodes from `HAS_ANNEX_C_PAGE_HEX` (Annex C table).
"""
const HAS_ANNEX_C_HEADER =
    (HAS_status = 0, message_type = 1, message_id = 15, message_size = 15, page_id = 55)

"""
The 53 HAS Encoded Page octets the ICD lists for `HAS_ANNEX_C_PAGE_HEX`
(identical to Page ID 55 of `HAS_EXAMPLE_1`).
"""
const HAS_ANNEX_C_ENCODED_PAGE_HEX = "847bc749eb7d7174244788fb4546918c00272aebc15492cc6eb55a5880e261bae3171a23dd0be562fc8d6fd88e6229c29e7d8c99df"

# ---- Annex B: generator-matrix spot checks ----------------------------------

"""
Rows 33, 34, 254 and 255 of the Annex B RS(255, 32) systematic generator matrix.
"""
const HAS_ANNEX_B_ROWS_HEX = Dict(
    33 => "138fb43bdd1d312de70949499f029e88d4da0e71d714bb3789b5cb7161870efb",
    34 => "1b1b0132ff6dfb9c949755154a74fa4d3ccb71c4d517ca7d1ffc5a01b0e22cfc",
    254 => "f478562a6ecbd19e7773cf05688c8a7119993bab694388461e0acb500dc8acd8",
    255 => "744034ae367e10c2a221219db0c5e10c3b37fde4942fb3b9188afd148e37ac58",
)

"""
First 14 coefficients g_0 … g_13 of the RS generator polynomial (ICD Table 42).
"""
const HAS_GENERATOR_POLYNOMIAL_HEAD =
    UInt8[88, 216, 195, 23, 111, 82, 79, 81, 62, 120, 249, 250, 11, 134]

"""
Leading coefficient g_223 of the RS generator polynomial — it is monic (ICD Table 42).
"""
const HAS_GENERATOR_POLYNOMIAL_TAIL = UInt8(1)

# ---- Annex D: two complete HPVRS decoding examples --------------------------

"""
Annex D example 1: a Message Type 1 of `MS = 15` pages carrying a mask,
orbit corrections, code biases and phase biases (flag sequence "110011"), with
the 15 received encoded pages and the message they reassemble to.
"""
const HAS_EXAMPLE_1 = (
    message_id = 15,
    message_type = 1,
    message_size = 15,
    page_ids = Int[55, 56, 57, 58, 59, 174, 175, 176, 187, 188, 239, 240, 241, 252, 253],
    encoded_pages_hex = [
        "847bc749eb7d7174244788fb4546918c00272aebc15492cc6eb55a5880e261bae3171a23dd0be562fc8d6fd88e6229c29e7d8c99df",   # Page ID 55
        "349ae3634d210bad3293a67fb62101e9dd54307bc679ed699bd50caeaec56485f3f816540caecea4c61692ee5b18caabb5bda27939",   # Page ID 56
        "55011d910ee6e155c2f28c4dd7fad628c8e26a05abd787974de2e16f8ef6b09c00d712e42908229718aeec691c0527f3c23f80b513",   # Page ID 57
        "2ca31b231553ee6a9c7a3bfffa842b2d0cf3080910b9c2027e8873dced2f8da7d423a42fd9ce58c3ee447d2caf31b18a04d5a5ba78",   # Page ID 58
        "37be60d823798db61a1c9822eef84b7ad5ed63d5223d98ad91cc858f4075775ce04cbb24a0d0b15f7fd53ad6862c79f8523fa9bf4b",   # Page ID 59
        "bb1c451d5904a0e416b92b589a0c56ce2bc7739828ef0bc049e491189a293f312824e0b0645e1f64986d6f87b976cf3a12f73b9021",   # Page ID 174
        "7519489afbc26f45cabffd9f78b2f644ab29fba37ccafeef98190205ccdfc0e7fa78c1b3ea506ca6a6a7d2c363879f76848fa48024",   # Page ID 175
        "8f0c9c348bcbc13d590335540ea865c2cf3d713bbc27c8631a2958ded386b275470f88969641587ccc80171c33a6ccddfb3f352cbe",   # Page ID 176
        "cbe2240a911b3681f38e2b3ff239f362e53b4ac9292c60c77c61c546764e86426a8a44c5408cbb5bc90a8a8710fe6d7190dc80cc5d",   # Page ID 187
        "1d379ea7c3df909e9e7457db6524471cbd34d711c75cb08b4a846c03197e2ebfe2ef0ea12c46f7fdcaf63a24231d4d90340ed98bdd",   # Page ID 188
        "7a392815304163154d32cc1ee9a67503300373fae04e8f6cf590ffc79372a12691296bac84525fcaa6984b53588f1919baca979fde",   # Page ID 239
        "7d1338cf705cb893efb571d118f5ad39ad3303a094ffb65c8ca892c2ea3d35be890f5be4e7096fde343ecdbd5ab981de4a139a5e1d",   # Page ID 240
        "a1cc75defd3dc942cf6a15a67595e0a4f9322dac47cd1d577051b15fd782d6a2532bb609bc70b76f05aee7b067977507e8a71321ea",   # Page ID 241
        "cf93cd158cf41fb295ad9d21a1558282ed74883336896a7b7eead039912274e5d1e21a563feff5d215d33dbd2b55d767a0aaeaa338",   # Page ID 252
        "d7c8a713d2a61260e04d05916a94de679dc4e9846d3de5bba398113e1bd22a43b502176c44cebd4c3a27a42bfe09572912e487d4a5",   # Page ID 253
    ],
    message_octets_hex = [
        "000cc00b20ffdfffff008100f7ffff7df55ffdfe0beee8a79a41241000a6000a01a01280400200200113fbc041febbf00080080042",   # non-encoded page M1
        "ff6822fea21807c193f7598035fd7f6a2f00080080016ff90287e7967f702580587fee217a10c9dfcc0e7f651df577d981603ffe41",   # non-encoded page M2
        "47f903ff9df7805c15ff9fdcff8008004004000a002407ff9d7c07df7ffe2b5fdcee305519011fd7fd24479f00500e8e7edc31401c",   # non-encoded page M3
        "43fdb02304007fe5030ff1ac40020020000200100100077fec06e00141feb02afcb2c400200200043ff5f6c022097f7c0e3f4412ff",   # non-encoded page M4
        "4fe1ff8825fe8ffcff0048081fe3fda097f4c04bf3812fe5ff27f0025fc6ff5ff40480edfa601c08ffe8023fcc0f00b00b80a825fd",   # non-encoded page M5
        "f00fff704bf71ffffdc097fb400c00812fe781a7f8025fe602203204801001a01607ffd006404012fec00e000825fc7fe500c04bff",   # non-encoded page M6
        "405605c08804004403012fe27feffbf0bb23dc94458ef0420afe1fa61544abda77c130444320a1104303d3f76f65fbbee7ccf5fe6b",   # non-encoded page M7
        "ddf8bfcff479b7a5f1dc3bf3fce1243b44e90d1784ac350b2f29f2bd607b1a1e7bb207519201003807069f8feb7cf00c0d42d85b06",   # non-encoded page M8
        "1f33d2fa7fa00fc3506a02015c4b09409bf07cbf950400641582a04fc8f40e88d2dd9f73efbdc40080400407c198588ad0e9f43d67",   # non-encoded page M9
        "aef9009c220420cdefbc9f90f920f0338660401a45a0b411a0841c8380c206c1882d0121243e87d02bf27d1fa2fc6184518a50dcb0",   # non-encoded page M10
        "0080040020010008004002001000800400200100080040020010008004002001000800400200100080040020010008004002001000",   # non-encoded page M11
        "8004002001000800400200100080040020010008004002001000800400200100080040020010008004002001000800400200100080",   # non-encoded page M12
        "0400200100080040020010008004002001000800400200100080040020010008004002001000800400200100080040020010008004",   # non-encoded page M13
        "0020010008004002001000800400200100080040020010008004002001000800400200100080040020010008004002001000800400",   # non-encoded page M14
        "2001000800400200100080040020010008002aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",   # non-encoded page M15
    ],
)

"""
Annex D example 2: a Message Type 1 of `MS = 2` pages carrying only clock
full-set corrections (flag sequence "001000"), which therefore references the
mask of another message by Mask ID. Exercises the smallest possible RS erasure
decode (k = 2) and the Delta Clock Multiplier semantics. The ICD does not state
this example's Message ID; 0 is used here.
"""
const HAS_EXAMPLE_2 = (
    message_id = 0,
    message_type = 1,
    message_size = 2,
    page_ids = Int[61, 151],
    encoded_pages_hex = [
        "cb56b07e135478006537f8cc42c9803df5cbc4313a55781a8d6f818d1129b5db164cf7a9a49f77e5e106e54b33a55c8bcd74fbfd94",   # Page ID 61
        "f8f0ecfe02681e00b071d818a911e134a2f84c56d3f59a8a8025f1ecdbbc9887df0d33f3ec39f17483de74e5cdbfe8a56cf07e7c5b",   # Page ID 151
    ],
    message_octets_hex = [
        "0072000b58afe4002d03000acd5826ae3000aaa5532b15581aaa572aa175b8800516e941454a28550ebd5556aa8c002001546a92c0",   # non-encoded page M1
        "02c08020fd6ff200bbfe4fe2fec41020210207ff7f85ff8007002bfe202d000ffbc052044febaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",   # non-encoded page M2
    ],
)
