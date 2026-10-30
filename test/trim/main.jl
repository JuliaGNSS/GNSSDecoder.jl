# Trim smoke test: an executable that runs every decoder, compiled with
#
#     juliac --output-exe gnssdecoder_trim --bundle build --trim=safe --experimental --project test/trim test/trim/main.jl
#
# `--trim=safe` fails the build on any call it cannot resolve statically, so a
# successful build is the actual check; running the binary then confirms the
# trimmed image decodes. It exits non-zero if the GPS L1 C/A capture (the one
# stream here that carries a full navigation fix) no longer decodes to it.
using GNSSDecoder, GNSSSignals

# Deterministic ±1 symbols from a xorshift generator, so the driver needs no
# Random and its output is reproducible. They exercise each decoder's sync
# search, FEC and CRC paths.
function pseudo_random_symbols(n)
    state = 0x9e3779b97f4a7c15
    out = Vector{Float32}(undef, n)
    for i = 1:n
        state ⊻= state << 13
        state ⊻= state >> 7
        state ⊻= state << 17
        out[i] = isodd(state) ? 1.0f0 : -1.0f0
    end
    out
end

function run_decoder(name, make_decoder, prn, symbols)
    # `copy` and `reset_decoder_state!` go into the trimmed image as well: the
    # copy of a decoded state is reset, then decodes the stream afresh (a
    # replay without the reset would, rightly, fail the time-of-week
    # plausibility screens).
    state = decode!(make_decoder(prn), symbols)
    decode!(reset_decoder_state!(copy(state)), symbols, length(symbols))
    println(
        Core.stdout,
        name,
        " PRN ",
        prn,
        ": healthy=",
        is_sat_healthy(state),
        " complete=",
        has_positioning_data(state),
    )
    return state
end

function (@main)(args::Vector{String})::Cint
    # Construct decoders without allocating spreading-code lookup tables.
    # E1B below also exercises generic construction from both tracking replicas.
    symbols = pseudo_random_symbols(8000)
    gps_l1ca =
        GNSSDecoder._precompile_soft_symbols(GNSSDecoder._PRECOMPILE_GPS_L1CA_SYMBOLS)
    galileo_e1b =
        GNSSDecoder._precompile_soft_symbols(GNSSDecoder._PRECOMPILE_GALILEO_E1B_SYMBOLS)

    l1ca = run_decoder("GPS L1 C/A", GPSL1CADecoderState, 25, gps_l1ca)
    run_decoder("GPS L1C-D", GPSL1C_DDecoderState, 7, symbols)
    run_decoder("GPS L2C", GPSL2CMDecoderState, 25, symbols)
    run_decoder("GPS L5I", GPSL5IDecoderState, 25, symbols)
    e1b = run_decoder("Galileo E1B", GalileoE1BDecoderState, 2, galileo_e1b)
    run_decoder("Galileo E5a", GalileoE5aDecoderState, 21, symbols)
    run_decoder("Galileo E5b", GalileoE5bDecoderState, 2, galileo_e1b)
    run_decoder("Galileo E6B", GalileoE6BDecoderState, 2, symbols)
    run_decoder("BeiDou B1I", BeiDouB1IDecoderState, 7, symbols)
    run_decoder("BeiDou B1C", BeiDouB1CDecoderState, 7, symbols)
    run_decoder("BeiDou B2a", BeiDouB2aDecoderState, 7, symbols)
    run_decoder("BeiDou B2b", BeiDouB2bDecoderState, 7, symbols)
    run_decoder("BeiDou B3I", BeiDouB3IDecoderState, 7, symbols)

    # Replica construction must also trim safely. Both tracking variants feed
    # the same decoder and retain only the broadcast component identity.
    for signal in (GalileoE1B(), GalileoE1B_CBOC())
        decoded = decode!(GNSSDecoderState(signal, 2), galileo_e1b)
        if decoded != e1b ||
           get_signal_id(decoded) !== :GalileoE1B ||
           get_signal_name(decoded) != "Galileo E1B"
            println(Core.stderr, "Galileo E1B replica changed decoder data or identity")
            return 1
        end
    end

    if !has_positioning_data(l1ca) || get_broadcast_time_of_week(l1ca) != 259278
        println(Core.stderr, "GPS L1 C/A capture did not decode to its navigation fix")
        return 1
    end
    if isnothing(get_decoded_time_of_week(l1ca))
        println(Core.stderr, "GPS L1 C/A decoded time-of-week metadata is missing")
        return 1
    end
    if isnothing(get_broadcast_ephemeris_validity(l1ca; approximate_year = 2020))
        println(Core.stderr, "GPS L1 C/A ephemeris validity metadata is missing")
        return 1
    end
    if !isnothing(get_broadcast_time_offset_parameters(l1ca, GST()))
        println(Core.stderr, "GPS L1 C/A unexpectedly has a broadcast GST offset")
        return 1
    end
    return 0
end
