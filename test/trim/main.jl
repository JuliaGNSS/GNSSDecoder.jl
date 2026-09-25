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

function run_decoder(name, system, prn, symbols)
    # Both entry points go into the trimmed image: the overwriting `decode!` and
    # `decode`, which decodes into a copy of its argument. Each gets a fresh
    # state — a replay of the same stream on a decoded state would (rightly)
    # fail the time-of-week plausibility screens.
    state = decode!(GNSSDecoderState(system, prn), symbols, length(symbols))
    copied = decode(GNSSDecoderState(system, prn), symbols, length(symbols))
    reset_decoder_state!(copied)
    println(
        Core.stdout,
        name,
        " PRN ",
        prn,
        ": healthy=",
        is_sat_healthy(state),
        " complete=",
        is_decoding_completed_for_positioning(state),
    )
    return state
end

function (@main)(args::Vector{String})::Cint
    symbols = pseudo_random_symbols(8000)
    gps_l1ca =
        GNSSDecoder._precompile_soft_symbols(GNSSDecoder._PRECOMPILE_GPS_L1CA_SYMBOLS)
    galileo_e1b =
        GNSSDecoder._precompile_soft_symbols(GNSSDecoder._PRECOMPILE_GALILEO_E1B_SYMBOLS)

    l1ca = run_decoder("GPS L1 C/A", GPSL1CA(), 25, gps_l1ca)
    run_decoder("GPS L1C-D", GPSL1C_D(), 7, symbols)
    run_decoder("GPS L2C", GPSL2CM(), 25, symbols)
    run_decoder("GPS L5I", GPSL5I(), 25, symbols)
    run_decoder("Galileo E1B", GalileoE1B(), 2, galileo_e1b)
    run_decoder("Galileo E5a", GalileoE5aI(), 21, symbols)
    run_decoder("Galileo E5b", GalileoE5bI(), 2, galileo_e1b)
    run_decoder("Galileo E6B", GalileoE6B(), 2, symbols)
    run_decoder("BeiDou B1I", BeiDouB1I(), 7, symbols)
    run_decoder("BeiDou B1C", BeiDouB1C_D(), 7, symbols)
    run_decoder("BeiDou B2a", BeiDouB2aI(), 7, symbols)
    run_decoder("BeiDou B2b", BeiDouB2bI(), 7, symbols)
    run_decoder("BeiDou B3I", BeiDouB3I(), 7, symbols)

    if !is_decoding_completed_for_positioning(l1ca) || get_time_of_week(l1ca) != 259278
        println(Core.stderr, "GPS L1 C/A capture did not decode to its navigation fix")
        return 1
    end
    return 0
end
