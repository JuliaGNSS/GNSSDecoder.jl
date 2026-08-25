# Opt-in real-data integration test.
#
# Decodes the Fraunhofer IIS Flexiband "III-7a" reference capture — the same
# real-sky recording (Hanoi, 2017-10-17, CC BY-NC) that Tracking.jl uses for its
# own integration test — end to end: Acquisition.jl + Tracking.jl produce soft
# symbols, which `GNSSDecoder.decode` turns into navigation data. Five decoders
# run against it: GPS L1 C/A, Galileo E1B and BeiDou B1I on the L1 front-end,
# GPS L5I and Galileo E5a-I on the L5 one. This is the real-signal counterpart
# to the synthetic / golden-fixture unit tests, and it exercises the multi-band
# `track!` path (two front-ends at different rates, three band keys).
#
# WHICH SIGNALS THIS CAPTURE CAN AND CANNOT REACH. The three sampled bands are
# L1 at 1580.0 ± 40.5 MHz, L5 at 1173.546875 ± 20.25 MHz and S at 2492.028 ±
# 10.125 MHz, so by carrier frequency alone the capture excludes GPS L2C
# (1227.60), Galileo E5b-I and BeiDou B2b (both 1207.14), Galileo E6-B (1278.75)
# and BeiDou B3I (1268.52) — those decoders cannot be exercised here at any
# epoch. The capture date rules out three more that *are* in band: GPS L1C-D
# (no satellite broadcast L1C until GPS III SV01, December 2018) and BeiDou B1C
# / B2a (BDS-3 MEO deployment began 2017-11-05, three weeks after this
# recording). That leaves exactly the five decoded below.
#
# Capture format (multiplexed `.usb`, see Tracking.jl's flexiband test): three
# bands interleaved in 6-byte base-clock (20.25 MHz) chunks, in lump order:
#
#     byte:  1   |  2     3     4     5   |    6
#     band:  S   | L1_0  L1_1  L1_2  L1_3 | L5 (2 samples)
#     rate:  1×  |          4×            |    2×
#
#   • S  (IRNSS, 20.25 MHz) — 1 sample/chunk. Not used here.
#   • L1 (81 MHz, GPS L1 C/A + Galileo E1) — 4 samples/chunk, 4-bit signed I/Q
#     packed I=high nibble, Q=low nibble, two's complement.
#     IF = 1575.42 - 1580.0 = -4.58 MHz.
#   • L5 (40.5 MHz, GPS L5) — 2 samples/chunk packed into one byte: 2-bit
#     signed I/Q, sample 0 in the high bits, sign-magnitude levels
#     {00,01,10,11} → {+1,+3,-3,-1}. IF = 1176.45 - 1173.546875 = +2.903125 MHz.
#
# The 81:40.5 MHz rate ratio is exactly 2:1, so the L1 and L5 buffers (filled
# from the same frames) span the same wall clock sample-for-sample, as the
# multi-band `track!` boundary check requires.
#
# WHY WE ASSERT ON `raw_data`, NOT `data`:
#   `data` is only published once a *complete, cross-validated* positioning set
#   is assembled — for GPS L1 C/A that needs subframes 1+2+3 (~18 s) with
#   consistent IODC, Galileo E1B needs its full word-type set (~30 s), and GPS
#   L5I needs message types 10+11+a clock message (~18 s). The `III-7a_short`
#   capture is only ~15.7 s long (its uncompressed `.usb` is ~1.93 GB; despite
#   the website calling it "one minute", the actual sample stream is ~16 s),
#   which is below that span. `raw_data`, however, is populated from
#   individually parity-/CRC-validated subframes/words/messages (GPS L1 C/A HOW
#   TOW gated by the 6-equation parity check; E1B TOW/WN gated by CRC-24Q over
#   the page pair; L5I TOW gated by CRC-24Q over the 300-bit message). So a
#   decoded `raw_data.TOW` is a genuine, error-checked end-to-end decode —
#   exactly what this test verifies — even though the capture is too short to
#   publish `data`. (Use the L125 210 s captures for a full-ephemeris `data`
#   assertion.)
#
# OPT-IN: skipped unless ENV["GNSSDECODER_RUN_INTEGRATION_TEST"] == "true", so a
# plain `]test` never downloads. CI enables it on one job and caches the capture
# (set GNSSDECODER_TESTDATA_DIR to the cached directory).

if get(ENV, "GNSSDECODER_RUN_INTEGRATION_TEST", "false") != "true"
    @info "Skipping Flexiband III-7a integration test (downloads a 1.6 GB capture). " *
          "Set ENV[\"GNSSDECODER_RUN_INTEGRATION_TEST\"] = \"true\" to run it."
else
    using Downloads: Downloads
    using ZipFile: ZipFile
    using Unitful: Hz, ustrip
    using GNSSSignals: GPSL1CA, GPSL5I, GalileoE1B, GalileoE5aI, BeiDouB1I
    using Acquisition: acquire, is_detected
    using Tracking: TrackState, BandMeasurement, add_satellite!, track!, get_soft_bits

    const III7A_URL =
        "https://www.iis.fraunhofer.de/content/dam/iis/de/doc/lv/los/lokalisierung/" *
        "SatNAV/Flexiband%20reference%20Data%2020171017_09-51-43_III-7a_short.zip"
    const III7A_SIZE = 1_589_435_475       # exact byte count, for integrity
    const III7A_FS = 81.0e6Hz              # L1 band sampling rate
    const III7A_FS_L5 = 40.5e6Hz           # L5 band sampling rate
    const III7A_IF = -4.58e6Hz             # GPS L1 / Galileo E1 within the 1580 MHz band
    const III7A_IF_L5 = 2.903125e6Hz       # GPS L5 / Galileo E5a in the 1173.546875 MHz band
    # BeiDou B1I sits 18.9 MHz below the L1 front-end's 1580 MHz centre: inside
    # the 81 MHz sampled span, but far enough out that the analog filter, not the
    # sample rate, decides whether anything survives there. Hence the soft
    # handling below — this band edge is the one part of the capture we cannot
    # predict from the format alone.
    const III7A_IF_B1I = 1561.098e6Hz - 1580.0e6Hz

    # `.usb` framing (Flexiband / ION SDR-metadata format).
    const III7A_BLOCK = 1024;
    const III7A_HEADER = 6;
    const III7A_CHUNK = 6
    const III7A_CYCLES = 169;
    const III7A_PAYLOAD = III7A_CYCLES * III7A_CHUNK
    const III7A_FREQBASE = 20.25e6
    const III7A_L5_LEVELS = (1.0f0, 3.0f0, -3.0f0, -1.0f0)   # 2-bit sign-magnitude
    iii7a_nframes(sec) = ceil(Int, sec * III7A_FREQBASE / III7A_CYCLES)

    # Download the capture in verified byte-range chunks (the Fraunhofer server
    # drops long-running connections), assembling only once every chunk is exact.
    function iii7a_fetch()
        dir = get(ENV, "GNSSDECODER_TESTDATA_DIR", tempdir());
        mkpath(dir)
        path = joinpath(dir, "Flexiband_III-7a_short.zip")
        (isfile(path) && filesize(path) == III7A_SIZE) && return path
        chunk = 200_000_000;
        parts = String[];
        start = 0;
        idx = 0
        tmpdir = mktempdir(dir)
        try
            while start < III7A_SIZE
                stop = min(start + chunk - 1, III7A_SIZE - 1);
                want = stop - start + 1
                part = joinpath(tmpdir, "p$(idx)");
                tries = 0
                while !(isfile(part) && filesize(part) == want)
                    (tries += 1) > 30 && error("III-7a download: chunk $idx failed")
                    try
                        Downloads.download(
                            III7A_URL,
                            part;
                            headers = ["Range" => "bytes=$start-$stop"],
                        )
                    catch err
                        @warn "III-7a chunk $idx retry $tries" err
                    end
                end
                push!(parts, part);
                start = stop + 1;
                idx += 1
            end
            open(path * ".part", "w") do out
                for p in parts
                    ;
                    write(out, read(p));
                end
            end
            filesize(path * ".part") == III7A_SIZE ||
                error("III-7a assembled size mismatch")
            mv(path * ".part", path; force = true)
        finally
            rm(tmpdir; recursive = true, force = true)
        end
        path
    end

    iii7a_open(zippath) = (
        r = ZipFile.Reader(zippath);
        (r, only(filter(f -> endswith(f.name, ".usb"), r.files)))
    )

    # Demux the L1 (81 MHz; GPS L1 C/A + Galileo E1) and L5 (40.5 MHz; GPS L5)
    # bands from the next `nblocks` USB frames. Both are filled from the same
    # frames, so the returned buffers are sample-for-sample time-aligned with
    # length(l1) == 2 * length(l5) (the exact 81:40.5 = 2:1 ratio the multi-band
    # `track!` boundary check needs). L1 is 4-bit two's-complement I/Q (I = high
    # nibble); the L5 byte packs two 2-bit sign-magnitude I/Q samples, sample 0
    # in the MSBs.
    #
    # Reading is windowed (one call per `track!` chunk) and demuxes straight from
    # the decompressing zip stream, so only one window of samples is ever live —
    # the full ~1.93 GB `.usb` is never materialized. The stream also
    # self-terminates at the first malformed frame: only the first ~13.68 s of
    # this capture is well-formed; the trailing ~2 s has broken `0x55 0xAA`
    # preambles and demuxes to noise (cf. Tracking.jl #157/#158). We stop at the
    # first bad preamble and return only the valid prefix rather than feeding
    # garbage tail frames to `track!`/`decode`.
    function iii7a_read_bands!(usb, nblocks)
        frame = Vector{UInt8}(undef, III7A_BLOCK)
        l1 = Vector{ComplexF32}(undef, III7A_CYCLES * nblocks * 4)
        l5 = Vector{ComplexF32}(undef, III7A_CYCLES * nblocks * 2)
        got = 0
        nib(x) = (v = Int(x & 0x0f); v ≥ 8 ? v - 16 : v)
        lvl(x) = @inbounds III7A_L5_LEVELS[(x&0x03)+1]
        for _ = 1:nblocks
            eof(usb) && break
            read!(usb, frame)
            (frame[1] == 0x55 && frame[2] == 0xAA) || break
            got += 1
            @inbounds for c = 0:(III7A_CYCLES-1)
                base = III7A_HEADER + c * III7A_CHUNK
                for k = 0:3
                    b = frame[base+2+k]
                    l1[(got-1)*III7A_CYCLES*4+c*4+k+1] = ComplexF32(nib(b >> 4), nib(b))
                end
                b5 = frame[base+6]
                l5[(got-1)*III7A_CYCLES*2+c*2+1] = ComplexF32(lvl(b5 >> 6), lvl(b5 >> 4))
                l5[(got-1)*III7A_CYCLES*2+c*2+2] = ComplexF32(lvl(b5 >> 2), lvl(b5))
            end
        end
        resize!(l1, got * III7A_CYCLES * 4), resize!(l5, got * III7A_CYCLES * 2)
    end

    @testset "Flexiband III-7a real-data decode (L1: GPS C/A + Galileo E1B + BeiDou B1I; L5: GPS L5I + Galileo E5a-I)" begin
        fshz = ustrip(Hz, III7A_FS);
        fshz_l5 = ustrip(Hz, III7A_FS_L5);
        nblk(sec) = iii7a_nframes(sec)
        zippath = iii7a_fetch()
        @test filesize(zippath) == III7A_SIZE

        # Acquire on the first 40 ms of each band.
        r, usb = iii7a_open(zippath);
        l1h, l5h = iii7a_read_bands!(usb, nblk(0.045));
        close(r)
        @test length(l1h) == 2 * length(l5h)
        nacq = round(Int, fshz * 0.040);
        nacq_l5 = round(Int, fshz_l5 * 0.040)
        det(a) = filter(x -> is_detected(x; pfa = 1e-8), a)
        topN(a, n) = sort(a; by = x -> -x.CN0)[1:min(n, length(a))]
        acq_gps = topN(
            det(
                acquire(
                    GPSL1CA(),
                    (@view l1h[1:nacq]),
                    III7A_FS,
                    1:32;
                    interm_freq = III7A_IF,
                    num_coherently_integrated_code_periods = 4,
                    num_noncoherent_accumulations = 2,
                ),
            ),
            4,
        )
        acq_gal = det(
            acquire(
                GalileoE1B(),
                (@view l1h[1:nacq]),
                III7A_FS,
                1:36;
                interm_freq = III7A_IF,
                num_coherently_integrated_code_periods = 1,
                num_noncoherent_accumulations = 4,
            ),
        )
        # GPS L5 is only carried by Block IIF satellites (Oct 2017), a handful in
        # view: integrate one full 10 ms NH10 secondary-code period (the longest
        # coherent window before the data symbol can flip).
        acq_l5 = det(
            acquire(
                GPSL5I(),
                (@view l5h[1:nacq_l5]),
                III7A_FS_L5,
                1:32;
                interm_freq = III7A_IF_L5,
                num_coherently_integrated_code_periods = 10,
                num_noncoherent_accumulations = 1,
            ),
        )
        # Galileo E5a-I shares the L5 front-end with GPS L5I — same carrier, same
        # buffer, same IF. It is *not* integrated the way the L5I line above is,
        # though. Acquisition.jl requires the coherent length to be a whole
        # multiple of the secondary code (CS20 here, NH10 there) so it can search
        # the code's rotation, and that search resolves for NH10 but never
        # declares a detection for either 20-chip secondary code — checked on
        # synthetic signals up to ~71 dB-Hz, where the CN0 estimate is correct and
        # `is_detected` still returns false, for E5a-I and BeiDou B1I alike. So
        # both take the fully non-coherent route instead: 20 single-period
        # accumulations, which detects at ~53 dB-Hz synthetically and needs no
        # secondary-phase hypothesis at all.
        acq_e5a = topN(
            det(
                acquire(
                    GalileoE5aI(),
                    (@view l5h[1:nacq_l5]),
                    III7A_FS_L5,
                    1:36;
                    interm_freq = III7A_IF_L5,
                    num_coherently_integrated_code_periods = 1,
                    num_noncoherent_accumulations = 20,
                ),
            ),
            4,
        )
        # BeiDou B1I, from the same L1 samples at its own IF, non-coherently for
        # the reason just given. PRNs 6+ only: 1-5 are the BDS-2 GEOs, which
        # broadcast D2 NAV *without* the NH20 secondary code that GNSSSignals
        # models for B1I uniformly, so the replica would not match them.
        acq_b1i = topN(
            det(
                acquire(
                    BeiDouB1I(),
                    (@view l1h[1:nacq]),
                    III7A_FS,
                    6:37;
                    interm_freq = III7A_IF_B1I,
                    num_coherently_integrated_code_periods = 1,
                    num_noncoherent_accumulations = 20,
                ),
            ),
            4,
        )
        @test !isempty(acq_gps)
        @test !isempty(acq_gal)
        @test !isempty(acq_l5)
        # Galileo reached initial services in December 2016, so E5a-I is
        # unconditionally on the air at this epoch and in the middle of a band
        # this capture already tracks — assert it like the other three.
        @test !isempty(acq_e5a)
        # B1I is not asserted: see `III7A_IF_B1I`. BDS-2 was operational and Hanoi
        # sits inside its regional service area, so the signal is certainly
        # transmitted; whether the L1 front-end passes 1561 MHz is a property of
        # hardware we have no specification for. An empty acquisition here is a
        # fact about the recording, not a regression.
        isempty(acq_b1i) &&
            @info "Flexiband III-7a: no BeiDou B1I acquired 18.9 MHz off L1 centre"
        gps_prns = sort!([a.prn for a in acq_gps]);
        gal_prns = sort!([a.prn for a in acq_gal])
        l5_prns = sort!([a.prn for a in acq_l5])
        e5a_prns = sort!([a.prn for a in acq_e5a]);
        b1i_prns = sort!([a.prn for a in acq_b1i])

        # Multi-band TrackState. The band each group tracks is derived from its
        # signal (`GNSSSignals.get_band_id`), so `track!` takes one
        # `BandMeasurement` per *band*, not per group: GPS L1 C/A and Galileo E1B
        # share `:L1`, GPS L5I and Galileo E5a-I share `:L5`, and BeiDou B1I is
        # alone on `:B1I` — three keys over two front-ends.
        ts = TrackState(;
            signals = (
                gps = (GPSL1CA(),),
                gal = (GalileoE1B(),),
                bds = (BeiDouB1I(),),
                gps_l5 = (GPSL5I(),),
                gal_e5a = (GalileoE5aI(),),
            ),
        )
        for a in acq_gps
            ;
            ts = add_satellite!(ts, a; group = :gps);
        end
        for a in acq_gal
            ;
            ts = add_satellite!(ts, a; group = :gal);
        end
        for a in acq_l5
            ;
            ts = add_satellite!(ts, a; group = :gps_l5);
        end
        for a in acq_e5a
            ;
            ts = add_satellite!(ts, a; group = :gal_e5a);
        end
        for a in acq_b1i
            ;
            ts = add_satellite!(ts, a; group = :bds);
        end
        dec_gps = Dict(p => GPSL1CADecoderState(p) for p in gps_prns)
        dec_gal = Dict(p => GalileoE1BDecoderState(p) for p in gal_prns)
        dec_l5 = Dict(p => GPSL5IDecoderState(p) for p in l5_prns)
        dec_e5a = Dict(p => GalileoE5aDecoderState(p) for p in e5a_prns)
        dec_b1i = Dict(p => BeiDouB1IDecoderState(p) for p in b1i_prns)

        # Track the whole ~15.7 s capture on both front-ends together, feeding
        # each band's soft symbols to its decoders. L1 and L5 come from the same
        # frames, so they stay time-aligned across every multi-band `track!`.
        r, usb = iii7a_open(zippath)
        while true
            l1, l5 = iii7a_read_bands!(usb, nblk(0.2));
            isempty(l1) && break
            # Three band keys, two front-ends: `get_band_id(BeiDouB1I) === :B1I`
            # is distinct from `:L1`, so B1I takes its own `BandMeasurement` over
            # the *same* L1 sample vector at its own IF. Galileo E5a-I reports
            # `:L5` and so rides the L5 measurement already there.
            ts = track!(
                (
                    L1 = BandMeasurement(l1, III7A_FS, III7A_IF),
                    B1I = BandMeasurement(l1, III7A_FS, III7A_IF_B1I),
                    L5 = BandMeasurement(l5, III7A_FS_L5, III7A_IF_L5),
                ),
                ts,
            )
            for p in gps_prns
                s = get_soft_bits(ts, :gps, p);
                isempty(s) || (dec_gps[p] = decode(dec_gps[p], s, length(s)))
            end
            for p in gal_prns
                s = get_soft_bits(ts, :gal, p);
                isempty(s) || (dec_gal[p] = decode(dec_gal[p], s, length(s)))
            end
            for p in l5_prns
                s = get_soft_bits(ts, :gps_l5, p);
                isempty(s) || (dec_l5[p] = decode(dec_l5[p], s, length(s)))
            end
            for p in e5a_prns
                s = get_soft_bits(ts, :gal_e5a, p);
                isempty(s) || (dec_e5a[p] = decode(dec_e5a[p], s, length(s)))
            end
            for p in b1i_prns
                s = get_soft_bits(ts, :bds, p);
                isempty(s) || (dec_b1i[p] = decode(dec_b1i[p], s, length(s)))
            end
        end
        close(r)

        # GPS L1 C/A: at least one satellite yields a parity-validated HOW time of week.
        gps_tows = [
            dec_gps[p].raw_data.TOW for p in gps_prns if !isnothing(dec_gps[p].raw_data.TOW)
        ]
        @test !isempty(gps_tows)
        @test all(t -> 0 <= t <= 604_800, gps_tows)

        # Galileo E1B: at least one satellite yields a CRC-24Q-validated word with TOW + WN.
        gal_ok = [
            (dec_gal[p].raw_data.TOW, dec_gal[p].raw_data.WN) for p in gal_prns if
            !isnothing(dec_gal[p].raw_data.TOW) && !isnothing(dec_gal[p].raw_data.WN)
        ]
        @test !isempty(gal_ok)
        @test all(((t, wn),) -> 0 <= t <= 604_800 && wn >= 0, gal_ok)

        # GPS L5I (CNAV): at least one satellite yields a CRC-24Q-validated
        # message time of week (every message type carries TOW in its header).
        l5_tows =
            [dec_l5[p].raw_data.TOW for p in l5_prns if !isnothing(dec_l5[p].raw_data.TOW)]
        @test !isempty(l5_tows)
        @test all(t -> 0 <= t <= 604_800, l5_tows)

        # Galileo E5a-I (F/NAV): reported, not asserted, and the reason is the
        # page length rather than the signal. An F/NAV page is 500 symbols at
        # 50 sps = 10 s, and page sync needs the *next* page's 12-symbol pattern
        # on top, so a decode costs 10.24 s of continuously tracked symbols. The
        # capture is ~15.7 s and tracking has to converge first, which leaves room
        # for exactly one page window — and only if the 10 s page boundary happens
        # to fall early. Galileo page epochs are common to the constellation, so
        # the extra satellites do not make it independent tries: it is one roll of
        # the dice for all of them. (E1B decodes reliably here because its page
        # parts are 1 s. Word types 1-4 of the five in an F/NAV subframe carry
        # TOW, so a page that *is* decoded very likely yields one.)
        e5a_tows = [
            dec_e5a[p].raw_data.TOW for p in e5a_prns if !isnothing(dec_e5a[p].raw_data.TOW)
        ]
        @test all(t -> 0 <= t <= 604_800, e5a_tows)

        # BeiDou B1I (D1 NAV): reported, not asserted, because acquisition itself
        # is not asserted. If it did acquire, the odds of a decode are good — a D1
        # subframe is 300 bits at 50 bps = 6 s, so two fit inside the capture, and
        # every subframe carries a BCH-checked SOW.
        b1i_sows = [
            dec_b1i[p].raw_data.SOW for p in b1i_prns if !isnothing(dec_b1i[p].raw_data.SOW)
        ]
        @test all(t -> 0 <= t <= 604_800, b1i_sows)

        @info "Decoded real Flexiband III-7a" gps_tow = first(gps_tows) gal = first(gal_ok) l5_tow =
            first(l5_tows) e5a_prns_tracked = e5a_prns e5a_tows = e5a_tows b1i_prns_acquired =
            b1i_prns b1i_sows = b1i_sows
    end
end
