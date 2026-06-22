# Opt-in real-data integration test.
#
# Decodes the Fraunhofer IIS Flexiband "III-7a" reference capture — the same
# real-sky recording (Hanoi, 2017-10-17, CC BY-NC) that Tracking.jl uses for its
# own integration test — end to end: Acquisition.jl + Tracking.jl produce soft
# symbols, which `GNSSDecoder.decode` turns into navigation data for GPS L1 C/A
# *and* Galileo E1B. This is the real-signal counterpart to the synthetic /
# golden-fixture unit tests, and it exercises both decoders on genuine
# (operational-constellation) signals.
#
# Capture format (multiplexed `.usb`, see Tracking.jl's flexiband test): three
# bands interleaved in 6-byte base-clock chunks; we use only the L1 band — GPS
# L1 / Galileo E1 at 81 MHz, IF = 1575.42 - 1580.0 = -4.58 MHz, 4-bit signed
# I/Q packed I=high nibble, Q=low nibble.
#
# WHY WE ASSERT ON `raw_data`, NOT `data`:
#   `data` is only published once a *complete, cross-validated* positioning set
#   is assembled — for GPS L1 C/A that needs subframes 1+2+3 (~18 s) with
#   consistent IODC, and Galileo E1B needs its full word-type set (~30 s). The
#   `III-7a_short` capture is only ~15.7 s long (its uncompressed `.usb` is
#   ~1.93 GB; despite the website calling it "one minute", the actual sample
#   stream is ~16 s), which is below that span. `raw_data`, however, is
#   populated from individually parity-/CRC-validated subframes and words
#   (GPS HOW TOW gated by the 6-equation parity check; E1B TOW/WN gated by
#   CRC-24Q over the page pair). So a decoded `raw_data.TOW` (+ `raw_data.WN`
#   for E1B) is a genuine, error-checked end-to-end decode — exactly what this
#   test verifies — even though the capture is too short to publish `data`.
#   (Use the L125 210 s captures if a full-ephemeris `data` assertion is wanted.)
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
    using GNSSSignals: GPSL1CA, GalileoE1B
    using Acquisition: acquire, is_detected
    using Tracking: TrackState, BandMeasurement, add_satellite!, track!, get_soft_bits

    const III7A_URL =
        "https://www.iis.fraunhofer.de/content/dam/iis/de/doc/lv/los/lokalisierung/" *
        "SatNAV/Flexiband%20reference%20Data%2020171017_09-51-43_III-7a_short.zip"
    const III7A_SIZE = 1_589_435_475       # exact byte count, for integrity
    const III7A_FS = 81.0e6Hz              # L1 band sampling rate
    const III7A_IF = -4.58e6Hz             # GPS L1 / Galileo E1 within the 1580 MHz band

    # `.usb` framing (Flexiband / ION SDR-metadata format).
    const III7A_BLOCK = 1024; const III7A_HEADER = 6; const III7A_CHUNK = 6
    const III7A_CYCLES = 169; const III7A_PAYLOAD = III7A_CYCLES * III7A_CHUNK
    const III7A_FREQBASE = 20.25e6
    iii7a_nframes(sec) = ceil(Int, sec * III7A_FREQBASE / III7A_CYCLES)

    # Download the capture in verified byte-range chunks (the Fraunhofer server
    # drops long-running connections), assembling only once every chunk is exact.
    function iii7a_fetch()
        dir = get(ENV, "GNSSDECODER_TESTDATA_DIR", tempdir()); mkpath(dir)
        path = joinpath(dir, "Flexiband_III-7a_short.zip")
        (isfile(path) && filesize(path) == III7A_SIZE) && return path
        chunk = 200_000_000; parts = String[]; start = 0; idx = 0
        tmpdir = mktempdir(dir)
        try
            while start < III7A_SIZE
                stop = min(start + chunk - 1, III7A_SIZE - 1); want = stop - start + 1
                part = joinpath(tmpdir, "p$(idx)"); tries = 0
                while !(isfile(part) && filesize(part) == want)
                    (tries += 1) > 30 && error("III-7a download: chunk $idx failed")
                    try
                        Downloads.download(III7A_URL, part; headers = ["Range" => "bytes=$start-$stop"])
                    catch err
                        @warn "III-7a chunk $idx retry $tries" err
                    end
                end
                push!(parts, part); start = stop + 1; idx += 1
            end
            open(path * ".part", "w") do out
                for p in parts; write(out, read(p)); end
            end
            filesize(path * ".part") == III7A_SIZE || error("III-7a assembled size mismatch")
            mv(path * ".part", path; force = true)
        finally
            rm(tmpdir; recursive = true, force = true)
        end
        path
    end

    iii7a_open(zippath) = (r = ZipFile.Reader(zippath);
        (r, only(filter(f -> endswith(f.name, ".usb"), r.files))))
    function iii7a_read_l1!(usb, nblocks)            # -> ComplexF32 L1 samples
        frame = Vector{UInt8}(undef, III7A_BLOCK); chunks = III7A_CYCLES * nblocks
        out = Vector{ComplexF32}(undef, chunks * 4); got = 0; nib(x) = (v = Int(x & 0x0f); v ≥ 8 ? v - 16 : v)
        for _ in 1:nblocks
            eof(usb) && break
            read!(usb, frame); got += 1
            @inbounds for c in 0:(III7A_CYCLES - 1), k in 0:3
                b = frame[III7A_HEADER + c * III7A_CHUNK + 2 + k]
                out[(got - 1) * III7A_CYCLES * 4 + c * 4 + k + 1] = ComplexF32(nib(b >> 4), nib(b))
            end
        end
        resize!(out, got * III7A_CYCLES * 4)
    end

    @testset "Flexiband III-7a real-data decode (GPS L1 C/A + Galileo E1B)" begin
        fshz = ustrip(Hz, III7A_FS); nblk(sec) = iii7a_nframes(sec)
        zippath = iii7a_fetch()
        @test filesize(zippath) == III7A_SIZE

        # Acquire on the first 40 ms.
        r, usb = iii7a_open(zippath); head = iii7a_read_l1!(usb, nblk(0.045)); close(r)
        nacq = round(Int, fshz * 0.040)
        det(a) = filter(x -> is_detected(x; pfa = 1e-8), a)
        topN(a, n) = sort(a; by = x -> -x.CN0)[1:min(n, length(a))]
        acq_gps = topN(det(acquire(GPSL1CA(), (@view head[1:nacq]), III7A_FS, 1:32;
            interm_freq = III7A_IF, num_coherently_integrated_code_periods = 4,
            num_noncoherent_accumulations = 2)), 4)
        acq_gal = det(acquire(GalileoE1B(), (@view head[1:nacq]), III7A_FS, 1:36;
            interm_freq = III7A_IF, num_coherently_integrated_code_periods = 1,
            num_noncoherent_accumulations = 4))
        @test !isempty(acq_gps)
        @test !isempty(acq_gal)
        gps_prns = sort!([a.prn for a in acq_gps]); gal_prns = sort!([a.prn for a in acq_gal])

        ts = TrackState(; signals = (gps = (GPSL1CA(),), gal = (GalileoE1B(),)))
        for a in acq_gps; ts = add_satellite!(ts, a; group = :gps); end
        for a in acq_gal; ts = add_satellite!(ts, a; group = :gal); end
        dec_gps = Dict(p => GPSL1CADecoderState(p) for p in gps_prns)
        dec_gal = Dict(p => GalileoE1BDecoderState(p) for p in gal_prns)

        # Track the whole ~15.7 s capture, feeding soft symbols to the decoders.
        r, usb = iii7a_open(zippath)
        while true
            l1 = iii7a_read_l1!(usb, nblk(0.2)); isempty(l1) && break
            ts = track!(BandMeasurement(l1, III7A_FS, III7A_IF), ts)
            for p in gps_prns
                s = get_soft_bits(ts, :gps, p); isempty(s) || (dec_gps[p] = decode(dec_gps[p], s, length(s)))
            end
            for p in gal_prns
                s = get_soft_bits(ts, :gal, p); isempty(s) || (dec_gal[p] = decode(dec_gal[p], s, length(s)))
            end
        end
        close(r)

        # GPS L1 C/A: at least one satellite yields a parity-validated HOW time of week.
        gps_tows = [dec_gps[p].raw_data.TOW for p in gps_prns if !isnothing(dec_gps[p].raw_data.TOW)]
        @test !isempty(gps_tows)
        @test all(t -> 0 <= t <= 604_800, gps_tows)

        # Galileo E1B: at least one satellite yields a CRC-24Q-validated word with TOW + WN.
        gal_ok = [(dec_gal[p].raw_data.TOW, dec_gal[p].raw_data.WN) for p in gal_prns
                  if !isnothing(dec_gal[p].raw_data.TOW) && !isnothing(dec_gal[p].raw_data.WN)]
        @test !isempty(gal_ok)
        @test all(((t, wn),) -> 0 <= t <= 604_800 && wn >= 0, gal_ok)

        @info "Decoded real Flexiband III-7a" gps_tow = first(gps_tows) gal = first(gal_ok)
    end
end
