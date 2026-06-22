# Opt-in real-data integration test.
#
# Decodes the Fraunhofer IIS Flexiband "L125 III-1b" reference capture (a Spirent
# GPS/Galileo simulation, CC BY-NC) end to end: it acquires and tracks the signal
# with Acquisition.jl + Tracking.jl, feeds Tracking's soft symbols into
# `GNSSDecoder.decode`, and asserts that validated GPS L1 C/A navigation data
# (TOW, week, ephemeris) is recovered — the real-signal counterpart to the
# synthetic/golden-fixture tests.
#
# The capture is a 3.45 GB zip. To keep CI light we never download all of it:
# we HTTP range-fetch only the first ~680 MB, then raw-DEFLATE-stream the ZIP
# entry's prefix (the local file header at offset 0 is sufficient — no central
# directory needed) and stop after ~35 s of samples.
#
# Demuxed L1/E1 format (per the Fraunhofer page): complex baseband (IF = 0),
# FS = 20 MHz, interleaved signed 8-bit I then Q.
#
# OPT-IN: skipped unless `ENV["GNSSDECODER_RUN_INTEGRATION_TEST"] == "true"`, so
# a plain `]test` never triggers the download. CI enables it on one matrix job
# and caches the fetched prefix (set `GNSSDECODER_TESTDATA_DIR` to the cached
# directory).

if get(ENV, "GNSSDECODER_RUN_INTEGRATION_TEST", "false") != "true"
    @info "Skipping Flexiband L125 integration test (range-fetches ~680 MB). " *
          "Set ENV[\"GNSSDECODER_RUN_INTEGRATION_TEST\"] = \"true\" to run it."
else
    using Downloads: Downloads
    using CodecZlib: DeflateDecompressorStream
    using Unitful: Hz, ustrip
    using GNSSSignals: GPSL1CA
    using Acquisition: acquire, is_detected
    using Tracking: TrackState, BandMeasurement, add_satellite!, track!, get_soft_bits

    const L125_URL =
        "https://www2.iis.fraunhofer.de/flexiband/reference-data/L125_III1b_210s_L1.bin.zip"
    const L125_FS = 20.0e6Hz       # demuxed L1/E1 sampling rate
    const L125_IF = 0.0Hz          # complex baseband
    const L125_TRACK_SECONDS = 35.0
    const L125_PREFIX_BYTES = 680_000_000   # ~40 s of samples after DEFLATE

    _l125_fshz() = ustrip(Hz, L125_FS)

    # Range-fetch just the prefix of the zip (cached; CI points the dir at a cache).
    function l125_fetch_prefix()
        dir = get(ENV, "GNSSDECODER_TESTDATA_DIR", tempdir())
        mkpath(dir)
        path = joinpath(dir, "L125_III1b_210s_L1_prefix.zip")
        (isfile(path) && filesize(path) >= L125_PREFIX_BYTES) && return path
        tmp = path * ".part"
        Downloads.download(L125_URL, tmp; headers = ["Range" => "bytes=0-$(L125_PREFIX_BYTES-1)"])
        mv(tmp, path; force = true)
        path
    end

    # Open the single .bin entry as a streaming raw-DEFLATE decompressor.
    function l125_open_deflate(path)
        io = open(path)
        read(io, UInt32); skip(io, 4)         # signature + version/flags
        method = read(io, UInt16); skip(io, 16)  # method, then time/date/crc/sizes
        nlen = read(io, UInt16); elen = read(io, UInt16)
        skip(io, Int(nlen) + Int(elen))        # filename + extra -> compressed data
        method == 8 || error("expected DEFLATE-compressed zip entry, got method $method")
        DeflateDecompressorStream(io)
    end

    function l125_read_samples!(ds, nsamp)
        raw = Vector{UInt8}(undef, 2 * nsamp)
        n = readbytes!(ds, raw, 2 * nsamp)
        m = n ÷ 2
        out = Vector{ComplexF32}(undef, m)
        @inbounds for k in 1:m
            out[k] = ComplexF32(reinterpret(Int8, raw[2k - 1]), reinterpret(Int8, raw[2k]))
        end
        out
    end

    @testset "Flexiband L125 real-data GPS L1 C/A decode (streamed)" begin
        fshz = _l125_fshz()
        nchunk(sec) = round(Int, fshz * sec)

        path = l125_fetch_prefix()
        @test filesize(path) >= L125_PREFIX_BYTES

        # Acquire on the first 40 ms.
        ds = l125_open_deflate(path)
        head = l125_read_samples!(ds, nchunk(0.045))
        close(ds)
        nacq = round(Int, fshz * 0.040)
        detected = filter(a -> is_detected(a; pfa = 1e-8),
            acquire(GPSL1CA(), (@view head[1:nacq]), L125_FS, 1:32; interm_freq = L125_IF,
                num_coherently_integrated_code_periods = 4, num_noncoherent_accumulations = 2))
        strongest = sort(detected; by = a -> -a.CN0)[1:min(4, length(detected))]
        @test length(strongest) >= 3        # the Spirent sim broadcasts many strong SVs
        prns = sort!([a.prn for a in strongest])

        track_state = TrackState(; signals = (gps = (GPSL1CA(),),))
        for a in strongest
            track_state = add_satellite!(track_state, a; group = :gps)
        end
        decoders = Dict(p => GPSL1CADecoderState(p) for p in prns)

        # Stream + track in 0.2 s chunks, feeding each satellite's soft bits to its decoder.
        ds = l125_open_deflate(path)
        elapsed = 0.0
        while elapsed < L125_TRACK_SECONDS
            samples = l125_read_samples!(ds, nchunk(0.2))
            isempty(samples) && break
            track_state = track!(BandMeasurement(samples, L125_FS, L125_IF), track_state)
            for p in prns
                soft = get_soft_bits(track_state, :gps, p)
                isempty(soft) || (decoders[p] = decode(decoders[p], soft, length(soft)))
            end
            elapsed += length(samples) / fshz
        end
        close(ds)

        # At least one strong satellite must yield validated nav data with a
        # physically sane GPS ephemeris (sqrt_A ≈ 5153.7 √m, the ~26,560 km orbit).
        decoded = [decoders[p].data for p in prns if !isnothing(decoders[p].data.TOW)]
        @test !isempty(decoded)
        for d in decoded
            @test 0 <= d.TOW <= 604_800
            @test d.trans_week !== nothing
            @test 5150 < d.sqrt_A < 5160
        end
        @info "Decoded GPS L1 C/A from Flexiband L125" satellites = length(decoded) TOW =
            first(decoded).TOW week = first(decoded).trans_week sqrt_A = first(decoded).sqrt_A
    end
end
