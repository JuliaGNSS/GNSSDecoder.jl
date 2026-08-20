using Documenter
using GNSSDecoder

DocMeta.setdocmeta!(GNSSDecoder, :DocTestSetup, :(using GNSSDecoder); recursive=true)

makedocs(
    sitename = "GNSSDecoder.jl",
    modules = [GNSSDecoder],
    authors = "Soeren Schoenbrod, Michael Niestroj, Erik Deinzer",
    format = Documenter.HTML(
        canonical = "https://JuliaGNSS.github.io/GNSSDecoder.jl",
        edit_link = "master",
        # api.md documents every per-signal data container on one page and has
        # outgrown Documenter's 200 KiB default with the BeiDou decoders.
        size_threshold = 400 * 1024,
        size_threshold_warn = 300 * 1024,
    ),
    pages = [
        "Home" => "index.md",
        "API Reference" => "api.md",
    ],
    doctest = true,
    checkdocs = :exports,
)

deploydocs(
    repo = "github.com/JuliaGNSS/GNSSDecoder.jl.git",
    devbranch = "master",
    push_preview = true,
)
