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
