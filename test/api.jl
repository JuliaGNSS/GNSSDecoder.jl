# The API is what the reference documents: its bare names are the exports, its
# qualified ones the unexported names a caller binds by hand.
@testset "the API reference lists the exports bare and the rest qualified" begin
    bare, qualified = Symbol[], Symbol[]
    in_docs = false
    for line in eachline(joinpath(@__DIR__, "..", "docs", "src", "api.md"))
        line == "```@docs" && (in_docs = true; continue)
        line == "```" && (in_docs = false; continue)
        (in_docs && !startswith(line, "Base.")) || continue
        if startswith(line, "GNSSDecoder.")
            push!(qualified, Symbol(chopprefix(line, "GNSSDecoder.")))
        else
            push!(bare, Symbol(line))
        end
    end
    # Without `all`, `names` lists the exports.
    exported = setdiff(names(GNSSDecoder), [:GNSSDecoder])
    @test sort(bare) == sort(exported)
    @test all(name -> isdefined(GNSSDecoder, name) && name ∉ exported, qualified)
end
