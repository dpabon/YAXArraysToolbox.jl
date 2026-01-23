using Documenter
using YAXArraysToolbox

# Build documentation
makedocs(
    sitename = "YAXArraysToolbox.jl",
    authors = "Daniel E. Pabon-Moreno and contributors",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://dpabon.github.io/YAXArraysToolbox.jl",
        assets = [
            "assets/custom.css",
        ],
        sidebar_sitename = true,
        collapselevel = 2,
        highlights = ["julia"],
        ansicolor = true,
    ),
    modules = [YAXArraysToolbox],
    pages = [
        "Home" => "index.md",
        "Tutorials" => [
            "Basic Operations" => "tutorials/basic_operations.md",
            "Space-for-Time Method" => "tutorials/space4time_proof_of_concept.md",
        ],
        "API Reference" => "api.md",
    ],
    warnonly = [:missing_docs, :cross_references],
    checkdocs = :exports,
)

# Deploy documentation
deploydocs(
    repo = "github.com/dpabon/YAXArraysToolbox.jl.git",
    devbranch = "main",
    push_preview = true,
)