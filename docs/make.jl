using YAXArraysToolbox
using Documenter

DocMeta.setdocmeta!(YAXArraysToolbox, :DocTestSetup, :(using YAXArraysToolbox); recursive=true)

makedocs(;
    modules=[YAXArraysToolbox],
    authors="Pabon-Moreno, Daniel E.; Duveiller, Gregory; Gans, Fabian; Winkler, Alexander",
    repo="https://github.com/dpabon/YAXArraysToolbox.jl/blob/{commit}{path}#{line}",
    sitename="YAXArraysToolbox.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://dpabon.github.io/YAXArraysToolbox.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/dpabon/YAXArraysToolbox.jl",
    devbranch="main",
)
