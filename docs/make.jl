using FastaProcessing
using Documenter

DocMeta.setdocmeta!(FastaProcessing, :DocTestSetup, :(using FastaProcessing); recursive=true)

makedocs(;
    modules=[FastaProcessing],
    authors="Shane Kuei-Hsien Chu (skchu@wustl.edu)",
    sitename="FastaProcessing.jl",
    format=Documenter.HTML(;
        canonical="https://kchu25.github.io/FastaProcessing.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/kchu25/FastaProcessing.jl",
    devbranch="main",
)
