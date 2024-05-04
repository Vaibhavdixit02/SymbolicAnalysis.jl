using Documenter, DocumenterVitepress

using SymbolicAnalysis

makedocs(;
    modules=[SymbolicAnalysis],
    authors="Vaibhav Dixit, Shashi Gowda",
    repo="https://github.com/Vaibhavdixit02/SymbolicAnalysis.jl",
    sitename="SymbolicAnalysis.jl",
    format=DocumenterVitepress.MarkdownVitepress(
        repo = "https://github.com/Vaibhavdixit02/SymbolicAnalysis.jl",
        devurl = "dev",
        deploy_url = "vaibhavdixit02.github.io/SymbolicAnalysis.jl",
    ),
    pages=[
        "Home" => "index.md",
        "Atoms" => "atoms.md",
        "Functions" => "functions.md",
    ],
    warnonly = true,
)

deploydocs(;
    repo="github.com/Vaibhavdixit02/SymbolicAnalysis.jl",
    push_preview=true,
)