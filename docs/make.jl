using Documenter
using Tapir

makedocs(
    sitename = "Tapir.jl",
    format = Documenter.HTML(prettyurls = false),
    pages = [
        "Introduction" => "index.md",
        "API" => "api.md"
    ]
)

deploydocs(
    repo = "https://github.com/LucianNita/Tapir.jl.git",
    devbranch = "main"
)