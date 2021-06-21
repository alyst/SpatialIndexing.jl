using Documenter, SpatialIndexing

makedocs(
    source = "source",
    format = Documenter.HTML(prettyurls=false),
    sitename = "SpatialIndexing.jl",
    modules = [SpatialIndexing],
    pages = [
        "Introduction" => "index.md",
        "API" => [
            "regions.md",
            "abstract.md",
            "rtree.md",
            "simple.md",
            "query.md"
        ],
    ],
)

deploydocs(
    repo = "github.com/alyst/SpatialIndexing.jl.git",
)
