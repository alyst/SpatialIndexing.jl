using SpatialIndexing
using Test, Random

const SI = SpatialIndexing

@testset "SpatialIndexing" begin
    include("regions.jl")
    include("rtree.jl")
    include("examples.jl")
end