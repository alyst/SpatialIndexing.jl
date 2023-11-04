include("../examples/pareto.jl")
include("../examples/spiral.jl")

using ReferenceTests, FileIO 
@test_reference "spiral_rtree_seq.png" load("../examples/spiral_rtree_seq.png")
@test_reference "pareto3d_rtree_bulk.jpg" load("../examples/pareto3d_rtree_bulk.jpg")