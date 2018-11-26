using Revise
using Random, DataFrames, PlotlyJS, SpatialIndexing
using Printf: @sprintf

const SI = SpatialIndexing

Random.seed!(1)
bulk_tree = RTree{Float64, 2}(Int, String, leaf_capacity = 20, branch_capacity = 20)
mbrs = Vector{SI.mbrtype(bulk_tree)}()
for i in 1:100000
    t = 10*randn()
    u, v = 0.01 .* randn(2)
    x = sin(t)^3 + u
    y = cos(t)^3 + v
    rmbr = SI.Rect((x, y), (x, y))
    push!(mbrs, rmbr)
end
seq_tree = RTree{Float64, 2}(Int, String, leaf_capacity = 20, branch_capacity = 20)
for (i, br) in enumerate(mbrs)
    insert!(seq_tree, br, i, string(i))
end
SI.check(seq_tree)

rect = SI.Rect((-1.0, -2.0), (2.0, -0.5))
n_in_rect = sum(br -> in(br, rect), mbrs)
SI.subtract!(seq_tree, rect);
SI.check(seq_tree)
@show seq_tree.nelems seq_tree.nelem_deletions seq_tree.nnode_reinsertions seq_tree.nnode_splits

SI.load!(bulk_tree, enumerate(mbrs), getid=x->x[1], getmbr=x->x[2], getval=x->string(x[1]));
SI.check(bulk_tree)

# insert a point `pt` into pareto frontier
# (i.e. only if `tree` does not contain points that dominate it)
# remove all points that are dominated by `pt`
# return true if `pt` was inserted
function pareto_insert!(tree::RTree{T,N}, pt::SI.Point{T,N}, key, val) where {T,N}
    betterbr = SI.Rect(ntuple(_ -> typemin(T), N), pt.coord)
    # check if there are points that dominate pt
    isempty(tree, betterbr) || return false
    # remove the dominated points
    worsebr = SI.Rect(pt.coord, ntuple(_ -> typemax(T), N))
    SI.subtract!(tree, worsebr)
    # insert pt
    insert!(tree, pt, key, val)
end

seq_tree = RTree{Float64, 2}(Int, String, leaf_capacity = 8, branch_capacity = 8)
for (i, br) in enumerate(mbrs)
    #insert!(seq_tree, i, br, string(i))
    pareto_insert!(seq_tree, SI.Point(br.low), i, string(i))
end
SI.check(seq_tree)

include(joinpath(@__DIR__, "plot_utils.jl"))

seq_tree_plot = plot(seq_tree);
open(joinpath(@__DIR__, "pareto_rtree_seq.html"), "w") do io
    PlotlyJS.savehtml(io, seq_tree_plot, :embed)
end
