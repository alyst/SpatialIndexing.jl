using Revise
using Random, DataFrames, PlotlyJS, SpatialIndexing
using Printf: @sprintf

const SI = SpatialIndexing

Random.seed!(1)
pts = Vector{SI.Point{Float64, 3}}()
for i in 1:100000
    r = 1 - 1/(1+abs(randn()))
    a = rand()*pi/2
    b = rand()*pi/2
    x = (sin(a) * cos(b))^3 * r
    y = (cos(a) * cos(b))^3 * r
    z = (sin(b))^3 * r
    push!(pts, SI.Point((x, y, z)))
end
seq_tree = RTree{Float64, 3}(Int, String, leaf_capacity = 20, branch_capacity = 20)
for (i, pt) in enumerate(pts)
    insert!(seq_tree, pt, i, string(i))
end
SI.check(seq_tree)

#=
rect = SI.Rect((-1.0, -2.0, -3.0), (2.0, -0.5, 1.0))
n_in_rect = sum(br -> in(br, rect), mbrs)
SI.subtract!(seq_tree, rect);
SI.check(seq_tree)
@show seq_tree.nelems seq_tree.nelem_deletions seq_tree.nnode_reinsertions seq_tree.nnode_splits
=#

# insert a point `pt` into pareto frontier
# (i.e. only if `tree` does not contain points that dominate it)
# remove all points that are dominated by `pt`
# return true if `pt` was inserted
function pareto_insert!(tree::RTree{T,N}, pt::SI.Point{T,N}, key, val) where {T,N}
    betterbr = SI.Rect(pt.coord, ntuple(_ -> typemax(T), N))
    # check if there are points that dominate pt
    isempty(tree, betterbr) || return false
    # remove the dominated points
    worsebr = SI.Rect(ntuple(_ -> typemin(T), N), pt.coord)
    SI.subtract!(tree, worsebr)
    # insert pt
    insert!(tree, pt, key, val)
end

pareto_tree = RTree{Float64, 3}(Int, String, leaf_capacity = 8, branch_capacity = 8)
for (i, pt) in enumerate(pts)
    pareto_insert!(pareto_tree, pt, i, string(i))
end
SI.check(pareto_tree)

include(joinpath(@__DIR__, "plot_utils.jl"))

pareto_tree_plot = plot(pareto_tree);
open(joinpath(@__DIR__, "pareto3d_rtree_seq.html"), "w") do io
    PlotlyJS.savehtml(io, pareto_tree_plot, :embed)
end

bulk_pareto_tree = RTree{Float64, 3}(Int, String, leaf_capacity = 8, branch_capacity = 8)
SI.load!(bulk_pareto_tree, pareto_tree)
SI.check(bulk_pareto_tree)

bulk_pareto_tree_plot = plot(bulk_pareto_tree);
open(joinpath(@__DIR__, "pareto3d_rtree_bulk.html"), "w") do io
    PlotlyJS.savehtml(io, bulk_pareto_tree_plot, :embed)
end

using ORCA
savefig(bulk_pareto_tree_plot, joinpath(@__DIR__, "pareto3d_rtree_bulk.png"), width=900, height=1000)
