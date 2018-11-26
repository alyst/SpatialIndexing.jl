# scripts for displaying R-trees using PlotlyJS
using PlotlyJS
using Printf: @sprintf

# webcolor palette for R-tree levels
const LevelsPalette = Dict(
    1 => "#228B22", # Forest Green
    2 => "#DAA520", # Goldenrod
    3 => "#FF6347", # Tomato
    4 => "#B22222", # Firebrick
    5 => "#800080", # Purple
    6 => "#4169E1", # Royal blue
    7 => "#008080",) # Teal

# create plotly trace for a single R-tree `node` (rectangle edges)
function node_trace(node::SI.Node{<:Any,2}, ix::Union{Integer, Nothing};
                    showlegend::Bool)
    nbr = SI.mbr(node)
    res = scatter(x = [nbr.low[1], nbr.high[1], nbr.high[1], nbr.low[1], nbr.low[1]],
                  y = [nbr.low[2], nbr.low[2], nbr.high[2], nbr.high[2], nbr.low[2]],
                  mode=:lines, line_color=get(LevelsPalette, SI.level(node), "#708090"),
                  line_width=SI.level(node),
                  hoverinfo=:text, hoveron=:points,
                  text="lev=$(SI.level(node)) ix=$(ix !== nothing ? ix : "none")<br>nchildren=$(length(node)) nelems=$(SI.nelements(node)) area=$(@sprintf("%.2f", SI.area(SI.mbr(node))))",
                  name="level $(SI.level(node))",
                  legendgroup="level $(SI.level(node))", showlegend=showlegend)
    return res
end

# create plotly trace for a single R-tree `node` (cube edges)
function node_trace(node::SI.Node, ix::Union{Integer, Nothing};
                    showlegend::Bool)
    nbr = SI.mbr(node)
    res = scatter3d(x = [nbr.low[1],  nbr.high[1], nbr.high[1], nbr.low[1], nbr.low[1],
                         nbr.low[1],  nbr.high[1], nbr.high[1], nbr.low[1], nbr.low[1],
                         NaN, nbr.low[1], nbr.low[1],
                         NaN, nbr.high[1], nbr.high[1],
                         NaN, nbr.high[1], nbr.high[1]],
                    y = [nbr.low[2],  nbr.low[2],  nbr.high[2], nbr.high[2], nbr.low[2],
                         nbr.low[2],  nbr.low[2],  nbr.high[2], nbr.high[2], nbr.low[2],
                         NaN, nbr.high[2], nbr.high[2],
                         NaN, nbr.high[2], nbr.high[2],
                         NaN, nbr.low[2], nbr.low[2]],
                    z = [nbr.low[3],  nbr.low[3],  nbr.low[3],  nbr.low[3], nbr.low[3],
                         nbr.high[3], nbr.high[3], nbr.high[3], nbr.high[3], nbr.high[3],
                         NaN, nbr.low[3], nbr.high[3],
                         NaN, nbr.low[3], nbr.high[3],
                         NaN, nbr.low[3], nbr.high[3]],
                  mode=:lines, line_color=get(LevelsPalette, SI.level(node), "#708090"),
                  line_width=SI.level(node),
                  hoverinfo=:text, hoveron=:points,
                  text="lev=$(SI.level(node)) ix=$(ix !== nothing ? ix : "none")<br>nchildren=$(length(node)) nelems=$(SI.nelements(node)) area=$(@sprintf("%.2f", SI.area(SI.mbr(node))))",
                  name="level $(SI.level(node))",
                  legendgroup="level $(SI.level(node))", showlegend=showlegend)
    return res
end

# create plotly traces for the nodes in a subtree with the `node` root
# and append them to `node_traces`
function append_subtree_traces!(node_traces::Vector{PlotlyBase.AbstractTrace},
                                node::SI.Node, ix::Union{Integer, Nothing},
                                levels::Set{Int})
    push!(node_traces, node_trace(node, ix, showlegend = SI.level(node) ∉ levels))
    push!(levels, SI.level(node)) # show each level once in legend
    if node isa SI.Branch
        for (i, child) in enumerate(SI.children(node))
            append_subtree_traces!(node_traces, child, i, levels)
        end
    end
    return nothing
end

data_trace(tree::RTree{<:Any, 2}) =
    scatter(mode=:markers, name=:data, marker_color = "#333333", marker_size = 2,
            x=[SI.center(SI.mbr(x)).coord[1] for x in tree],
            y=[SI.center(SI.mbr(x)).coord[2] for x in tree],
            text=["id=$(SI.id(x))" for x in tree], hoverinfo=:text)

data_trace(tree::RTree) =
    scatter3d(mode=:markers, name=:data, marker_color = "#333333", marker_size = 2,
              x=[SI.center(SI.mbr(x)).coord[1] for x in tree],
              y=[SI.center(SI.mbr(x)).coord[2] for x in tree],
              z=[SI.center(SI.mbr(x)).coord[3] for x in tree],
              text=["id=$(SI.id(x))" for x in tree], hoverinfo=:text)

# create Plotly plot of the given tree
function PlotlyJS.plot(tree::RTree)
    ndims(tree) == 1 && throw(ArgumentError("1-D R-trees not supported"))
    ndims(tree) > 3 && @warn("Only 1st-3rd dimensions would be shown for $(ndims(tree))-D trees")
    node_traces = Vector{PlotlyBase.AbstractTrace}()
    append_subtree_traces!(node_traces, tree.root, nothing, Set{Int}())
    PlotlyJS.plot([node_traces; [data_trace(tree)]], Layout(hovermode=:closest))
end
