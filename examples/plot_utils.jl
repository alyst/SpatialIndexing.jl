# scripts for displaying R-trees using PlotlyJS
using PlotlyJS
using Printf: @sprintf

# webcolor palette for R-tree levels
const LevelsPalette = Dict(
    0 => "#228B22", # Forest Green
    1 => "#DAA520", # Goldenrod
    2 => "#FF6347", # Tomato
    3 => "#B22222", # Firebrick
    4 => "#800080", # Purple
    5 => "#4169E1", # Royal blue
    6 => "#008080",) # Teal

# create plotly trace for a single R-tree `node`
function plotly_trace(node::SI.Node,
                      ix::Union{Integer, Nothing};
                      showlegend::Bool)
    nbr = SI.mbr(node)
    res = scatter(x = [nbr.low[1], nbr.high[1], nbr.high[1], nbr.low[1], nbr.low[1]],
                  y = [nbr.low[2], nbr.low[2], nbr.high[2], nbr.high[2], nbr.low[2]],
                  mode=:lines, line_color=get(LevelsPalette, SI.level(node), "#708090"),
                  line_width=SI.level(node) + 1,
                  hoverinfo=:text, hoveron=:points,
                  text="lev=$(SI.level(node)) ix=$(ix !== nothing ? ix : "none")<br>nchildren=$(length(node)) nelems=$(SI.nelements(node)) area=$(@sprintf("%.2f", SI.area(SI.mbr(node))))",
                  name="level $(SI.level(node))",
                  legendgroup="level $(SI.level(node))", showlegend=showlegend)
    return res
end

# create plotly traces for the nodes in a subtree with the `node` root
# and append them to `node_traces`
function append_subtree_traces!(node_traces::Vector{PlotlyBase.AbstractTrace},
                                node::SI.Node,
                                ix::Union{Integer, Nothing},
                                elems_data::NamedTuple,
                                levels::Set{Int})
    push!(node_traces, plotly_trace(node, ix, showlegend = SI.level(node) ∉ levels))
    push!(levels, SI.level(node)) # show each level once in legend
    if node isa SI.Leaf
        for (i, child) in enumerate(SI.children(node)) # elements are just points
            push!(elems_data.id, SI.id(child))
            push!(elems_data.x, SI.mbr(child).low[1])
            push!(elems_data.y, SI.mbr(child).low[2])
        end
    else
        for (i, child) in enumerate(SI.children(node))
            append_subtree_traces!(node_traces, child, i, elems_data, levels)
        end
    end
    return nothing
end

# create Plotly plot of the given tree
function PlotlyJS.plot(tree::RTree)
    ndims(tree) == 1 && throw(ArgumentError("1-D R-trees not supported"))
    ndims(tree) > 2 && @warn("Only the 1st and 2nd dimensions would be shown for $(ndims(tree))-D trees")
    node_traces = Vector{PlotlyBase.AbstractTrace}()
    elems_data = (id = Vector{SI.idtype(tree)}(), x = Vector{Float64}(), y = Vector{Float64}())
    append_subtree_traces!(node_traces, tree.root, nothing, elems_data, Set{Int}())
    elems_trace = scatter(mode=:markers, name=:data, marker_color = "#333333", marker_size = 2,
                          x=elems_data.x, y=elems_data.y, text=["id=$id" for id in elems_data.id],
                          hoverinfo=:text)
    PlotlyJS.plot([node_traces; [elems_trace]], Layout(hovermode=:closest))
end
