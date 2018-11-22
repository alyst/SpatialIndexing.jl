# scripts for converting R-tree into `DataFrame`
using DataFrames

# appends node/element information into the `dfcols`
function append_node_info!(dfcols::NamedTuple, node::Any,
                           parent::Union{SI.Node, Missing},
                           ix::Union{Integer, Missing})
    push!(dfcols.type, string(typeof(node).name.name))
    push!(dfcols.level, isa(node, SI.Node) ? SI.level(node) : missing)
    push!(dfcols.ix, ix)
    push!(dfcols.id, isa(node, SI.Node) ? string(Int(pointer_from_objref(node)), base=16) :
         (SI.idtrait(typeof(node)) !== SI.HasNoID ? repr(SI.id(node)) : missing))
    push!(dfcols.parent, !ismissing(parent) ? string(Int(pointer_from_objref(parent)), base=16) : missing)
    push!(dfcols.area, SI.area(SI.mbr(node)))
    push!(dfcols.xmin, SI.mbr(node).low[1])
    push!(dfcols.xmax, SI.mbr(node).high[1])
    push!(dfcols.ymin, SI.mbr(node).low[2])
    push!(dfcols.ymax, SI.mbr(node).high[2])
    return dfcols
end

# appends information about all the nodes in the `node` subtree
# into the `dfcols`
function append_subtree_info!(dfcols::NamedTuple, node::Any,
                              parent::Union{SI.Node, Missing},
                              ix::Union{Integer, Missing})
    append_node_info!(dfcols, node, parent, ix)
    if node isa SI.Node
        for (i, child) in enumerate(SI.children(node))
            append_subtree_info!(dfcols, child, node, i)
        end
    end
    return node
end

# convert `RTree` into `DataFrame`
function Base.convert(::Type{DataFrame}, tree::RTree)
    dfcols = (type = Vector{String}(), level = Vector{Union{Int, Missing}}(),
              id = Vector{Union{String, Missing}}(),
              ix = Vector{Union{Int, Missing}}(),
              parent = Vector{Union{String, Missing}}(),
              area =Vector{Float64}(),
              xmin = Vector{Float64}(), ymin = Vector{Float64}(),
              xmax = Vector{Float64}(), ymax = Vector{Float64}())
    append_subtree_info!(dfcols, tree.root, missing, missing)
    return DataFrame(dfcols)
end
