"""
Base class for `RTree` node (`Branch` or `Leaf`).
"""
abstract type Node{T,N,V} end

dimtype(::Type{<:Node{T}}) where T = T
dimtype(node::Node) = dimtype(typeof(node))
Base.ndims(::Type{<:Node{<:Any, N}}) where N = N
Base.ndims(node::Node) = ndims(typeof(node))
mbrtype(::Type{<:Node{T,N}}) where {T,N} = Rect{T,N}
mbrtype(node::Node) = mbrtype(typeof(node))
mbr(node::Node) = node.mbr
parent(node::Node) = node.parent
hasparent(node::Node) = node.parent !== nothing
children(node::Node) = node.children
Base.isempty(node::Node) = isempty(children(node))
Base.length(node::Node) = length(children(node))
Base.getindex(node::Node, i::Integer) = getindex(children(node), i)

# `node` position among its siblings in `children` vector of the parent or nothing
pos_in_parent(node::Node) =
    hasparent(node) ? findfirst(c -> c === node, children(parent(node))) : nothing

# update mbr to match the combined MBR of its children
syncmbr!(node::Node) =
    node.mbr = !isempty(node) ? mapreduce(mbr, combine, node.children) : empty(mbrtype(node))

# check whether type E could be used as the data element type in R-tree
# FIXME also need to check whether E could be converted to V of Node
check_eltype_rtree(::Type{E}, ::Type{N}) where {N<:Node, E} =
    check_hasmbr(mbrtype(N), E)

check_eltype_rtree_eltype(el::Any, node::Node) =
    check_eltype_rtree(typeof(el), typeof(node))

# equality check for R-tree elements (by MBR and, optionally, ID)
isequal_rtree(el::Any, reg::Region, id::Any = nothing) =
    isequal_rtree(idtrait(typeof(el)), el, reg, id)

isequal_rtree(::Type{<:HasID}, el::Any, reg::Region, key::Any) =
    isequal(id(el), convert(idtype(idtrait(typeof(el))), key)) && (mbr(el) == reg)

isequal_rtree(::Type{HasNoID}, el::Any, reg::Region, ::Nothing = nothing) =
    (mbr(el) == reg)

"""
R-Tree leaf (level 1 node).
Its children are data elements of type `V`.
"""
mutable struct Leaf{T,N,V} <: Node{T,N,V}
    parent::Union{Node{T,N,V}, Nothing} # actually Branch{T,N,V}, but cannot do without forward decl
    mbr::Rect{T,N}
    children::Vector{V}

    function Leaf{T,N,V}(parent::Union{Nothing, Node{T,N,V}},
                         ::Type{V}) where {T,N,V}
        new{T,N,V}(parent, empty(Rect{T,N}), Vector{V}())
    end
end

level(node::Leaf) = 1 # always

function Base.setindex!(leaf::Leaf, el::Any, i::Integer)
    check_eltype_rtree(el)
    setindex!(children(node), child, i)
end

"""
R-Tree node for levels above 1 (non-`Leaf`).
"""
mutable struct Branch{T,N,V} <: Node{T,N,V}
    parent::Union{Branch{T,N,V}, Nothing}
    level::Int  # level (thickness) in the R-tree (0=leaves)
    mbr::Rect{T,N}
    children::Union{Vector{Branch{T,N,V}},Vector{Leaf{T,N,V}}}

    function Branch{T,N,V}(parent::Union{Nothing, Branch{T,N,V}}, ::Type{C}) where {T,N,V, C<:Node{T,N,V}}
        new{T,N,V}(parent, 0, empty(Rect{T,N}), Vector{C}())
    end
end

level(node::Branch) = node.level

function Base.setindex!(node::Branch, child::Node, i::Integer)
    (level(node) == level(child) + 1) ||
        throw(SpatialIndexException("Parent ($(level(node)) and child ($(level(child))) levels don't match"))
    setindex!(children(node), child, i)
    child.parent = node
end
