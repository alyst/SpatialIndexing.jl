"""
R-Tree variants.
"""
@enum RTreeVariant::Int RTreeLinear = 0 RTreeQuadratic = 1 RTreeStar = 2

"""
Pool of the `T1` R-tree nodes (`Leaf` or `Branch`) with `T2` children
(`Branch`, `Leaf` or `Elem`).
It allows reusing the deleted nodes and reduce the stress on GC.
"""
mutable struct NodePool{T1<:Node, T2} <: AbstractPool{T1}
    objs::Vector{T1}

    capacity::Int       # pool capacity
    node_capacity::Int  # node capacity (max number of children)

    NodePool{T1,T2}(pool_capacity::Integer, node_capacity::Integer) where {T1<:Node, T2} =
        new{T1,T2}(sizehint!(Vector{T1}(), pool_capacity), pool_capacity, node_capacity)
end

const BranchPool{T,N,V} = NodePool{Branch{T,N,V}, Branch{T,N,V}}
const TwigPool{T,N,V} = NodePool{Branch{T,N,V}, Leaf{T,N,V}}
const LeafPool{T,N,V} = NodePool{Leaf{T,N,V}, V}

node_capacity(pool::NodePool) = pool.node_capacity
newelem(pool::NodePool{T1,T2}) where {T1, T2} = T1(nothing, T2)

"""
R-Tree: `N`-dimensional spatial data index [guttman84].

R-tree groups data elements (`V`) into leaves (`Leaf`) and leaves into branches (`Branch`).
It uses various heuristics to ensure that the minimal bounding rectangles (MBRs)
of the nodes (`Rect{T,N}` rectangles that encompass the data elements attached to these nodes)
stay compact and that the MBRs of the nodes that are on the same level of
R-tree hierarchy have minimal overlap with each other.
This property makes R-trees efficient for spatial queries.

To facilitate spatial indexing, the `V` data elements need to support `HasMBR`
trait (i.e. define `mbrtype(V)` and `mbr(v::V)` methods) and, optionally,
`HasID` trait (via `idtype(V)` and `id(v::V)` methods). `mbr(v::V)` should
return minimal bounding rectangle (MBR) of type `Rect{T,N}` that contains `v`.
`SpatialElem{T,N,D}` type provides generic implementation of spatial data element
that explicitly stores `id`, `mbr` and data object of type `D` and implements
`HasMBR` and `HasID` traits.

# Parameters

The behaviour of `RTree` is defined by the parameters supplied at its creation:
  * `T`: the numeric type for the spatial coordinate
  * `N`: the number of spatial dimensions
  * `variant`: one of `RTreeLinear`, `RTreeQuadratic`, or `RTreeStar` (default)
  * `tight_mbrs`: recalculate node MBR when the child is removed (default is `true`)
  * `branch_capacity`: capacity of branch nodes (default is `100`)
  * `leaf_capacity`: capacity of leaf nodes (default is `100`)
  * `leafpool_capacity`: How many detached zero level nodes (leaves) should be kept for reuse (default is `100`)
  * `twigpool_capacity`: How many detached first level nodes should be kept for reuse (default is `100`)
  * `branchpool_capacity`: How many other (level > 1) detached branch nodes should be kept for reuse (default is `100`)
  * `nearmin_overlap`: How many candidates to consider when identifying the node with minimal overlap (default is `32`)
  * `fill_factor`: How much should the node be filled (fraction of its capacity) after splitting (default is `0.7`)
  * `splitdistribution_factor`: How much can the sizes of the two nodes differ after splitting (default is `0.4`)
  * `reinsert_factor`: How much should the node be underfilled (fraction of its capacity)
    to consider removing it and redistributing its children to other nodes (default is `0.3`)

# Performance

The nodes in R-tree have limited capacity (maximual number of children)
specified at `RTree` creation (`leaf_capacity` and `branch_capacity`).
Larger capacities results in shorter trees, but they time required to locate
the specific spatial region grows linearly with the capacity.

# References
[guttman84] “R-Trees: A Dynamic Index Structure for Spatial Searching”
    A. Guttman, Proc. 1984 ACM-SIGMOD Conference on Management of
    Data (1985), 47-57.
[beckmann90] "The R*-tree: an efficient and robust access method for points and rectangles"
    N. Beckmann, H.P. Kriegel, R. Schneider, B. Seeger, Proc. 1990 ACM SIGMOD
    international conference on Management of data (1990), p.322
"""
mutable struct RTree{T,N,V} <: SpatialIndex{T,N,V}
    # R-tree insert/update policy
    variant::RTreeVariant
    tight_mbrs::Bool        # whether MBR are always updated upon child change/removal

    nearmin_overlap::Int
    fill_factor::Float64
    splitdistribution_factor::Float64
    reinsert_factor::Float64

    root::Node{T,N}     # root node

    nelems::Int         # number of data elements in a tree (`length()`)
    nnodes_perlevel::Vector{Int}    # number of nodes (branches and leaves) at each tree level
    nelem_insertions::Int           # number of elements insertions
    nelem_deletions::Int            # number of elements deletions
    nnode_splits::Int               # number of node splits
    nnode_reinsertions::Int         # number of node reinsertions

    branchpool::BranchPool{T,N,V}
    twigpool::TwigPool{T,N,V}
    leafpool::LeafPool{T,N,V}

    function RTree{T,N,V}(;
        variant::RTreeVariant = RTreeStar,
        tight_mbrs::Bool = true,
        branch_capacity::Integer = 100,
        leaf_capacity::Integer = 100,
        leafpool_capacity::Integer = 100,
        branchpool_capacity::Integer = 100,
        nearmin_overlap::Integer = floor(Int, 0.32 * leaf_capacity),
        fill_factor::Real = 0.7,
        splitdistribution_factor::Real = 0.4,
        reinsert_factor::Real = 0.3
    ) where {T<:Number,N,V}
        check_eltype_rtree(V, Leaf{T,N,V})
        leafpool = LeafPool{T,N,V}(leafpool_capacity, leaf_capacity)
        new{T,N,V}(variant, tight_mbrs,
                   nearmin_overlap, fill_factor,
                   splitdistribution_factor, reinsert_factor,
                   acquire!(leafpool),
                   0, [1], # init with 1 leaf == parent
                   0, 0, 0, 0,
                   BranchPool{T,N,V}(branchpool_capacity, branch_capacity),
                   TwigPool{T,N,V}(branchpool_capacity, branch_capacity),
                   leafpool)
    end
end

RTree{T,N}(::Type{K}, ::Type{V}; kwargs...) where {T,N,K,V} =
    RTree{T,N,SpatialElem{T,N,K,V}}(; kwargs...)

mbrtype(::Type{<:RTree{T,N}}) where {T,N} = Rect{T,N}
mbrtype(tree::RTree) = mbrtype(typeof(tree))
regiontype(R::Type{<:RTree}) = mbrtype(R)
leaftype(::Type{RTree{T,N,V}}) where {T,N,V} = Leaf{T,N,V}
leaftype(tree::RTree) = leaftype(typeof(tree))
branchtype(::Type{RTree{T,N,V}}) where {T,N,V} = Branch{T,N,V}
branchtype(tree::RTree) = branchtype(typeof(tree))
nodetype(::Type{RTree{T,N,V}}) where {T,N,V} = Node{T,N,V} # not a concrete type
nodetype(tree::RTree) = nodetype(typeof(tree))

check_eltype_rtree(el::Any, tree::RTree) =
    check_eltype_rtree(typeof(el), leaftype(tree))

height(tree::RTree) = level(tree.root) + 1
variant(tree::RTree) = tree.variant
mbr(tree::RTree) = mbr(tree.root)

capacity(::T, tree::RTree) where T<:Node = capacity(T, tree)
capacity(::Type{<:Leaf}, tree::RTree) = node_capacity(tree.leafpool)
capacity(::Type{<:Branch}, tree::RTree) = node_capacity(tree.branchpool)

function acquire(tree::RTree, ::Type{<:Leaf}, lev::Int = 0, br::Rect = empty(mbrtype(tree)))
    @assert lev == 0
    leaf = acquire!(tree.leafpool)
    @assert isempty(leaf)
    leaf.mbr = br
    return leaf
end

function release(tree::RTree, node::Leaf)
    node.parent = nothing # don't refer to the parent (so it could be GCed)
    empty!(node.children) # don't refer to the children (so they could be GCed)
    release!(tree.leafpool, node)
end

function acquire(tree::RTree, ::Type{<:Branch}, lev::Int, br::Rect = empty(regiontype(tree)))
    node = acquire!(lev == 1 ? tree.twigpool : tree.branchpool)
    @assert isempty(node)
    node.mbr = br
    node.level = lev
    return node
end

function release(tree::RTree, node::Branch)
    empty!(node.children) # don't refer to the children (so they could be GCed)
    release!(eltype(node.children) <: Leaf ? tree.twigpool : tree.branchpool, node)
end

# low-level insertion of a child into the back of the list of the node children
# `child` parent is set to `node`, MBR of `node` is extended with child MBR
# no subtree balancing or anything like that
function _attach!(node::Node, child::Any, tree::RTree; force::Bool = false)
    force || (length(node) < capacity(node, tree)) ||
        throw(SpatialIndexException("Node is full ($(capacity(node, tree)) items)"))
    if child isa Node
        (level(node) == level(child) + 1) ||
            throw(SpatialIndexException("Parent ($(level(node)) and child ($(level(child))) levels don't match"))
        child.parent = node
    end
    push!(node.children, child)
    node.mbr = length(node) > 1 ? combine(node.mbr, mbr(child)) : mbr(child)
    return node
end

# low-level removal of a child node from the specified position, no subtree balancing
# return true if node MBR was updated
function _detach!(node::Node, pos::Integer, tree::RTree; updatembr = tree.tight_mbrs)
    if length(node) > 1
        delmbr = mbr(node[pos])   # cache it, since its required for "touches" later
        # copy the last child into i-th position and pop it
        # this should be more efficient that deleting the i-th node
        if pos < length(node)
            node.children[pos] = pop!(node.children)
        else
            pop!(node.children)
        end
        if updatembr && touches(mbr(node), delmbr)
            syncmbr!(node)
            return true
        end
    else # node becomes empty
        @assert pos == 1
        pop!(node.children)
        node.mbr = empty(mbrtype(node)) # reset br
        return true
    end
    return false
end

# The context for R-tree node insertion/balancing.
struct SubtreeContext{T,N,V,O}
    tree::RTree{T,N,V}
    level_overflow::O       # nothing or an indicator of the R-tree level overflow

    function SubtreeContext(tree::RTree{T,N,V}) where {T,N,V}
        O = variant(tree) == RTreeStar ? BitVector : Nothing
        new{T,N,V,O}(tree, O == BitVector ? falses(height(tree)) : nothing)
    end
end

isoverflow(con::SubtreeContext, lev::Int) = con.level_overflow[lev+1]
setoverflow!(con::SubtreeContext, lev::Int) = con.level_overflow[lev+1] = true
