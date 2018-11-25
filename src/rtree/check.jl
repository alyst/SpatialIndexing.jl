nelements(node::Leaf) = length(node)
nelements(node::Branch) = sum(nelements, children(node))

function check(node::Node{T,N,V}, tree::RTree{T,N,V},
               nnodes_perlevel::AbstractVector{Int}, ids::Union{Set, Nothing}) where {T,N,V}
    nnodes_perlevel[level(node) + 1] += 1
    for (i, child) in enumerate(children(node))
        if node isa Branch
            child.parent === node ||
                throw(SpatialIndexException("Node (lev=$(level(node)) len=$(length(node))): child #$i references different parent"))
            level(child) + 1 == level(node) ||
                throw(SpatialIndexException("Node (lev=$(level(node)) len=$(length(node))): child #$i level ($(level(child))) doesn't match parent's ($(level(node)))"))
            check(child, tree, nnodes_perlevel, ids)
        elseif idtrait(typeof(child)) !== HasNoID
            cid = id(child)
            (cid ∈ ids) && throw(SpatialIndexException("Duplicate data id=$cid"))
            push!(ids, cid)
        end
    end
    br = !isempty(node) ? mapreduce(mbr, combine, node.children) : empty(Rect{T,N})
    if isempty(node)
        isequal(br, mbr(node)) ||
            throw(SpatialIndexException("Node (lev=$(level(node)) len=$(length(node))) MBR ($(mbr(node))) should be empty"))
    elseif tree.tight_mbrs && br != mbr(node) || !in(br, mbr(node))
        throw(SpatialIndexException("Node (lev=$(level(node)) len=$(length(node))) MBR ($(mbr(node))) doesn't match its children MBR ($(br))"))
    end
    return true
end

# check the tree:
# * the parent
# * the height
# * the number of elements
# * the number of nodes per level
# * element id uniqueness
# * MBRs of the nodes
function check(tree::RTree)
    tree.root.parent === nothing ||
        throw(SpatialIndexException("RTree: root has a parent"))
    level(tree.root) + 1 == height(tree) ||
        throw(SpatialIndexException("RTree: root level ($(level(tree.root))) doesn't match tree height ($(height(tree)))"))
    length(tree.nnodes_perlevel) == height(tree) ||
        throw(SpatialIndexException("RTree: nnodes_perlevel length ($(length(tree.nnodes_perlevel))) doesn't match tree height ($(height(tree)))"))
    ids = idtype(tree) !== Union{} ? Set{idtype(tree)}() : nothing
    nnodes_perlevel = zeros(Int, height(tree))
    check(tree.root, tree, nnodes_perlevel, ids)
    nelms = nelements(tree.root)
    length(tree) == nelms ||
        throw(SpatialIndexException("RTree: actual ($(nelms)) and reported ($(length(tree))) number of elements do not match"))
    ids === nothing || length(tree) == length(ids) ||
        throw(SpatialIndexException("RTree: the number of ids ($(length(ids))) doesn't match the reported number of elements ($(length(tree)))"))
    tree.nnodes_perlevel == nnodes_perlevel ||
        throw(SpatialIndexException("RTree: actual ($(nnodes_perlevel)) and reported ($(tree.nnodes_perlevel)) number of nodes per level do not match"))
    return true
end
