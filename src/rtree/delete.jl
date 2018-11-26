"""
    delete!(tree::RTree, pt::Point, [id])
    delete!(tree::RTree, br::Rect, [id])

Deletes the value identified by `br` bounding box (or point `pt`) and
the `id` (if tree elements support `HasID` trait) from the `tree`.
"""
function Base.delete!(tree::RTree{T,N}, br::Rect{T,N}, id::Any = nothing) where {T,N}
    leafx = findleaf(tree, br, id)
    leafx === nothing && __spatial_keyerror(eltype(tree), br, id)
    _detach!(leafx[1], leafx[2], tree)
    tmpdetached = Vector{nodetype(tree)}() # FIXME use pool, Union{Leaf,Branch} ?
    _condense!(leafx[1], tree, tmpdetached)
    _reinsert!(tree, tmpdetached)
    tree.nelems -= 1
    tree.nelem_deletions += 1
    return tree
end

Base.delete!(tree::RTree{T,N}, pt::Point{T,N}, id::Any = nothing) where {T,N} =
    delete!(tree, Rect(pt), id)

# deletes the subtree with the `node` root from the `tree`
# (does not update the parent MBR)
function delete_subtree!(tree::RTree, node::Node)
    #@debug "delete_subtree(): lev=$(level(node))"
    _delete_subtree!(node, tree, level(node))
end

function _delete_subtree!(node::Node, tree::RTree, height::Integer)
    if node isa Branch
        for child in children(node)
            _delete_subtree!(child, tree, height)
        end
    elseif node isa Leaf
        tree.nelems -= length(node)
        tree.nelem_deletions += length(node)
    end
    #@debug "_delete_subtree(): empty children ($(length(node.children)) objs) of node lev=$(level(node)) parent=$(hasparent(node))"
    empty!(node.children)
    if hasparent(node) # remove the node
        tree.nnodes_perlevel[level(node)] -= 1 # don't count this node anymore
        if level(node) == height # this is the root of the deleted subtree, detach it from the parent
            _detach!(parent(node), pos_in_parent(node), tree)
            # FIXME propagate parent mbr updates (if required)
        end
        release(tree, node)
    else # root node just stays empty
        node.mbr = empty(mbrtype(node))
    end
    #@debug "_delete_subtree(): done node lev=$(level(node)) parent=$(hasparent(node))"
    return nothing
end

# FIXME replace `region` with a more generic `filter`
"""
    subtract!(tree::RTree, reg::Region)

Subtracts the `region` from the `tree`, i.e. removes all elements within
`region`.
"""
function subtract!(tree::RTree{T,N}, reg::Region{T,N}) where {T,N}
    #@debug "subtract!(): region=$(reg)"
    isempty(tree) && return tree

    tmpdetached = Vector{nodetype(tree)}()
    status = _subtract!(tree.root, 0, reg, tree, tmpdetached)
    #@debug "subtract!(): status=$(status) tmpdetached=$(length(tmpdetached))"
    if status == 1 # tree changed, but not removed
        _condense!(tree.root, tree, tmpdetached) # try to condense the root
    elseif status == 2 # whole tree removed (contained in reg)
        # should have no items and nothing to reinsert
        @assert isempty(tree)
        @assert isempty(tmpdetached)
    end
    _reinsert!(tree, tmpdetached)
    return tree
end

function _subtract!(node::Node, node_ix::Int, reg::Region, tree::RTree,
                    tmpdetached::AbstractVector{<:Node})
    nodembr = mbr(node)
    #@debug "_subtract!(): lev=$(level(node)) i=$node_ix len=$(length(node))"
    # TODO juxtaposition() method that combines in() and intersects()?
    if in(nodembr, reg)
        #@debug "_subtract!(): delete subtree lev=$(level(node)) i=$node_ix len=$(length(node))"
        delete_subtree!(tree, node)
        return 2 # node removed
    elseif intersects(nodembr, reg)
        mbr_dirty = false
        i = 1
        while i <= length(node)
            #@debug "1: lev=$(level(node)) i=$i len=$(length(node))"
            child = node[i]
            oldmbr = mbr(child)
            if node isa Branch
                child_res = _subtract!(child, i, reg, tree, tmpdetached)
            else
                @assert node isa Leaf
                if in(oldmbr, reg)
                    #@debug "_subtract!(): detach elem lev=$(level(node)) i=$i len=$(length(node))"
                    _detach!(node, i, tree, updatembr = false) # don't update MBR, we do it later
                    tree.nelems -= 1
                    tree.nelem_deletions += 1
                    child_res = 2
                    # FIXME intersect is not considered (would be important if elem mbr is not a point), should be handled by `filter`
                else
                    child_res = 0
                end
            end
            if !mbr_dirty && (child_res != 0) && tree.tight_mbrs
                mbr_dirty = touches(nodembr, oldmbr)
            end
            if child_res != 2 # don't increment if i-th node removed
                i += 1
            end
        end
        if hasparent(node) && length(node) < floor(Int, tree.reinsert_factor * capacity(node, tree))
            _detach!(parent(node), node_ix, tree)
            tree.nnodes_perlevel[level(node)] -= 1
            push!(tmpdetached, node)
            return 2 # node removed
        elseif mbr_dirty
            syncmbr!(node)
            return 1 # mbr changed
        else
            return 0 # no mbr change
        end
    else
        return 0 # no change
    end
end
