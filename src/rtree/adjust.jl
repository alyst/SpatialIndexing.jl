# condense the tree after the child was removed from the node
# accumulate the nodes that are temporary detached and need to be reinserted in `tmpdetached`
function _condense!(node::Node, tree::RTree, tmpdetached::AbstractVector{<:Node})
    min_load = floor(Int, capacity(node, tree) * tree.fill_factor)
    @debug "_condense!() lev=$(level(node)) len=$(length(node)) to_reinsert=$(length(to_reinsert)) min_load=$min_load"

    if node === tree.root
        # eliminate root if it has only one child.
        if level(node) > 1
            if length(node) == 1
                @debug "_condense!(): shrinking single-child root"
                tree.root = node[1]
                tree.root.parent = nothing
                release(tree, node)
                pop!(tree.nnodes_perlevel) # tree became 1 level shorter
            elseif isempty(node)
                @debug "_condense!(): resetting empty root"
                # reset the root to a node of minimal level to accomodate the children of tmpdetached (or leaf)
                root_level = mapreduce(level, max, tmpdetached, init=1)
                tree.root = acquire(tree, root_level == 1 ? Leaf : Branch, root_level)
                resize!(tree.nnodes_perlevel, root_level)
                fill!(tree.nnodes_perlevel, 0)
                tree.nnodes_perlevel[end] = 1
            else
                # due to data removal.
                tree.tight_mbrs && syncmbr!(node)
            end
        else
            # due to data removal.
            tree.tight_mbrs && syncmbr!(node)
        end
        return node
    else
        # find the entry in the parent, that points to this node
        node_ix = pos_in_parent(node)
        @assert node_ix !== nothing

        if length(node) < min_load
            # used space less than the minimum
            # 1. eliminate node entry from the parent. deleteEntry will fix the parent's MBR.
            _detach!(parent, node_ix, tree)
            # 2. add this node to the stack in order to reinsert its entries.
            push!(tmpdetached, node)
        else
            # global recalculation necessary since the MBR can only shrink in size,
            # due to data removal.
            tree.tight_mbrs && syncmbr(parent(node))
        end

        return _condense!(parent(node), tree, tmpdetached)
    end
end

# reinsert the *children* of the `detached` nodes back to the tree
# the `detached` node themselves are released back to the pool
function _reinsert!(tree::RTree, detached::AbstractVector{<:Node})
    @debug "_reinsert!(): reinsert children of $(length(detached)) detached nodes"
    isempty(detached) && return tree
    con = SubtreeContext(tree)
    for node in detached
        for child in children(node)
            _insert!(con, child, level(node))
            tree.nnode_reinsertions += 1
        end
        release(tree, node)
    end
    empty!(detached)
    return tree
end

# update the MBR of node if required due to the child_ix MBR changing from
# child_oldmbr to its current state or if forced,
# propagate the MBR update to the higher tree levels
# return true if MBR was updated
function _updatembr!(node::Branch, child_ix::Integer, child_oldmbr::Rect,
                     con::SubtreeContext; force::Bool = false)
    @debug "_updatembr!() tree_height=$(height(con.tree)) tree_len=$(length(con.tree)) node_lev=$(level(node)) node_len=$(length(node)) force=$(force) child_ix=$(child_ix) oflow=$(isoverflow(con, level(node)))"
    child = node[child_ix]
    # MBR needs recalculation if either:
    #   1. the NEW child MBR is not contained.
    #   2. the OLD child MBR is touching.
    node_oldmbr = mbr(node)
    mbr_dirty = force || !in(mbr(child), node_oldmbr) ||
                (con.tree.tight_mbrs && touches(node_oldmbr, child_oldmbr))
    mbr_dirty && syncmbr!(node)

    if mbr_dirty && hasparent(node)
        _updatembr!(parent(node), pos_in_parent(node), node_oldmbr, con, force=force)
    end
    return mbr_dirty
end

# replace the node (with oldmbr MBR) with the newnode
# return true if node's parent MBR update was necessary
function _replace!(node::Node, newnode::Node, oldmbr::Rect,
                   con::SubtreeContext)
    @debug "_replace!() lev=$(level(node)) len=$(length(node)) newlen=$(length(newnode))"
    # find an entry pointing to the old child
    @assert level(node) == level(newnode)
    par = parent(node)
    @assert par !== nothing
    node_ix = pos_in_parent(node)
    @assert node_ix !== nothing
    par[node_ix] = newnode
    return _updatembr!(par, node_ix, oldmbr, con)
end

# replace node with n1 and n2
# return true if node's parent MBR update or any other tree
# restructurings were necessary
function _replace!(node::Node, n1::Node, n2::Node,
                   oldmbr::Rect, con::SubtreeContext)
    @debug "_replace!() lev=$(level(node)) len=$(length(node)) newlens=($(length(n1)), $(length(n2)))"
    @assert level(node) == level(n1) == level(n2)
    # find entry pointing to old node
    par = parent(node)
    @assert par !== nothing
    node_ix = pos_in_parent(node)
    @assert node_ix !== nothing

    # MBR needs recalculation if either:
    #   1. MBRs of n1 and n2 are not contained in parent
    #   2. the OLD node MBR is touching parent (FIXME should it be recalced though?)
    par_oldmbr = mbr(par)
    mbr_dirty = !in(mbr(n1), par_oldmbr) || !in(mbr(n2), par_oldmbr) ||
                (con.tree.tight_mbrs && touches(par_oldmbr, oldmbr))
    # replace the node with n1
    par[node_ix] = n1
    mbr_dirty && syncmbr!(par)

    con.tree.nnode_splits += 1
    con.tree.nnodes_perlevel[level(node)] += 1

    # No registering necessary. insert!() will write the node if needed.
    #_register!(tree, node)
    adjusted = _insert!(par, n2, con)

    # if n2 is contained in the node and there was no split or reinsert,
    # we need to adjust only if recalculation took place.
    # In all other cases insertData above took care of adjustment.
    if !adjusted && mbr_dirty && hasparent(par)
        _updatembr!(parent(par), pos_in_parent(par), par_oldmbr, con)
    end
    return adjusted || mbr_dirty
end
