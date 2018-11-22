"""
    insert!(tree::RTree, key, pt::Point, val)
    insert!(tree::RTree, key, br::Rect, val)
    insert!(tree::RTree, elem::SpatialElem)

Inserts `val` value identified by `key` and `br` bounding box (or point `pt`) into the `tree`.
"""
function Base.insert!(tree::RTree, el::Any)
    check_eltype_rtree(el, tree)
    _insert!(tree, el, 0)
    tree.nelems += 1
    tree.nelem_insertions += 1
    return tree
end

Base.insert!(tree::RTree{T,N,SpatialElem{T,N,K,V}}, br::Rect{T,N}, key::Any, val::Any) where {T,N,K,V} =
    insert!(tree, SpatialElem{T,N,K,V}(br, key, val))

Base.insert!(tree::RTree{T,N,SpatialElem{T,N,K,V}},
             pt::Point{T,N}, key::Any, val::Any) where {T,N,K,V} =
    insert!(tree, SpatialElem{T,N,K,V}(Rect(pt), key, val))

Base.insert!(tree::RTree{T,N,SpatialElem{T,N,Nothing,V}},
             br::Union{Rect{T,N}, Point{T,N}}, val::Any) where {T,N,V} =
    insert!(tree, nothing, br, val)

# inserts the child into the node rebalancing R-tree as needed
# returns true if MBR was changed and upstream R-tree readjustments are required
function _insert!(node::Node, child::Any, con::SubtreeContext)
    if length(node) < capacity(node, con.tree) # underfull node
        oldmbr = mbr(node) # this has to happen before _attach!() modifies node's MBR
        _attach!(node, child, con.tree)
        if !in(mbr(child), oldmbr) && hasparent(node)
            #@debug "_insert!(): mbr change"
            _updatembr!(parent(node), pos_in_parent(node), oldmbr, con)
            return true
        else
            return false
        end
    elseif variant(con.tree) == RTreeStar && (node !== con.tree.root) && !isoverflow(con, level(node))
        return _insert!_fullnode_rstar(node, child, con)
    else
        return _insert!_fullnode(node, child, con)
    end
end

function _insert!(con::SubtreeContext, node::Any, lev::Integer)
    #@debug "_insert!(con): len=$(length(con.tree)) node_lev=$(node isa Node ? level(node) : "data") node_len=$(node isa Node ? length(node) : "data")"
    par = choose_subtree(con.tree.root, lev, mbr(node), con.tree)
    @assert level(par) == lev
    _insert!(par, node, con)
end

_insert!(tree::RTree, node::Any, lev::Int) =
    _insert!(SubtreeContext(tree), node, lev)

# index of the child that would have the least MBR enlargement by br
function find_least_enlargement(node::Node, br::Rect)
    min_enl = Inf
    best_ix = 0
    best_area = Inf
    for (i, child) in enumerate(children(node))
        child_area = area(mbr(child))
        enl = combined_area(mbr(child), br) - child_area
        if enl < min_enl || (min_enl == enl && child_area < best_area)
            min_enl = enl
            best_area = child_area
            best_ix = i
        end
    end
    return best_ix
end

# index of the child that would have the least MBR overlap with br
function find_least_overlap(node::Node, br::Rect, tree::RTree)
    min_enl = Inf
    best_area = Inf
    best_ix = 0
    enls = Vector{Float64}(undef, length(node)) # FIXME use pool?
    # find combined region and enlargement of every entry and store it.
    for (i, child) in enumerate(children(node))
        child_area = area(mbr(child))
        enls[i] = enl = combined_area(br, mbr(child)) - child_area
        if enl < min_enl || (enl == min_enl && child_area < best_area)
            min_enl = enl
            best_area = child_area
            best_ix = i
        end
    end

    abs(min_enl) <= eps(min_enl) && return best_ix # no enlargement

    if length(node) > tree.nearmin_overlap
        # sort entries in increasing order of enlargement
        enl_order = sortperm(enls) # FIXME reuse from pool
        @assert min_enl == enls[enl_order[1]] <= enls[enl_order[end]]
        niter = tree.nearmin_overlap
    else
        niter = length(node)
        enl_order = collect(1:niter) # FIXME reuse from pool
    end

    # calculate overlap of most important original entries (near minimum overlap cost).
    min_delta_olap = Inf
    best_ix = 0
    best_area = Inf
    best_enl = Inf
    for i in 1:niter
        node_ix = enl_order[i]
        node_mbr = mbr(node[node_ix])
        comb_mbr = combine(node_mbr, br)
        delta_olap = 0.0

        for (j, child) in enumerate(children(node))
            if node_ix != j
                olap = overlap_area(comb_mbr, mbr(child))
                if olap > 0.0
                    delta_olap += olap - overlap_area(node_mbr, mbr(child))
                end
            end
        end

        node_area = area(node_mbr)
        # keep the one with the least, in the order:
        #   delta overlap, enlargement, area
        if (delta_olap < min_delta_olap) || (delta_olap == min_delta_olap &&
           (enls[node_ix] < best_enl || (enls[node_ix] == best_enl &&
            node_area < best_area)))
            min_delta_olap = delta_olap
            best_ix = node_ix
            best_area = node_area
            best_enl = enls[node_ix]
        end
    end
    return best_ix
end

# choose the node (for inserting the given mbr at specified level)
function choose_subtree(node::Node, lev::Int, br::Rect, tree::RTree)
    while level(node) != lev
        if variant(tree) == RTreeLinear || variant(tree) == RTreeQuadratic
            child_ix = find_least_enlargement(node, br)
        elseif variant(tree) == RTreeStar
            child_ix = level(node) == 1 ?
                    find_least_overlap(node, br, tree) : # for the node pointing to leaves
                    find_least_enlargement(node, br)
        else
            throw(SpatialIndexException("RTree variant not supported"))
        end
        node = node[child_ix]
    end
    @assert level(node) == lev
    return node
end

# R*-tree-specific insertion into the full node
# splits the node and rebalances the R-tree
function _insert!_fullnode_rstar(node::Node, child::Any, con::SubtreeContext)
    #@debug "_insert!_fullnode_rstar(): lev=$(level(node)) len=$(length(node))"
    @assert isa(con.level_overflow, BitVector)
    @assert hasparent(node)

    setoverflow!(con, level(node)) # the R* insert could only be done once per level
    oldmbr = mbr(node)
    _attach!(node, child, con.tree, force=true) # force-push the child exceeding the capacity

    node_center = center(mbr(node))
    # calculate relative distance of every entry from the node MBR (ignore square root.)
    children_dist = sqrdistance.(Ref(node_center), center.(mbr.(children(node)))) # FIXME reuse array?
    # return the indices sorted by the increasing order of distances.
    # Since children_dist is sorted in ascending order
    # it would match the reinsertion order suggested in the paper.
    children_order = sortperm(children_dist) # FIXME reuse array?
    nkeep = length(node) - floor(Int, length(node) * con.tree.reinsert_factor)

    newnode = acquire(con.tree, typeof(node), level(node))   # the replacement node
    for i in 1:nkeep
        _attach!(newnode, node[children_order[i]], con.tree)
    end

    # Divertion from R*-Tree algorithm here. First update
    # the path to the root, then start reinserts, to avoid complicated handling
    # of changes to the same node from multiple insertions.
    #@debug "_insert!_fullnode_rstar(): replace the node"
    _replace!(node, newnode, oldmbr, con) # replace the node (so (nkeep+1):... children are detached from the tree)
    #@debug "_insert!_fullnode_rstar(): reinsert $(length((nkeep+1):length(children_order))) children of lev=$(level(node)) node"
    # don't use par during/after reinsertion, it might have changed
    for i in (nkeep+1):length(children_order)
        child = node[children_order[i]]
        #@debug "_insert!_fullnode_rstar(): reinsert #$i child (type=$(typeof(child)) lev=$(child isa Node ? level(child) : "data") parent=$(isa(child, Node) && hasparent(child)))"
        _insert!(con, child, level(node))
        con.tree.nnode_reinsertions += 1
    end

    release(con.tree, node) # return node to the pool
    return true
end

# R-tree-specific insertion into the full node
# splits the node and rebalances the R-tree
function _insert!_fullnode(node::Node, child::Any, con::SubtreeContext)
    #@debug "_insert!_fullnode(): lev=$(level(node)) len=$(length(node)) parent=$(hasparent(node))"
    oldmbr = mbr(node)
    _attach!(node, child, con.tree, force=true) # temporary overflow the node capacity
    n1, n2 = _split!(node, con.tree)

    if node === con.tree.root
        #@debug "_insert!_fullnode(): newroot tree_height=$(height(con.tree)) tree_len=$(length(con.tree))"
        # the node must be the root, make the new root and attach n1, n2 to it
        newroot = acquire(con.tree, Branch, level(node) + 1)
        _attach!(newroot, n1, con.tree)
        _attach!(newroot, n2, con.tree)
        # tree grows
        con.tree.root = newroot
        con.tree.nnodes_perlevel[level(node)+1] = 2
        push!(con.tree.nnodes_perlevel, 1)
    else # non-root, n1 and n2 replace the node in its parent
        _replace!(node, n1, n2, oldmbr, con)
    end
    return true
end
