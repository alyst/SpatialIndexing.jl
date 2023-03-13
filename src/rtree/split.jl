# select the two seeds maximizing separation, O(N)
function _splitseeds_linear(node::Node)
    @assert length(node) > 1

    seed1 = seed2 = 0
    max_separation = -Inf
    for n in 1:ndims(node)
        argmax_low = argmin_high = 1
        min_low = max_low = mbr(node[1]).low[n]
        min_high = max_high = mbr(node[1]).high[n]

        for i in 2:length(node)
            @inbounds lowi = mbr(node[i]).low[n]
            @inbounds highi = mbr(node[i]).high[n]
            if lowi > max_low
                max_low = lowi
                argmax_low = i
            elseif lowi < min_low
                min_low = lowi
            end
            if highi < min_high
                min_high = highi
                argmin_high = i
            elseif highi > max_high
                max_high = highi
            end
        end

        sepn = max_low - min_high
        if max_high > min_low
            sepn /= (max_high - min_low)
        end

        if sepn > max_separation
            max_separation = sepn
            seed1 = argmin_high
            seed2 = argmax_low
        end
    end

    @assert seed1 > 0 && seed2 > 0
    if seed1 == seed2 # avoid collisions
        if seed2 == 1
            seed2 += 1
        else
            seed2 -= 1
        end
    end
    return seed1, seed2
end

# select the two seeds maximizing inefficiency, O(N²)
function _splitseeds_quadratic(node::Node)
    @assert length(node) > 1

    seed1 = 0
    seed2 = 0
    inefficiency = -Inf
    # for each pair of Regions (account for overflow Region too!)
    for i in eachindex(node)
        @inbounds mbri = mbr(node[i])
        areai = area(mbri)

        for j in (i+1):length(node)
            # get the combined MBR of those two entries.
            @inbounds mbrj = mbr(node[j])
            @inbounds r = combine(mbri, mbrj)

            # find the inefficiency of grouping these entries together.
            d = area(r) - areai - area(mbrj)
            if d > inefficiency
                inefficiency = d;
                seed1 = i
                seed2 = j
            end
        end
    end
    return seed1, seed2
end

# pick the two seeds (children indices) for splitting the node into two
function _splitseeds(node::Node, tree::RTree)
    length(node) > 1 ||
        throw(SpatialIndexError("Cannot split the node with less than 2 children"))

    if variant(tree) == RTreeLinear || variant(tree) == RTreeStar
        return _splitseeds_linear(node)
    elseif variant(tree) == RTreeQuadratic
        return _splitseeds_quadratic(node)
    else
        throw(SpatialIndexException("RTree variant not supported"))
    end
end

# R-tree split node (select the node by minimal enlargment and area)
function _split!_rtree(node::Node, tree::RTree)
    nnodes_min = floor(Int, capacity(node, tree) * tree.fill_factor)
    #@debug "_split!_rtree(): lev=$(level(node)) len=$(length(node)) nmin=$nnodes_min"

    # initialize each group with the seed entries.
    seed1, seed2 = _splitseeds(node, tree)
    n1 = _attach!(acquire(tree, typeof(node), level(node)), node[seed1], tree)
    n2 = _attach!(acquire(tree, typeof(node), level(node)), node[seed2], tree)

    # use this mask for marking visited entries.
    available = trues(length(node))
    available[seed1] = available[seed2] = false
    navail = length(node) - 2
    while navail > 0
        if length(n1) + navail == nnodes_min || length(n2) + navail == nnodes_min
            # all remaining entries must be assigned to n1 or n2 to comply with minimum load requirement.
            sink = length(n1) + navail == nnodes_min ? n1 : n2
            i = findfirst(available)
            while i !== nothing
                _attach!(sink, node[i], tree)
                available[i] = false
                i = findnext(available, i+1)
            end
            navail = 0
        else
            # For all remaining entries compute the difference of the cost of grouping an
            # entry in either group. When done, choose the entry that yielded the maximum
            # difference. In case of linear split, select any entry (e.g. the first one.)
            area1 = area(mbr(n1))
            area2 = area(mbr(n2))

            max_cost_diff = -Inf
            sel_cost1 = sel_cost2 = NaN
            sel_child_ix = 0
            i = findfirst(available)
            while i !== nothing
                child_mbr = mbr(node[i])
                enl1 = area(combine(mbr(n1), child_mbr)) - area1
                enl2 = area(combine(mbr(n2), child_mbr)) - area2
                enl_diff = abs(enl1 - enl2)

                if enl_diff > max_enl_diff
                    max_cost_diff = enl_diff
                    sel_enl1 = enl1
                    sel_enl2 = enl2
                    sel_child_ix = i
                    (variant(tree) != RTreeQuadratic) && break
                end
                i = findnext(available, i+1)
            end

            # select the node where we should add the new entry.
            sink = sel_enl1 < sel_enl2 || (sel_enl1 == sel_enl2 &&
                  (area1 < area2 || (area1 == area2 &&
                   length(n1) <= length(n2)))) ? n1 : n2
            available[sel_child_ix] = false
            navail -= 1
            _attach!(sink, node[sel_child_ix], tree)
        end
    end
    @assert !any(available)

    return n1, n2
end

# get low/high of a specific MBR dimension
# accessor functions to speedup sortperm!()
_mbr_low(node, dim::Integer) = @inbounds mbr(node).low[dim]
_mbr_high(node, dim::Integer) = @inbounds mbr(node).high[dim]

# R*-star node split
function _split!_rstar(node::Node, tree::RTree)
    nsplit = floor(Int, length(node) * tree.split_factor)
    nsplit_distr = length(node) - 2 * nsplit + 1
    @assert nsplit_distr > 0

    #@debug "_split!_rstar(): lev=$(level(node)) len=$(length(node)) nsplit=$nsplit nsplit_distr=$nsplit_distr"

    # find split_dim that minimizes the perimeter in all splits
    min_perim = Inf
    split_dim = 0
    loworder = zeros(Int, length(node))
    highorder = similar(loworder)
    use_low = false
    for dim in 1:ndims(node)
        sortperm!(loworder, children(node), by=Base.Fix2(_mbr_low, dim))
        sortperm!(highorder, children(node), by=Base.Fix2(_mbr_high, dim))

        # calculate the sum of perimiters for all splits
        low_perim = 0.0
        high_perim = 0.0
        for i in nsplit:(nsplit + nsplit_distr)
            @inbounds br_low1 = mapreduce(j -> @inbounds(mbr(node[j])), combine, view(loworder, 1:i))
            @inbounds br_high1 = mapreduce(j -> @inbounds(mbr(node[j])), combine, view(highorder, 1:i))
            @inbounds br_low2 = mapreduce(j -> @inbounds(mbr(node[j])), combine, view(loworder, (i+1):length(node)))
            @inbounds br_high2 = mapreduce(j -> @inbounds(mbr(node[j])), combine, view(highorder, (i+1):length(node)))

            low_perim += perimeter(br_low1) + perimeter(br_low2)
            high_perim += perimeter(br_high1) + perimeter(br_high2)
        end

        perim = min(low_perim, high_perim)
        if perim < min_perim
            min_perim = perim
            split_dim = dim
            use_low = low_perim < high_perim
        end
    end

    # final sorting
    selorder = use_low ?
        sortperm!(loworder, children(node), by=Base.Fix2(_mbr_low, split_dim)) :
        sortperm!(highorder, children(node), by=Base.Fix2(_mbr_high, split_dim))

    # find the best split point (minimizes split overlap and area)
    min_overlap = Inf
    best_area_sum = Inf
    split_ix = 0
    for i in nsplit:(nsplit + nsplit_distr)
        @inbounds br1 = mapreduce(i -> mbr(node[selorder[i]]), combine, 1:i)
        @inbounds br2 = mapreduce(i -> mbr(node[selorder[i]]), combine, (i+1):length(node))

        overlap = overlap_area(br1, br2)
        area_sum = area(br1) + area(br2)
        if overlap < min_overlap || (overlap == min_overlap &&
           area_sum < best_area_sum)
            split_ix = i
            min_overlap = overlap
            best_area_sum = area_sum
        end
    end

    n1 = acquire(tree, typeof(node), level(node))
    for i in 1:split_ix
        _attach!(n1, node[selorder[i]], tree)
    end
    n2 = acquire(tree, typeof(node), level(node))
    for i in (split_ix+1):length(selorder)
        _attach!(n2, node[selorder[i]], tree)
    end
    return n1, n2
end

# split the node into two
function _split!(node::Node, tree::RTree)
    #@debug "_split!(): lev=$(level(node)) len=$(length(node))"

    if variant(tree) == RTreeLinear || variant(tree) == RTreeQuadratic
        return _split!_rtree(node, tree)
    elseif variant(tree) == RTreeStar
        return _split!_rstar(node, tree)
    else
        throw(SpatialIndexError("RTree variant not supported"))
    end
end
