# default branch fill strategy
# returns the tuple:
# the number of 1st dim slices and the number of child nodes in each slice
function omt_branch_fill(tree::RTree; fill_factor::Number=1.0) where T<:Node
    @assert 0.0 <= fill_factor <= 1.0
    maxlen = ceil(Int, capacity(Branch, tree) * fill_factor)
    slicelen = ceil(Int, 0.75*sqrt(maxlen))
    nslices = fld1(maxlen, slicelen)
    return (nslices, slicelen)
end

"""
Bulk-load the `data` into `tree`.
"""
function load!(tree::RTree{T,N,V}, data::Any;
               getmbr = first,
               getid = idtrait(V) === HasNoID ? x -> nothing : x -> x[2],
               getval = idtrait(V) === HasNoID ? x -> x[2] : x -> x[3],
               method::Symbol = :OMT,
               leaf_fill::Integer = capacity(Leaf, tree),
               branch_fill::Tuple{Integer, Integer} = omt_branch_fill(tree)) where {T, N, V}
    # load from arbitrary data into Elem-backed R-tree
    isempty(tree) || throw(ArgumentError("Cannot bulk-load into non-empty tree"))
    if isempty(data)
        @warn "No bulk-load data provided"
        return tree
    end
    # FIXME kwargs... doesn't work for some reason
    return load!(tree, V[V(getmbr(x), getid(x), getval(x)) for x in data],
                 method=method, leaf_fill=leaf_fill, branch_fill=branch_fill)
end

function load!(tree::RTree{T,N,V}, data::AbstractVector{V};
               method::Symbol = :OMT,
               leaf_fill::Integer = capacity(Leaf, tree),
               branch_fill::Tuple{Integer, Integer} = omt_branch_fill(tree)) where {T,N,V}
    isempty(tree) || throw(ArgumentError("Cannot bulk-load into non-empty tree"))
    if isempty(data)
        @warn "No bulk-load data provided"
        return tree
    end
    0 < leaf_fill <= capacity(Leaf, tree) ||
        throw(ArgumentError("Leaf fill should be positive and not exceed leaf capacity"))
    prod(branch_fill) > 1 || throw(ArgumentError("Branch fill should be > 1"))
    if method == :OMT
        return load_omt!(tree, data, leaf_fill=leaf_fill, branch_fill=branch_fill)
    else
        throw(ArgumentError("Unsupported bulk-loadaing method $method"))
    end
end

# See "OMT: Overlap Minimizing Top-down Bulk Loading Algorithm for R-tree"
# by Taewon Lee and Sukho Lee, http://ceur-ws.org/Vol-74/files/FORUM_18.pdf
# Modified to take into account that leaf and branch capacity is different
# and allow branch slices at any tree level (specified by `branch_fill`)
function load_omt!(tree::RTree, elems::AbstractVector;
                   leaf_fill::Integer = capacity(Leaf, tree),
                   branch_fill::Tuple{Integer, Integer} = omt_branch_fill(tree))
    # calculate the tree properties
    nbranch_subtrees = min(branch_fill[1] * branch_fill[2], capacity(Leaf, tree))
    height = ceil(Int, max(log(length(elems)) - log(leaf_fill), 0.0) / log(nbranch_subtrees)) + 1 # 1 for leaves
    maxelems_subtree = leaf_fill*nbranch_subtrees^(height-2) # max elements in each root subtree
    nroot_subtrees = fld1(length(elems), maxelems_subtree)
    nroot_slices = floor(Int, sqrt(nroot_subtrees))
    #@show length(elems) leaf_fill branch_fill height maxelems_subtree nsubtrees nslices
    # sort by the center of the first dim
    dim = mod1(height, ndims(tree)) # root level dimension
    elems_sorted = sort(elems, by = elem -> mbr(elem).low[dim] + mbr(elem).high[dim])
    resize!(tree.nnodes_perlevel, height)
    fill!(tree.nnodes_perlevel, 0)
    tree.root = omt_subtree(elems_sorted, tree, height,
                            nroot_slices, nroot_subtrees, leaf_fill, branch_fill)
    tree.nnodes_perlevel[end] = 1
    return tree
end

# put the elems into the subtree of height lev
# the root node of the subtree should have nsubtree children
# organized into nslices slices along the current elems order
# recursively calls itself
function omt_subtree(elems::AbstractVector, tree::RTree,
                     lev::Integer, nslices::Integer, nsubtrees::Integer,
                     leaf_fill::Integer, branch_fill::Tuple{Integer, Integer})
    @debug "omt_subtree(): lev=$lev nslices=$nslices nsubtrees=$nsubtrees nelems=$(length(elems))"
    @assert lev > 0
    if length(elems) <= nsubtrees
        # if fewer elements than the number of subtrees,
        # then all elements should be put into single leaf
        if lev == 1 # create a Leaf and attach all elements
            node = acquire(tree, Leaf, lev)
            for elem in elems
                _attach!(node, elem, tree)
            end
            tree.nelems += length(elems)
            return node
        else # not leaf level yet, create a branch with a single child
            child = omt_subtree(elems, tree, lev - 1, nslices, nsubtrees,
                                leaf_fill, branch_fill)
            node = _attach!(acquire(tree, Branch, lev), child, tree)
            tree.nnodes_perlevel[lev-1] += 1
            return node
        end
    end

    # subtree root
    @assert lev > 1
    node = acquire(tree, Branch, lev)

    # create subtrees
    nelems_subtree = fld1(length(elems), nsubtrees) # actual fill of the subtree
    nelems_slice = nelems_subtree * fld1(nsubtrees, nslices)
    dim = mod1(lev, ndims(tree)) # cycle sorting dimensions through tree levels
    for i in 1:nelems_slice:length(elems) # slice using the external elements order (sort-dim of the level above)
        # sort the elements of the slice by the current dim
        elems_slice = sort(view(elems, i:min(i+nelems_slice-1, length(elems))),
                           by = elem -> mbr(elem).low[dim] + mbr(elem).high[dim])
        # create slice subtrees and attach to the node
        for j in 1:nelems_subtree:length(elems_slice)
            subtree_elems_range = j:min(j+nelems_subtree-1, length(elems_slice))
            child = omt_subtree(view(elems_slice, subtree_elems_range), tree,
                                lev - 1, branch_fill[1],
                                min(branch_fill[1] * branch_fill[2], capacity(Branch, tree)),
                                leaf_fill, branch_fill)
            _attach!(node, child, tree)
            tree.nnodes_perlevel[lev - 1] += 1
        end
    end
    return node
end
