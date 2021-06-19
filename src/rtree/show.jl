function Base.show(io::IO, tree::RTree)
    print(io, typeof(tree))
    print(io, "(variant="); print(io, tree.variant)
    print(io, ", tight_mbrs="); print(io, tree.tight_mbrs)
    print(io, ", nearmin_overlap="); print(io, tree.nearmin_overlap)
    print(io, ", fill_factor="); print(io, tree.fill_factor)
    print(io, ", split_factor="); print(io, tree.split_factor)
    print(io, ", reinsert_factor="); print(io, tree.reinsert_factor)

    print(io, ", leaf_capacity="); print(io, capacity(Leaf, tree))
    print(io, ", branch_capacity="); print(io, capacity(Branch, tree)); println(io, ")")
    print(io, length(tree)); print(io, " element(s) in "); print(io, height(tree));
    print(io, " level(s) ("); print(io, join(reverse(tree.nnodes_perlevel), ", ")); println(io, " node(s) per level):")
    _show(io, tree.root, recurse=true, indent=1)
end

function _show(io::IO, node::Node; recurse::Bool=false, indent::Int=0)
    for _ in 1:indent; print(io, ' '); end
    print(io, "level="); print(io, level(node));
    print(io, " nchildren="); print(io, length(node));
    print(io, " mbr=("); print(io, mbr(node).low); print(io, ", "); print(io, mbr(node).high); print(io, ")")
    if recurse && !isempty(node)
        print(io, ':')
        if node isa Leaf
            for child in children(node)
                println(io)
                for _ in 1:(indent+1); print(io, ' '); end
                show(io, child)
            end
        else
            for child in children(node)
                println(io)
                _show(io, child, recurse=true, indent=indent+1)
            end
        end
    end
end
