module SpatialIndexing

export SpatialIndex, SpatialIndexException, SpatialElem,
    RTree, SimpleSpatialIndex,
    contained_in, intersects_with

include("regions.jl")
include("pool.jl")

include("abstract.jl")

## Spatial Indices

#   SimpleSpatialIndex
include("simple/simpleindex.jl")

#   R-tree
include("rtree/node.jl")
include("rtree/rtree.jl")
include("rtree/check.jl")
include("rtree/show.jl")
include("rtree/query.jl")
include("rtree/split.jl")
include("rtree/insert.jl")
include("rtree/adjust.jl")
include("rtree/delete.jl")
include("rtree/bulk.jl")

end # module
