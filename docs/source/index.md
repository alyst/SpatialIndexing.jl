# SpatialIndexing.jl package

`SpatialIndexing` package provides the tools for efficient in-memory storage and retrieval of
spatial data in Julia (http://julialang.org/).

## Installation
```julia
using Pkg; Pkg.add("SpatialIndexing")
```
from Julia REPL.

## Spatial Indices

  * [spatial primitives](@ref regions)
  * [basic types](@ref abstract)
  * [R-tree, R*-tree](@ref rtree)
  * [simple index](@ref simple_index)
  * [spatial queries](@ref query)

## See also

Other Julia packages for spatial data:

  * [LibSpatialIndex.jl](https://github.com/yeesian/LibSpatialIndex.jl)
    ([libspatialindex](https://github.com/libspatialindex/libspatialindex) wrapper)
  * [NearestNeighbors.jl](https://github.com/KristofferC/NearestNeighbors.jl)
  * [RegionTrees.jl](https://github.com/rdeits/RegionTrees.jl)
  * [LSH.jl](https://github.com/Keno/LSH.jl)

## References

* A.Guttman, _“R-Trees: A Dynamic Index Structure for Spatial Searching”_
  Proc. 1984 ACM-SIGMOD Conference on Management of Data (1985), pp.47-57.
* N. Beckmann, H.P. Kriegel, R. Schneider, B. Seeger,
  _"The R*-tree: an efficient and robust access method for points and rectangles"_
  Proc. 1990 ACM SIGMOD international conference on Management of data (1990), p.322
* T. Lee and S. Lee, _"OMT: Overlap Minimizing Top-down Bulk Loading Algorithm for R-tree"_,
  CAiSE Short Paper Proceedings (2003) [paper](http://ceur-ws.org/Vol-74/files/FORUM_18.pdf)
