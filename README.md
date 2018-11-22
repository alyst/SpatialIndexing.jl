SpatialIndexing.jl
==============

`SpatialIndexing` package provides the tools for efficient in-memory indexing of
spatial data in Julia (http://julialang.org/).

# Installation
```julia
using Pkg; Pkg.add("https://github.com/alyst/SpatialIndexing.jl")
```
from a Julia REPL.

# Features

## R-tree

[R-tree](https://en.wikipedia.org/wiki/R-tree) organizes the data into
hierarchical structure and ensures that the minimal bounding rectangles (MBRs)
of the nodes (rectangles that encompass all attached data elements) stay compact,
and that the MBRs of the nodes from the same R-tree level have minimal overlap
with each other. The key benefit of R-tree is its ability to rebalance itself
and maintain efficient structure while handling dynamic data (massive insertions
and deletions).

`RTree` implementation supports:
  * different R-tree variants (classic `R`-tree, [R*-tree](https://en.wikipedia.org/wiki/R*_tree),
linear and quadratic node splits)
* `insert!()`, `delete!()` for element-wise insertion and deletion
  * bulk-loading of data using Overlap-minimizing Top-down (OMT) approach (`load!()`)
  * `subtract!()` for removing data within specified region
  * TODO: spatial queries

## Simple Spatial Index

`SimpleSpatialIndex` stores all data elements in a vector. So, while insertion
of new data takes constant time, the time of spatial searches grows linearly
with the number of elements. This spatial index is intended as a reference
implementation for benchmarking and not recommended for production usage.

TODO

# Usage

`examples` folder contains `spiral.jl` and `pareto.jl` examples of using R-tree
for storing spatial data.

![R*-tree of 10000 random points (sequential insertions)](examples/spiral_rtree_seq.png)

# See also

Other Julia packages implementing spatial indexing:

  * [LibSpatialIndexing.jl](https://github.com/yeesian/LibSpatialIndex.jl)
  * [NearestNeighbors.jl](https://github.com/KristofferC/NearestNeighbors.jl)
  * [RegionTrees.jl](https://github.com/rdeits/RegionTrees.jl)
  * [LSH.jl](https://github.com/Keno/LSH.jl)
