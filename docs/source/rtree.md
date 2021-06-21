# [R-tree](@id rtree)

[R-tree](https://en.wikipedia.org/wiki/R-tree) organizes data into
hierarchical structure and ensures that:
  * minimal bounding rectangles (MBRs) of the nodes (rectangles that
    encompass all data elements in the subtree) stay compact,
  * MBRs of the nodes from the same R-tree level have minimal overlap
    with each other.

The key benefit of R-tree is its ability to rebalance itself
and maintain efficient structure while handling dynamic data (massive insertions
and deletions).

`SpatialIndexing` provides `RTree` type that supports:
  * different R-tree variants (classic [R-tree](https://en.wikipedia.org/wiki/R-tree),
    [R*-tree](https://en.wikipedia.org/wiki/R*_tree), linear and quadratic node splits)
  * `insert!(tree, item)`, `delete!(tree, item)` for element-wise insertion and deletion
  * bulk-loading of data using Overlap-minimizing Top-down (OMT) approach (`load!(tree, data)`)
  * `subtract!(tree, reg)` for removing data within specified region `reg`
  * `findfirst(tree, reg, [id])`, `contained_in(tree, reg)` and `intersects_with(tree, reg)` spatial queries

```@docs
RTree
SpatialIndexing.RTreeVariant
Base.insert!
Base.delete!
```
