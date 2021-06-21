# [Simple Spatial Index](@id simple_index)

`SimpleSpatialIndex` stores all data elements in a vector. So, while insertion
of new data takes constant time, the time of spatial searches grows linearly
with the number of elements. This spatial index is intended as a reference
implementation for benchmarking and not recommended for production usage.

```@docs
SimpleSpatialIndex
```
