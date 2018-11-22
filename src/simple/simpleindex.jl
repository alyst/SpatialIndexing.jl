# TODO
"""
Vector-based `SpatialIndex`.
While insertion is `O(1)`, the search is `O(N)`.

Generally should not be used, except for performance comparisons
or when the number of stored elements is expected to be very small (<100).
"""
struct SimpleSpatialIndex{T,N,V} <: SpatialIndex{T,N,V}
    elems::Vector{V}
end

Base.length(si::SimpleSpatialIndex) = length(si.elems)
Base.isempty(si::SimpleSpatialIndex) = isempty(si.elems)
