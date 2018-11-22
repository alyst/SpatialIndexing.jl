"""
Base class for implementating the pool of `T` objects.
The pool allows to reduce the stress on GC by collecting the unneeded objects
(`release!(pool, obj)`) and reusing them later (`acquire!(pool)`).
"""
abstract type AbstractPool{T} end

capacity(pool::AbstractPool) = pool.capacity
Base.eltype(::Type{AbstractPool{T}}) where T = T
Base.eltype(pool::AbstractPool) = eltype(typeof(pool))
Base.length(pool::AbstractPool) = length(pool.objs)

acquire!(pool::AbstractPool) = length(pool) > 0 ? pop!(pool.objs) : newelem(pool)

function release!(pool::AbstractPool{T}, obj::T) where T
    if length(pool) >= capacity(pool)
        #@warn "Pool capacity exceeded" FIXME add parameter to enable it selectively
    else
        push!(pool.objs, obj)
    end
end

"""
The default `AbstarctPool` implementation.
"""
struct Pool{T} <: AbstractPool{T}
    objs::Vector{T}
    capacity::Int
end
