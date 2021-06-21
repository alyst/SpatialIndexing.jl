"""
Base abstract class for implementing regions in `N`-dimensional space with
dimensions of type `T`.
"""
abstract type Region{T,N} end

dimtype(::Type{<:Region{T}}) where T = T
dimtype(r::Region) = dimtype(typeof(r))
Base.ndims(::Type{<:Region{<:Any, N}}) where N = N
Base.ndims(r::Region) = ndims(typeof(r))

"""
`N`-dimensional point.
"""
struct Point{T, N} <: Region{T, N}
    coord::NTuple{N, T}
end

"""
    empty(::Type{T}) where T<:Region

Generate empty (uninitialized) region of type `T`.
"""
empty(::Type{Point{T,N}}) where {T,N} =
    Point{T,N}(ntuple(_ -> convert(T, NaN), N))
empty(::Type{Point{T,N}}) where {T<:Integer,N} =
    Point{T,N}(ntuple(_ -> zero(T), N))

"""
    area(a::Region)

`N`-dimensional "area" (volume etc) of `a`.
"""
area(a::Point{T}) where T = zero(typeof(zero(T)*zero(T)))

"""
    perimeter(a::Region)

The sum of the `a` sides.
"""
perimeter(a::Point{T}) where T = zero(T)

"""
    isvalid(a::Region)

Check that the parameters of `a` are valid and it defines a proper region.
"""
isvalid(a::Point) = all(!isnan, a.coord)

"""
Rectangular region constrained by `low[i]`...`high[i]` in each of `N` dimensions.
"""
struct Rect{T, N} <: Region{T, N}
    low::NTuple{N, T}
    high::NTuple{N, T}
end

Rect(low::NTuple{N, T1}, high::NTuple{N, T2}) where {N, T1<:Number, T2<:Number} =
    Rect{Base.promote_type(T1, T2), N}(low, high)

Rect{T,N}(pt::Point{T,N}) where {T,N} = Rect{T,N}(pt.coord, pt.coord)
Rect(pt::Point) = Rect(pt.coord, pt.coord)

# "empty" (uninitialized) rectangles
empty(::Type{Rect{T,N}}) where {T,N} =
    Rect{T,N}(ntuple(_ -> convert(T, NaN), N), ntuple(_ -> convert(T, NaN), N))
empty(::Type{Rect{T,N}}) where {T<:Integer,N} =
    Rect{T,N}(ntuple(_ -> typemax(T), N), ntuple(_ -> typemax(T), N))

"""
    combine(a::Region, b::Region)

MBR that contains both `a` and `b` regions.
"""
combine(a::Rect{T,N}, b::Rect{T,N}) where {T,N} =
    @inbounds Rect{T,N}(ntuple(i -> min(a.low[i], b.low[i]), N),
                        ntuple(i -> max(a.high[i], b.high[i]), N))

# check that coords are not NaN and side lengths not negative.
@generated isvalid(a::Rect{T,N}) where {T,N} =
    quote
        Base.Cartesian.@nall $N i -> !isnan(a.low[i]) && !isnan(a.high[i]) && a.low[i] <= a.high[i]
    end

"""
    intersect(a::Region, b::Region)

MBR that is an intersection of `a` and `b` regions.
"""
function intersect(a::Rect{T,N}, b::Rect{T,N}) where {T,N}
    res = @inbounds Rect{T,N}(ntuple(i -> max(a.low[i], b.low[i]), N),
                              ntuple(i -> min(a.high[i], b.high[i]), N))
    return isvalid(res) ? res : empty(typeof(a))
end

@generated area(a::Rect{T,N}) where {T,N} =
    quote
        @inbounds res = a.high[1] - a.low[1]
        @inbounds Base.Cartesian.@nexprs $(N-1) i -> res *= (a.high[i+1] - a.low[i+1])
        return res
    end

@generated perimeter(a::Rect{T,N}) where {T,N} =
    quote
        @inbounds res = a.high[1] - a.low[1]
        @inbounds Base.Cartesian.@nexprs $(N-1) i -> res += (a.high[i+1] - a.low[i+1])
        return res*ndims(a)
    end

"""
    overlap_area(a::Region, b::Region)

The area of MBR for `a` and `b` intersection.
"""
@generated overlap_area(a::Rect{T,N}, b::Rect{T,N}) where {T,N} =
    quote
        @inbounds res = min(a.high[1], b.high[1]) - max(a.low[1], b.low[1])
        res <= 0.0 && return zero(typeof(res))
        @inbounds Base.Cartesian.@nexprs $(N-1) i -> begin
            width_i = min(a.high[i+1], b.high[i+1]) - max(a.low[i+1], b.low[i+1])
            width_i <= zero(typeof(width_i)) && return zero(typeof(res))
            res *= width_i
        end
        return res
    end

"""
    combined_area(a::Region, b::Region)

The area of MBR for the union of `a` and `b`.
"""
@generated combined_area(a::Rect{T,N}, b::Rect{T,N}) where {T,N} =
    quote
        @inbounds res = max(a.high[1], b.high[1]) - min(a.low[1], b.low[1])
        @inbounds Base.Cartesian.@nexprs $(N-1) i -> res *= max(a.high[i+1], b.high[i+1]) - min(a.low[i+1], b.low[i+1])
        return res
    end

"""
    enlargement(a::Region, b::Region)

How much `a` grows when combined with `b`.

The difference between the MBR area of the `a` and `b` union and the MBR area of `a`.
"""
enlargement(a::Rect{T,N}, b::Rect{T,N}) where {T,N} = combined_area(a, b) - area(a)

"""
    center(a::Region)

The center point of `a`.
"""
center(a::Rect{T,N}) where {T,N} =
    Point(ntuple(i -> 0.5*(a.low[i] + a.high[i]), N))

@generated sqrdistance(a::Point{T,N}, b::Point{T,N}) where {T,N} =
    quote
        @inbounds res = abs2(a.coord[1] - b.coord[1])
        @inbounds Base.Cartesian.@nexprs $(N-1) i -> res += abs2(a.coord[i+1] - b.coord[i+1])
        return res
    end

"""
    intersects(a::Region, b::Region)

Check whether `a` intersects with `b`.
"""
@generated intersects(a::Rect{T,N}, b::Rect{T,N}) where {T,N} =
    quote
        @inbounds Base.Cartesian.@nall $N i -> (a.low[i] <= b.low[i] <= a.high[i]) || (b.low[i] <= a.low[i] <= b.high[i])
    end

"""
    contains(a::Region, b::Region)

Check whether `a` contains `b`.
"""
@generated contains(a::Rect{T,N}, b::Rect{T,N}) where {T,N} =
    quote
        @inbounds Base.Cartesian.@nall $N i -> (a.low[i] <= b.low[i]) && (a.high[i] >= b.high[i])
    end

@generated contains(a::Rect{T,N}, b::Point{T,N}) where {T,N} =
    quote
        @inbounds Base.Cartesian.@nall $N i -> (a.low[i] <= b.coord[i] <= a.high[i])
    end

"""
    in(a::Region, b::Region)

Check whether `a` is contained inside `b`.
"""
Base.in(a::Rect, b::Rect) = contains(b, a)
Base.in(a::Point, b::Rect) = contains(b, a)

# point is equal to the MBR if its low and high corners coincides
Base.:(==)(a::Point, b::Rect) = a.coord == b.low == b.high
Base.:(==)(a::Rect, b::Point) = b == a

"""
    touches(a::Rect, b::Rect)

Check whether `a` and `b` rectangles touches from the inside
(i.e. any `low` side touches `low` or `high` touches `high`).
"""
@generated touches(a::Rect{T,N}, b::Rect{T,N}) where {T,N} =
    quote
        @inbounds Base.Cartesian.@nany $N i -> (a.low[i] == b.low[i]) || (a.high[i] == b.high[i])
    end
