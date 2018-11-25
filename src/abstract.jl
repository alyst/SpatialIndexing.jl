"""
`SpatialIndex`-related exception raised within Julia
"""
struct SpatialIndexException <: Exception
    msg::AbstractString     # error message
end

function Base.showerror(io::IO, err::SpatialIndexException)
    print(io, "SpatialIndexException: ")
    print(io, err.msg)
end

# throws KeyError with or without id (depending on eltype HasID trait)
__spatial_keyerror(eltype::Type, br::Region, id::Any = nothing) =
    idtrait(eltype) === HasNoID ? throw(KeyError(br)) : throw(KeyError((br, id)))

@enum QueryType::Int QueryContainedIn QueryIntersectsWith QueryPoint QueryNearestNeighbours

"""
Base abstract class for spatial indexing of elements of type `V`
in `N`-dimensional space with dimensions of type `T`.
"""
abstract type SpatialIndex{T<:Number, N, V} end

"""
Base abstract class for implementing spatial queries in `N`-dimensional space.
"""
abstract type SpatialQueryIterator{T<:Number,N,Q} end

# basic SpatialIndex API
Base.eltype(::Type{<:SpatialIndex{<:Any, <:Any, V}}) where V = V
Base.eltype(si::SpatialIndex) = eltype(typeof(si))
idtype(::Type{T}) where T<:SpatialIndex = idtype(idtrait(eltype(T)))
idtype(si::SpatialIndex) = idtype(typeof(si))
dimtype(::Type{<:SpatialIndex{T}}) where T = T
dimtype(si::SpatialIndex) = dimtype(typeof(si))
Base.ndims(::Type{<:SpatialIndex{<:Any, N}}) where N = N
Base.ndims(si::SpatialIndex) = ndims(typeof(si))
regiontype(::Type{<:SpatialIndex{T,N}}) where {T,N} = Region{T,N} # default unless overridden
regiontype(si::SpatialIndex) = regiontype(typeof(si))
Base.length(si::SpatialIndex) = si.nelems
Base.isempty(si::SpatialIndex) = length(si) == 0

# arbitrary spatial elements support

# Type trait for supporting id() method
abstract type HasID{T} end
abstract type HasNoID end
idtrait(::Type) = HasNoID
idtype(::Type{HasID{T}}) where T = T
idtype(::Type{HasNoID}) = Union{}

# Type trait for supporting mbr() method
abstract type HasMBR{T} end
abstract type HasNoMBR end
mbrtrait(::Type) = HasNoMBR
mbrtype(::Type{HasMBR{T}}) where T = T
mbrtype(::Type{HasNoMBR}) = Union{}

# check whether type V has HasMBR trait and that the returned object
# is `Rect{T,N}`
function check_hasmbr(::Type{Rect{T,N}}, ::Type{V}) where {T, N, V}
    R = mbrtrait(V)
    R <: HasMBR ||
        throw(ArgumentError("Element type $V doesn't have mbr() method"))
    mbrtype(R) === Rect{T,N} ||
        throw(ArgumentError("Element type $V: MBR type ($(mbrtype(R))) incompatible with ($(Rect{T,N}))"))
    return true
end

# check whether type V has HasID trait and that the returned object is `K`
function check_hasid(::Type{K}, ::Type{V}) where {K, V}
    R = idtrait(V)
    R <: HasID ||
        throw(ArgumentError("Element type $V doesn't have id() method"))
    idtype(R) === K ||
        throw(ArgumentError("Element type $V: ID type ($(idtype(R))) incompatible with ($K)"))
    return true
end

"""
Simple `N`-dimensional spatial data element that stores values of type `V`
and supports `HasMBR` and `HasID{K}` (if `K` is not `Nothing`) traits.
"""
struct SpatialElem{T,N,K,V}
    mbr::Rect{T,N}
    id::K
    val::V
end

idtrait(::Type{<:SpatialElem{<:Any,<:Any,K}}) where K = HasID{K}
idtrait(::Type{<:SpatialElem{<:Any,<:Any,Nothing}}) = HasNoID
id(el::SpatialElem) = el.id
mbrtrait(::Type{<:SpatialElem{T,N}}) where {T,N} = HasMBR{Rect{T,N}}
mbr(el::SpatialElem) = el.mbr
