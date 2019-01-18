#=
Arrays

Examples of how to define high-performance custom arrays. Arrays are implemented as a concrete type.

Required method definitions
 - `size`: size of the array as a tuple
 - `getindex`: access array element(s)

Optional method definitions:
 - `IndexStyle`: either IndexCartesian (default) or IndexLinear if linearly indexible
=#

"""
    KmerVector(sequence::AbstractSequence, k::Int)

Construct a vector of length `k` k-mers in a sequence `sequence`.
"""
struct KmerVector{T<:AbstractSequence} <: AbstractVector{T}
    sequence::T
    k::Int
    function KmerVector(sequence::T, k::Int) where T<:AbstractSequence
        k > 1 ? new{T}(sequence, k) : throw(DomainError(k, "`k` must be >= 2"))
    end
end

sequence(x::KmerVector) = x.sequence
kay(x::KmerVector) = x.k

Base.size(x::KmerVector) = (length(sequence(x)) - kay(x) + 1,)
Base.getindex(x::KmerVector, i::Integer) = getindex(sequence(x), i:i+kay(x)-1)
Base.getindex(x::KmerVector, inds) = map(i->getindex(x, i), inds)
Base.IndexStyle(::KmerVector) = IndexLinear()

#=
Locality-sensitive hashing.
=#

abstract type AbstractSketch{T} <: AbstractVector{T} end

"""
    hashes(x::AbstractSketch)

Returns the hashes that make up the sketch `x`.
"""
hashes(x::AbstractSketch) = x.hashes


"""
    kay(x::AbstractSketch)

Returns the paramater `k` that defines the sketch `x`.
"""
kay(x::AbstractSketch) = length(hashes(x))

Base.size(x::AbstractSketch) = size(hashes(x))
Base.getindex(x::AbstractSketch, inds...) = getindex(hashes(x), inds...)
Base.IndexStyle(::AbstractSketch) = IndexLinear()

"""
    jaccard(a, b)

Returns the Jaccard similarity of `a` and `b` in the range [0,1].
"""
function jaccard(a::T, b::T) where T<:AbstractArray
    length(a) == 0 && length(b) == 0 && return 1.
    a = Set(a)
    b = Set(b)
    length(a ∩ b) / length(a ∪ b)
end

"""
    MinHashSketch(collection, k)

Construct a `k` k-wise MinHash sketch of `collection`.

Elements are hashed using the 64-bit MurmurHash algorithm.
"""
struct MinHashSketch{T<:Unsigned} <: AbstractSketch{T}
    hashes::Vector{T}
end

function MinHashSketch(xs, k::Int)
    hashes = zeros(UInt, k)

    # Hash `xs` with k separate hashing functions
    for i = 1:k
        seed = UInt(i)
        hmin = typemax(UInt)

        # Find MinHash for each hashing function
        for j = 1:length(xs)
            h = hash(xs[j], seed)

            if h < hmin
                hmin = h
            end
        end

        hashes[i] = hmin
    end

    MinHashSketch{eltype(hashes)}(hashes)
end

"""
    BottomKSketch(collection, k[, seed])

Construct a `k` bottom k sketch of `collection`.

Elements are hashed using the 64-bit MurmurHash algorithm.
"""
struct BottomKSketch{T<:Unsigned} <: AbstractSketch{T}
    hashes::Vector{T}
    seed::T
end

function BottomKSketch(xs, k::Int, seed::Integer=0)
    n = length(xs)
    k ≤ n || throw(ArgumentError("`k` > `length(xs)`"))
    seed = UInt(seed)
    heap = BinaryMaxHeap([typemax(UInt)])

    # Hash all elements in `xs` with a single hashing function
    for i = 1:n
        h = hash(xs[i], seed)

        # Keep track of bottom `k` hashes in a heap
        if h < top(heap)
            # Length of heap is always ≤ k
            if length(heap) == k
                pop!(heap)
            end

            push!(heap, h)
        end
    end

    hashes = heap.valtree
    BottomKSketch{eltype(hashes)}(hashes, seed)
end

seed(x::BottomKSketch) = x.seed
