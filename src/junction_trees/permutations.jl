"""
    Permutation <: AbstractVector{Int}

A [permutation](https://en.wikipedia.org/wiki/Permutation) ``\\sigma`` of the set ``\\{1, \\dots, n\\}``.
This type implements the [abstract array interface](https://docs.julialang.org/en/v1/manual/interfaces/#man-interface-array).
"""
struct Permutation <: AbstractVector{Int}
    inner::Vector{Int} # permutation
    index::Vector{Int} # inverse permutation
end


"""
    Permutation(permutation::AbstractVector)

Construct a permutation ``\\sigma`` from a sequence ``(\\sigma(1), \\dots, \\sigma(n))``.
"""
function Permutation(permutation::AbstractVector)
    Permutation(permutation, invperm(permutation))
end


"""
    Permutation(graph[, ealg::EliminationAlgorithm])

Construct a fill-reducing permutation of the vertices of a graph.
"""
function Permutation(graph)
    Permutation(graph, DEFAULT_ELIMINATION_ALGORITHM)
end


function Base.invperm(permutation::Permutation)
    Permutation(permutation.index, permutation.inner)
end


function Base.permute!(permutation::Permutation, other::AbstractVector)
    permute!(permutation.inner, other)
    permutation.index[permutation.inner] = eachindex(permutation.inner)
    permutation
end


function Base.copy(permutation::Permutation)
    Permutation(permutation.inner, permutation.index)
end


function Base.deepcopy(permutation::Permutation)
    Permutation(copy(permutation.inner), copy(permutation.index))
end


#############################
# Abstract Vector Interface #
#############################


function Base.getindex(permutation::Permutation, i)
    permutation.inner[i]
end


function Base.IndexStyle(::Type{Permutation})
    IndexLinear()
end


function Base.size(permutation::Permutation)
    (length(permutation.inner),)
end
