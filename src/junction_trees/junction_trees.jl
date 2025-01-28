"""
    JunctionTree{I} <: AbstractVector{Bag{I}}

A [junction tree](https://en.wikipedia.org/wiki/Tree_decomposition).
This type implements the [indexed tree interface](https://juliacollections.github.io/AbstractTrees.jl/stable/#The-Indexed-Tree-Interface).
"""
struct JunctionTree{I} <: AbstractVector{Bag{I}}
    tree::SupernodeTree{I}
    sepptr::Vector{I}
    sepval::Vector{I}
    relval::Vector{I}

    function JunctionTree{I}(tree::SupernodeTree, sepptr::AbstractVector, sepval::AbstractVector, relval::AbstractVector) where I
        # validate arguments
        length(tree) != length(sepptr) - 1 && throw(ArgumentError("length(tree) != length(sepptr) - 1"))
        sepptr[1] != 1 && throw(ArgumentError("sepptr[1] != 1"))
        sepptr[end] != length(sepval) + 1 && throw(ArgumentError("sepptr[end] != length(sepval) + 1"))
        sepptr[end] != length(relval) + 1 && throw(ArgumentError("sepptr[end] != length(relval) + 1"))

        # construct tree
        new{I}(tree, sepptr, sepval, relval)
    end
end


function JunctionTree(tree::SupernodeTree{I}, sepptr::AbstractVector{I}, sepval::AbstractVector{I}, relval::AbstractVector{I}) where I
    JunctionTree{I}(tree, sepptr, sepval, relval)
end


function Tree(tree::JunctionTree)
    Tree(tree.tree)
end


function Tree{I}(tree::JunctionTree) where I
    Tree{I}(tree.tree)
end


"""
    junctiontree(graph;
        alg::PermutationOrAlgorithm=DEFAULT_ELIMINATION_ALGORITHM,
        snd::SupernodeType=DEFAULT_SUPERNODE_TYPE)

Construct a [tree decomposition](https://en.wikipedia.org/wiki/Tree_decomposition) of a simple graph.
The vertices of the graph are first ordered by a fill-reducing permutation computed by the algorithm `alg`.
The size of the resulting decomposition is determined by the supernode partition `snd`.
```julia
julia> using SparseArrays, StructuredDecompositions

julia> graph = sparse([
           0 1 1 0 0 0 0 0
           1 0 1 0 0 1 0 0
           1 1 0 1 1 0 0 0
           0 0 1 0 1 0 0 0
           0 0 1 1 0 0 1 1
           0 1 0 0 0 0 1 0
           0 0 0 0 1 1 0 1
           0 0 0 0 1 0 1 0
       ]);

julia> label, tree = junctiontree(graph);

julia> tree
6-element JunctionTree{Int64}:
[6, 7, 8]
├─ [1, 6, 7]
├─ [4, 6, 8]
│  └─ [3, 4, 6]
│     └─ [2, 3, 6]
└─ [5, 7, 8]
```
"""
function junctiontree(graph; alg::PermutationOrAlgorithm=DEFAULT_ELIMINATION_ALGORITHM, snd::SupernodeType=DEFAULT_SUPERNODE_TYPE)
    junctiontree(graph, alg, snd)
end


function junctiontree(graph, alg::PermutationOrAlgorithm, snd::SupernodeType)
    junctiontree!(supernodetree(graph, alg, snd)...)
end


"""
    junctiontree!(graph;
        alg::PermutationOrAlgorithm=DEFAULT_ELIMINATION_ALGORITHM,
        snd::SupernodeType=DEFAULT_SUPERNODE_TYPE)

A mutating version of [`junctiontree`](@ref).
"""
function junctiontree!(graph; alg::PermutationOrAlgorithm=DEFAULT_ELIMINATION_ALGORITHM, snd::SupernodeType=DEFAULT_SUPERNODE_TYPE)
    label, tree, lower, cache = junctiontree!(graph, alg, snd)
    label, tree
end


function junctiontree!(graph, alg::PermutationOrAlgorithm, snd::SupernodeType)
    junctiontree!(supernodetree!(graph, alg, snd)...)
end


function junctiontree!(label::Vector{I}, tree::SupernodeTree{I}, index::Vector{I}, sepptr::Vector{I}, lower::SparseMatrixCSC{Nothing, I}, cache::SparseMatrixCSC{Nothing, I}) where I
    lower = sympermute!(cache, lower, index, ReverseOrdering())
    sepval = sepvals(lower, tree, sepptr)
    relval = relvals(tree, sepptr, sepval)
    label, JunctionTree(tree, sepptr, sepval, relval), cache, lower
end


"""
    treewidth(tree::JunctionTree)

Compute the width of a junction tree.
"""
function treewidth(tree::JunctionTree{I}) where I
    width::I = maximum(length, tree; init=1) - 1
    width
end


"""
    treewidth(graph;
        alg::PermutationOrAlgorithm=DEFAULT_ELIMINATION_ALGORITHM)

Compute an upper bound to the [tree width](https://en.wikipedia.org/wiki/Treewidth) of a simple graph.
See [`junctiontree`](@ref) for the meaning of `alg`.
"""
function treewidth(graph; alg::PermutationOrAlgorithm=DEFAULT_ELIMINATION_ALGORITHM)
    treewidth(graph, alg)
end


function treewidth(graph, alg::PermutationOrAlgorithm)
    label, tree, upper = eliminationtree(graph, alg)
    cache = spzeros(Nothing, indtype(upper), size(upper))
    treewidth!(tree, upper, cache)
end


"""
    treewidth!(graph;
        alg::PermutationOrAlgorithm=DEFAULT_ELIMINATION_ALGORITHM)

A mutating version of [`treewidth`](@ref).
"""
function treewidth!(graph; alg::PermutationOrAlgorithm=DEFAULT_ELIMINATION_ALGORITHM)
    treewidth!(graph, alg)
end


function treewidth!(matrix::SparseMatrixCSC{<:Any, I}, alg::PermutationOrAlgorithm) where I
    label, tree, upper = eliminationtree(matrix, alg)
    cache = SparseMatrixCSC{Nothing, I}(size(matrix)..., getcolptr(matrix), rowvals(matrix), Vector{Nothing}(undef, nnz(matrix)))
    treewidth!(tree, upper, cache)
end


function treewidth!(tree::Tree{I}, upper::SparseMatrixCSC{Nothing, I}, cache::SparseMatrixCSC{Nothing, I}) where I
    lower = transpose!(cache, upper)
    rowcount, colcount = supcnt(lower, tree)
    maximum(colcount; init=one(I)) - 1
end


"""
    residual(tree::JunctionTree, i::Integer)

Get the residual at node `i`.
"""
function residual(tree::JunctionTree, i::Integer)
    tree.tree[i]
end


function seprange(tree::JunctionTree{I}, i::Integer) where I
    tree.sepptr[i]:tree.sepptr[i + 1] - one(I)
end


"""
    separator(tree::JunctionTree, i::Integer)

Get the separator at node `i`.
"""
function separator(tree::JunctionTree, i::Integer)
    @view tree.sepval[seprange(tree, i)]
end


"""
    relative(tree::JunctionTree, i::Integer)

Get the relative indices at node `i`.
"""
function relative(tree::JunctionTree, i::Integer)
    @view tree.relval[seprange(tree, i)]
end


##########################
# Indexed Tree Interface #
##########################


function AbstractTrees.rootindex(tree::JunctionTree)
    rootindex(tree.tree)
end


function AbstractTrees.parentindex(tree::JunctionTree, i::Integer)
    parentindex(tree.tree, i)
end


function firstchildindex(tree::JunctionTree, i::Integer)
    firstchildindex(tree.tree, i)
end


function AbstractTrees.nextsiblingindex(tree::JunctionTree, i::Integer)
    nextsiblingindex(tree.tree, i)
end


function rootindices(tree::JunctionTree)
    rootindices(tree.tree)
end


function AbstractTrees.childindices(tree::JunctionTree, i::Integer)
    childindices(tree.tree, i)
end


function ancestorindices(tree::JunctionTree, i::Integer)
    ancestorindices(tree.tree, i)
end


#############################
# Abstract Vector Interface #
#############################


function Base.getindex(tree::JunctionTree, i::Integer)
    Bag(residual(tree, i), separator(tree, i))
end


function Base.IndexStyle(::Type{JunctionTree})
    IndexLinear()
end


function Base.size(tree::JunctionTree)
    size(tree.tree)
end
