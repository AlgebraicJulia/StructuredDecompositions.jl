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

    function JunctionTree{I}(tree::SupernodeTree, sepptr::AbstractVector, sepval::AbstractVector) where I
        # validate arguments
        length(tree) != length(sepptr) - 1 && throw(ArgumentError("length(tree) != length(sepptr) - 1"))
        sepptr[1] != 1 && throw(ArgumentError("sepptr[1] != 1"))
        sepptr[end] != length(sepval) + 1 && throw(ArgumentError("sepptr[end] != length(sepval) + 1"))

        # construct tree
        tree = new{I}(tree, sepptr, sepval, I[])
    end
end


function JunctionTree(tree::SupernodeTree{I}, sepptr::AbstractVector{I}, sepval::AbstractVector{I}) where I
    JunctionTree{I}(tree, sepptr, sepval)
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
    label, tree, index, sepptr, lower, upper = supernodetree(graph, alg, snd)
    lower = sympermute!(upper, lower, index, ReverseOrdering())
    label, JunctionTree(tree, sepptr, sepvals(lower, tree, sepptr))
end


"""
    treewidth(tree::JunctionTree)

Compute the width of a junction tree.
"""
function treewidth(tree::JunctionTree)
    maximum(length, tree; init=1) - 1
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
    rowcount, colcount = supcnt(copy(transpose(upper)), tree)
    maximum(colcount; init=1) - 1
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
These are not cached by default; compute them by running `relative!(tree)`.
"""
function relative(tree::JunctionTree, i::Integer)
    @view tree.relval[seprange(tree, i)]
end


function relative!(tree::JunctionTree)
    resize!(tree.relval, length(tree.sepval))

    for (j, bag) in enumerate(tree)
        for i in childindices(tree, j)
            indexinsorted!(relative(tree, i), separator(tree, i), bag)
        end
    end

    tree
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
