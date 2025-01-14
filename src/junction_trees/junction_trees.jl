"""
    JunctionTree <: AbstractVector{Bag}

A [junction tree](https://en.wikipedia.org/wiki/Tree_decomposition).
This type implements the [indexed tree interface](https://juliacollections.github.io/AbstractTrees.jl/stable/#The-Indexed-Tree-Interface).
"""
struct JunctionTree <: AbstractVector{Bag}
    tree::SupernodeTree
    sepptr::Vector{Int}
    sepval::Vector{Int}
    relval::Vector{Int}

    function JunctionTree(tree::SupernodeTree, sepptr::AbstractVector, sepval::AbstractVector, relval::AbstractVector)
        # validate arguments
        length(tree) != length(sepptr) - 1 && throw(ArgumentError("length(tree) != length(sepptr) - 1"))
        sepptr[1] != 1 && throw(ArgumentError("sepptr[1] != 1"))
        sepptr[end] != length(sepval) + 1 && throw(ArgumentError("sepptr[end] != length(sepval) + 1"))
        sepptr[end] != length(relval) + 1 && throw(ArgumentError("sepptr[end] != length(relval) + 1"))

        # construct tree
        new(tree, sepptr, sepval, relval)
    end
end


function JunctionTree(tree::JunctionTree)
    JunctionTree(tree.tree, tree.sepptr, tree.sepval, tree.relval)
end


function SupernodeTree(tree::JunctionTree)
    SupernodeTree(tree.tree)
end


function Tree(tree::JunctionTree)
    Tree(tree.tree)
end


function Tree(tree::JunctionTree, root::Integer)
    Tree(tree.tree, root)
end


"""
    junctiontree(matrix::AbstractMatrix;
        alg::PermutationOrAlgorithm=DEFAULT_ELIMINATION_ALGORITHM,
        snd::SupernodeType=DEFAULT_SUPERNODE_TYPE)

A non-mutating version of [`junctiontree!`](@ref).
"""
function junctiontree(matrix::AbstractMatrix; alg::PermutationOrAlgorithm=DEFAULT_ELIMINATION_ALGORITHM, snd::SupernodeType=DEFAULT_SUPERNODE_TYPE)
    junctiontree!(sparse(matrix); alg, snd)
end


"""
    junctiontree!(matrix::SparseMatrixCSC;
        alg::PermutationOrAlgorithm=DEFAULT_ELIMINATION_ALGORITHM,
        snd::SupernodeType=DEFAULT_SUPERNODE_TYPE)

Construct a [tree decomposition](https://en.wikipedia.org/wiki/Tree_decomposition) of a simple graph.
The vertices of the graph are first ordered by a fill-reducing permutation computed by the algorithm `alg`.
The size of the resulting decomposition is determined by the supernode partition `snd`.
- `matrix`: adjacency matrix
- `alg`: ordering algorithm
- `snd`: supernode partition type
```julia
julia> using StructuredDecompositions.JunctionTrees

julia> graph = [
           0 1 1 0 0 0 0 0
           1 0 1 0 0 1 0 0
           1 1 0 1 1 0 0 0
           0 0 1 0 1 0 0 0
           0 0 1 1 0 0 1 1
           0 1 0 0 0 0 1 0
           0 0 0 0 1 1 0 1
           0 0 0 0 1 0 1 0
       ];

julia> label, tree = junctiontree(graph);

julia> tree
6-element JunctionTree:
[6, 7, 8]
├─ [1, 6, 7]
├─ [4, 6, 8]
│  └─ [3, 4, 6]
│     └─ [2, 3, 6]
└─ [5, 7, 8]
```
"""
function junctiontree!(matrix::SparseMatrixCSC; alg::PermutationOrAlgorithm=DEFAULT_ELIMINATION_ALGORITHM, snd::SupernodeType=DEFAULT_SUPERNODE_TYPE)
    label, tree, lower, cache = junctiontree!(matrix, alg, snd)
    label, tree
end


# Construct a junction tree. 
function junctiontree!(matrix::SparseMatrixCSC, alg::PermutationOrAlgorithm, snd::SupernodeType)
    label, tree, index, sepptr, lower, cache = supernodetree!(matrix, alg, snd)
    sepval = sepvals(sympermute!(cache, lower, index, ReverseOrdering()), tree, sepptr)
    relval = relvals(tree, sepptr, sepval)
    label, JunctionTree(tree, sepptr, sepval, relval), cache, lower
end


"""
    treewidth(tree::JunctionTree)

Compute the width of a junction tree.
"""
function treewidth(tree::JunctionTree)
    maximum(length, tree) - 1
end


"""
    treewidth(matrix::AbstractMatrix;
        alg::PermutationOrAlgorithm=DEFAULT_ELIMINATION_ALGORITHM)

A non-mutating version of [`treewidth!`](@ref).
"""
function treewidth(matrix::AbstractMatrix; alg::PermutationOrAlgorithm=DEFAULT_ELIMINATION_ALGORITHM)
    treewidth!(sparse(matrix); alg)
end


"""
    treewidth!(matrix::SparseMatrixCSC;
        alg::PermutationOrAlgorithm=DEFAULT_ELIMINATION_ALGORITHM)

Compute an upper bound to the [tree width](https://en.wikipedia.org/wiki/Treewidth) of a simple graph.
See [`junctiontree!`](@ref) for the meaning of `alg`.
"""
function treewidth!(matrix::SparseMatrixCSC; alg::PermutationOrAlgorithm=DEFAULT_ELIMINATION_ALGORITHM)
    label, tree, upper, cache = eliminationtree!(matrix, alg)
    rowcount, colcount = supcnt(transpose!(cache, upper), tree)
    maximum(colcount) - 1
end


"""
    residual(tree::JunctionTree, i::Integer)

Get the residual at node `i`.
"""
function residual(tree::JunctionTree, i::Integer)
    tree.tree[i]
end


"""
    separator(tree::JunctionTree, i::Integer)

Get the separator at node `i`.
"""
function separator(tree::JunctionTree, i::Integer)
    @view tree.sepval[tree.sepptr[i]:tree.sepptr[i + 1] - 1]
end


"""
    relative(tree::JunctionTree, i::Integer)

Get the indices in `tree[parentindex(tree, i)]` corresponding to the elements of `separator(tree, i)`.

```julia
julia> using AbstractTrees

julia> using StructuredDecompositions.JunctionTrees

julia> graph = [
           0 1 1 0 0 0 0 0
           1 0 1 0 0 1 0 0
           1 1 0 1 1 0 0 0
           0 0 1 0 1 0 0 0
           0 0 1 1 0 0 1 1
           0 1 0 0 0 0 1 0
           0 0 0 0 1 1 0 1
           0 0 0 0 1 0 1 0
       ];

julia> label, tree = junctiontree(graph);

julia> bag = tree[parentindex(tree, 1)]
3-element Bag:
 6
 7
 8

julia> sep = separator(tree, 1)
2-element view(::Vector{Int64}, 1:2) with eltype Int64:
 6
 7

julia> rel = relative(tree, 1)
2-element view(::Vector{Int64}, 1:2) with eltype Int64:
 1
 2

julia> bag[rel] == sep
true
```
"""
function relative(tree::JunctionTree, i::Integer)
    @view tree.relval[tree.sepptr[i]:tree.sepptr[i + 1] - 1]
end


function Base.show(io::IO, ::MIME"text/plain", tree::JunctionTree)
    println(io, "$(length(tree))-element JunctionTree:")
    print_tree(io, IndexNode(tree))
end


###########################
# Abstract Tree Interface #
###########################


function firstchildindex(tree::JunctionTree, i::Integer)
    firstchildindex(tree.tree, i)
end


function AbstractTrees.rootindex(tree::JunctionTree)
    rootindex(tree.tree)
end


function AbstractTrees.parentindex(tree::JunctionTree, i::Integer)
    parentindex(tree.tree, i)
end


function AbstractTrees.nextsiblingindex(tree::JunctionTree, i::Integer)
    nextsiblingindex(tree.tree, i)
end


function AbstractTrees.childindices(tree::JunctionTree, i::Integer)
    childindices(tree.tree, i)
end


function AbstractTrees.ParentLinks(::Type{IndexNode{JunctionTree, Int}})
    StoredParents()
end


function AbstractTrees.SiblingLinks(::Type{IndexNode{JunctionTree, Int}})
    StoredSiblings()
end


function AbstractTrees.NodeType(::Type{IndexNode{JunctionTree, Int}})
    HasNodeType()
end


function AbstractTrees.nodetype(::Type{IndexNode{JunctionTree, Int}})
    IndexNode{JunctionTree, Int}
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
