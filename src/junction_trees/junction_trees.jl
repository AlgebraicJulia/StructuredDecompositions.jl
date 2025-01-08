"""
    JunctionTree <: AbstractVector{Bag}

A [junction tree](https://en.wikipedia.org/wiki/Tree_decomposition).
This type implements the [indexed tree interface](https://juliacollections.github.io/AbstractTrees.jl/stable/#The-Indexed-Tree-Interface).
"""
struct JunctionTree <: AbstractVector{Bag}
    tree::Tree
    sndptr::Vector{Int}
    sepptr::Vector{Int}
    sepval::Vector{Int}
    relval::Vector{Int}
end


function Bag(tree::JunctionTree, i::Integer)
    Bag(residual(tree, i), separator(tree, i))
end


"""
    junctiontree(matrix::AbstractMatrix;
        alg::PermutationOrAlgorithm=AMDJL_AMD(),
        snd::SupernodeType=Maximal())

A non-mutating version of [`junctiontree!`](@ref).
"""
function junctiontree(matrix::AbstractMatrix; alg::PermutationOrAlgorithm=DEFAULT_ELIMINATION_ALGORITHM, snd::SupernodeType=DEFAULT_SUPERNODE_TYPE)
    junctiontree!(sparse(matrix); alg, snd)
end


function junctiontree(matrix::SparseMatrixCSC; alg::PermutationOrAlgorithm=DEFAULT_ELIMINATION_ALGORITHM, snd::SupernodeType=DEFAULT_SUPERNODE_TYPE)
    label, index = permutation(matrix, alg)
    cache = triu(matrix)
    label, lower, tree, cache = junctiontree!(label, sympermute(cache, index), snd, cache)
    label, tree
end


"""
    junctiontree!(matrix::SparseMatrixCSC;
        alg::PermutationOrAlgorithm=AMDJL_AMD(),
        snd::SupernodeType=Maximal())

Construct a [tree decomposition](https://en.wikipedia.org/wiki/Tree_decomposition) of a [simple graph](https://mathworld.wolfram.com/SimpleGraph.html), represented by its adjacency matrix `matrix`.
The vertices of the graph are first ordered by a [fill-reducing permutation](https://www.mathworks.com/help/matlab/math/sparse-matrix-reordering.html) computed by the algorithm `alg`.
The size of the resulting decomposition is determined by the supernode partition `snd`.

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
    label, index = permutation(matrix, alg)
    cache = triu(matrix)
    label, lower, tree, cache = junctiontree!(label, sympermute!(matrix, cache, index), snd, cache)
    label, tree
end


# Construct a junction tree. 
function junctiontree!(label::AbstractVector, upper::SparseMatrixCSC, snd::SupernodeType, cache::SparseMatrixCSC=similar(upper))
    label, lower, tree, sndptr, sepptr, cache = supernodetree!(label, upper, snd, cache)
    sepval = sepvals(lower, tree, sndptr, sepptr)
    relval = relvals(tree, sndptr, sepptr, sepval)
    label, lower, JunctionTree(tree, sndptr, sepptr, sepval, relval), cache
end


# Construct a postordered elimination tree.
function eliminationtree!(label::AbstractVector, upper::SparseMatrixCSC, cache::SparseMatrixCSC=similar(upper))
    tree = etree(upper)
    index = dfs(tree)
    invpermute!(label, index), sympermute!(cache, upper, index), invpermute!(tree, index), upper
end


# Construct a postordered supernodal elimination tree.
function supernodetree!(label::AbstractVector, upper::SparseMatrixCSC, snd::SupernodeType, cache::SparseMatrixCSC=similar(upper))
    label, upper, etree, cache = eliminationtree!(label, upper, cache)
    lower = transpose!(cache, upper)
    rowcount, colcount = supcnt(lower, etree)
    new, ancestor, tree = stree(etree, colcount, snd)
    index = dfs(tree)

    sndptr = Vector{Int}(undef, length(tree) + 1)
    sepptr = Vector{Int}(undef, length(tree) + 1)
    eindex = Vector{Int}(undef, size(lower, 1))
    p = sndptr[1] = sepptr[1] = 1

    for (i, j) in enumerate(invperm(index))
        v = new[j]

        while !isnothing(v) && v != ancestor[j]
            eindex[v] = p
            v = parentindex(etree, v)
            p += 1
        end

        sndptr[i + 1] = p
        sepptr[i + 1] = sndptr[i] + sepptr[i] + colcount[new[j]] - p
    end

    invpermute!(label, eindex), sympermute!(upper, lower, eindex, ReverseOrdering()), invpermute!(tree, index), sndptr, sepptr, lower
end


# Construct a postordered supernodal elimination tree.
function supernodetree!(label::AbstractVector, upper::SparseMatrixCSC, snd::Nodal, cache::SparseMatrixCSC=similar(upper))
    label, upper, tree, cache = eliminationtree!(label, upper, cache)
    lower = transpose!(cache, upper)
    rowcount, colcount = supcnt(lower, tree)
    sndptr = OneTo(size(lower, 1) + 1)
    sepptr = cumsum(vcat(1, colcount .- 1))
    label, lower, tree, sndptr, sepptr, upper
end


# Get the separators of every node of a supernodal elimination tree.
function sepvals(lower::SparseMatrixCSC, tree::Tree, sndptr::AbstractVector, sepptr::AbstractVector)
    stack = sizehint!(Int[], maximum(i -> sepptr[i + 1] - sepptr[i], tree))
    sepval = sizehint!(Int[], last(sepptr) - 1)

    for j in tree
        notinsnd(v) = sndptr[j + 1] <= v
        column = @view rowvals(lower)[nzrange(lower, sndptr[j])]
        append!(sepval, Iterators.filter(notinsnd, column))

        for i in childindices(tree, j)
            self = @view sepval[sepptr[j]:end]
            child = @view sepval[sepptr[i]:sepptr[i + 1] - 1]
            union = mergesorted!(empty!(stack), self, Iterators.filter(notinsnd, child))
            append!(resize!(sepval, sepptr[j] - 1), union)
        end 
    end

    sepval
end


# Get the relative indices of every node of a supernodal elimination tree.
function relvals(tree::Tree, sndptr::AbstractVector, sepptr::AbstractVector, sepval::AbstractVector)
    relval = Vector{Int}(undef, length(sepval))
    p = 1

    for i in tree[1:end - 1]
        j = parentindex(tree, i)
        q = sepptr[j]

        while p < sepptr[i + 1] && sepval[p] < sndptr[j + 1]
            relval[p] = sepval[p] - sndptr[j] + 1
            p += 1
        end

        while p < sepptr[i + 1] && q < sepptr[j + 1]
            if sepval[p] <= sepval[q]
                relval[p] = q - sepptr[j] - sndptr[j] + sndptr[j + 1] + 1
                p += 1
            end

            q += 1
        end
    end

    relval
end


"""
    treewidth(tree::JunctionTree, i::Integer)

Compute the width of a junction tree.
"""
function treewidth(tree::JunctionTree)
    maximum(length, tree) - 1
end


"""
    residual(tree::JunctionTree, i::Integer)

Get the residual at node i.
"""
function residual(tree::JunctionTree, i::Integer)
    tree.sndptr[i]:tree.sndptr[i + 1] - 1
end


"""
    separator(tree::JunctionTree, i::Integer)

Get the separator at node i.
"""
function separator(tree::JunctionTree, i::Integer)
    view(tree.sepval, tree.sepptr[i]:tree.sepptr[i + 1] - 1)
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

julia> tree[parentindex(tree, 3)][relative(tree, 3)] == separator(tree, 3)
true
```
"""
function relative(tree::JunctionTree, i::Integer)
    view(tree.relval, tree.sepptr[i]:tree.sepptr[i + 1] - 1)
end


function Base.show(io::IO, ::MIME"text/plain", tree::JunctionTree)
    print(io, "$(length(tree))-element JunctionTree:\n")
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
    Bag(tree, i)
end


function Base.IndexStyle(::Type{JunctionTree})
    IndexLinear()
end


function Base.size(tree::JunctionTree)
    size(tree.tree)
end
