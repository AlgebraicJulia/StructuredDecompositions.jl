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
    junctiontree(matrix::SparseMatrixCSC; alg::PermutationOrAlgorithm=AMDJL_AMD(), snd::SupernodeType=Maximal())

A non-mutating version of [`junctiontree!`](@ref).
"""
function junctiontree(matrix::SparseMatrixCSC; alg::PermutationOrAlgorithm=DEFAULT_ELIMINATION_ALGORITHM, snd::SupernodeType=DEFAULT_SUPERNODE_TYPE)
    label, index = permutation(matrix, alg)
    cache = triu(matrix)
    label, lower, tree, cache = junctiontree!(label, sympermute(cache, index), snd, cache)
    label, tree
end


"""
    junctiontree!(matrix::SparseMatrixCSC; alg::PermutationOrAlgorithm=AMDJL_AMD(), snd::SupernodeType=Maximal())

Construct a [tree decomposition](https://en.wikipedia.org/wiki/Tree_decomposition) of a [simple graph](https://mathworld.wolfram.com/SimpleGraph.html) ``G = (V, E)``, represented by its adjacency matrix.
The vertices of ``G`` are first ordered by a [fill-reducing permutation](https://www.mathworks.com/help/matlab/math/sparse-matrix-reordering.html) ``\\sigma: V \\to V`` computed by the algorithm `alg`.
The size of the resulting junction tree ``J = (T, B)`` is determined by the supernode partition `snd`.
The function returns the pair ``(\\sigma, J)``.
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
    sndval = Vector{Int}(undef, size(lower, 1))
    sndptr[1] = sepptr[1] = 1

    for (i, j) in enumerate(invperm(index))
        v = new[j]
        p = sndptr[i]

        while !isnothing(v) && v != ancestor[j]
            sndval[v] = p
            v = parentindex(etree, v)
            p += 1
        end

        sndptr[i + 1] = p
        sepptr[i + 1] = sndptr[i] + sepptr[i] + colcount[new[j]] - p
    end

    invpermute!(label, sndval), sympermute!(upper, lower, sndval, ReverseOrdering()), invpermute!(tree, index), sndptr, sepptr, lower
end


# Construct a postordered supernodal elimination tree.
function supernodetree!(label::AbstractVector, upper::SparseMatrixCSC, snd::Node, cache::SparseMatrixCSC=similar(upper))
    label, upper, tree, cache = eliminationtree!(label, upper, cache)
    lower = transpose!(cache, upper)
    rowcount, colcount = supcnt(lower, tree)
    sndptr = OneTo(size(lower, 1) + 1)
    sepptr = cumsum(vcat(1, colcount .- 1))
    label, lower, tree, sndptr, sepptr, upper
end


# Get the separators of every node of a supernodal elimination tree.
function sepvals(lower::SparseMatrixCSC, tree::Tree, sndptr::AbstractVector, sepptr::AbstractVector)
    index = zeros(Int, size(lower, 1))
    sepval = Vector{Int}(undef, last(sepptr) - 1)
    p = 1

    for j in tree
        for v in view(rowvals(lower), nzrange(lower, sndptr[j]))
            if sndptr[j + 1] <= v
                sepval[p] = v
                index[v] = j
                p += 1
            end
        end

        for i in childindices(tree, j), v in view(sepval, sepptr[i]:sepptr[i + 1] - 1)
            if sndptr[j + 1] <= v && index[v] != j
                sepval[p] = v
                index[v] = j
                p += 1
            end
        end

        sort!(sepval, sepptr[j], sepptr[j + 1] - 1, DEFAULT_STABLE, ForwardOrdering())
    end

    sepval
end


# Get the relative indices of every node of a supernodal elimination tree.
function relvals(tree::Tree, sndptr::AbstractVector, sepptr::AbstractVector, sepval::AbstractVector)
    relval = Vector{Int}(undef, length(sepval))
    p = 1

    for i in tree[1:end - 1]
        j = parentindex(tree, i)

        while p < sepptr[i + 1] && sepval[p] < sndptr[j + 1]
            relval[p] = sepval[p] - sndptr[j] + 1
            p += 1
        end

        q = sepptr[j]

        while p < sepptr[i + 1]
            q = searchsortedfirst(sepval, sepval[p], q, sepptr[j + 1] - 1, ForwardOrdering()) + 1
            relval[p] = q - sepptr[j] - sndptr[j] + sndptr[j + 1]
            p += 1
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

Get the residual ``\\mathrm{res}(i) := \\mathrm{bag}(i) - \\mathrm{sep}(i).``
"""
function residual(tree::JunctionTree, i::Integer)
    tree.sndptr[i]:tree.sndptr[i + 1] - 1
end


"""
    separator(tree::JunctionTree, i::Integer)

Get the separator ``\\mathrm{sep}(i) := \\mathrm{bag}(i) \\cap \\mathrm{bag}(parent(i)).``
"""
function separator(tree::JunctionTree, i::Integer)
    view(tree.sepval, tree.sepptr[i]:tree.sepptr[i + 1] - 1)
end


"""
    relative(tree::JunctionTree, i::Integer)

Get the inclusion mapping ``\\mathrm{sep}(i) \\to \\mathrm{bag}(parent(i)).`` 
"""
function relative(tree::JunctionTree, i::Integer)
    view(tree.relval, tree.sepptr[i]:tree.sepptr[i + 1] - 1)
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
