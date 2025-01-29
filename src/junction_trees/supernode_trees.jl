"""
    SupernodeTree{I} <: AbstractVector{UnitRange{I}}

A supernodal elimination tree.
This type implements the [indexed tree interface](https://juliacollections.github.io/AbstractTrees.jl/stable/#The-Indexed-Tree-Interface).
"""
struct SupernodeTree{I} <: AbstractVector{UnitRange{I}}
    tree::Tree{I}
    sndptr::Vector{I}

    function SupernodeTree{I}(tree::Tree, sndptr::AbstractVector) where I
        # validate parameters
        tree != eachindex(sndptr)[1:end - 1] && throw(ArgumentError("tree != eachindex(sndptr)[1:end - 1]"))
        sndptr[1] != 1 && throw(ArgumentError("sndptr[1] != 1"))

        # construct tree
        new{I}(tree, sndptr)
    end
end


function SupernodeTree(tree::Tree{I}, sndptr::AbstractVector{I}) where I
    SupernodeTree{I}(tree, sndptr)
end


function Tree(tree::SupernodeTree)
    Tree(tree.tree)
end


function Tree{I}(tree::SupernodeTree) where I
    Tree{I}(tree.tree)
end


"""
    supernodetree(graph;
        alg::PermutationOrAlgorithm=DEFAULT_ELIMINATION_ALGORITHM,
        snd::SupernodeType=DEFAULT_SUPERNODE_TYPE)

Construct a supernodal elimination tree.
"""
function supernodetree(graph; alg::PermutationOrAlgorithm=DEFAULT_ELIMINATION_ALGORITHM, snd::SupernodeType=DEFAULT_SUPERNODE_TYPE)
    supernodetree(graph, alg, snd)
end


function supernodetree(graph, alg::PermutationOrAlgorithm, snd::SupernodeType)
    label, tree, upper = eliminationtree(graph, alg)
    cache = spzeros(Nothing, indtype(upper), size(upper))
    supernodetree!(label, tree, upper, cache, snd)
end


"""
    suoernodetree!(graph;
        alg::PermutationOrAlgorithm=DEFAULT_ELIMINATION_ALGORITHM,
        snd::SupernodeType=DEFAULT_SUPERNODE_TYPE)

A mutating version of [`supernodetree!`](@ref).
"""
function supernodetree!(graph; alg::PermutationOrAlgorithm=DEFAULT_ELIMINATION_ALGORITHM, snd::SupernodeType=DEFAULT_SUPERNODE_TYPE)
    label, tree, index, sepptr, lower, cache = supernodetree!(graph, alg, snd)
    label, tree
end


function supernodetree!(matrix::SparseMatrixCSC{<:Any, I}, alg::PermutationOrAlgorithm, snd::SupernodeType) where I
    label, tree, upper = eliminationtree(matrix, alg)
    cache = SparseMatrixCSC{Nothing, I}(size(matrix)..., getcolptr(matrix), rowvals(matrix), Vector{Nothing}(undef, nnz(matrix)))
    supernodetree!(label, tree, upper, cache, snd)
end


function supernodetree!(label::Vector{I}, etree::Tree{I}, upper::SparseMatrixCSC{Nothing, I}, cache::SparseMatrixCSC{Nothing, I}, snd::SupernodeType) where I
    lower = transpose!(cache, upper)
    rowcount, colcount = supcnt(lower, etree)
    new, ancestor, tree = stree(etree, colcount, snd)
    index = postorder(tree)

    eindex = Vector{I}(undef, size(lower, 2))
    sepptr = sizehint!(I[], length(tree) + 1)
    sndptr = sizehint!(I[], length(tree) + 1)
    push!(sndptr, 1)
    push!(sepptr, 1)

    for j in invperm(index)
        u = new[j]
        w = ancestor[j]
        p = eindex[u] = sndptr[end]

        for v in ancestorindices(etree, u)
            if v != w
                p += 1
                eindex[v] = p
            else
                break
            end
        end

        push!(sepptr, sndptr[end] + sepptr[end] + colcount[u] - p - 1)
        push!(sndptr, p + 1)
    end

    invpermute!(label, eindex), SupernodeTree(invpermute!(tree, index), sndptr), eindex, sepptr, lower, upper
end


function sepdiff(column::AbstractVector{I}, residual::AbstractVector{I}) where I
    vertex = residual[begin]
    @view column[searchsortedlast(column, vertex; by=v -> v ∉ residual) + 1:end]
end


# Get the separators of every node of a supernodal elimination tree.
function sepvals(lower::SparseMatrixCSC{Nothing, I}, tree::SupernodeTree{I}, sepptr::Vector{I}) where I
    cache = sizehint!(I[], maximum(i -> sepptr[i + 1] - sepptr[i], eachindex(tree); init=zero(I)))
    sepval = sizehint!(I[], sepptr[end] - 1)

    for (j, residual) in enumerate(tree)
        vertex = residual[begin]
        column = @view rowvals(lower)[nzrange(lower, vertex)]
        append!(sepval, sepdiff(column, residual))

        for i in childindices(tree, j)
            state = @view sepval[sepptr[j]:end]
            child = @view sepval[sepptr[i]:sepptr[i + 1] - 1]
            mergesorted!(empty!(cache), state, sepdiff(child, residual))
            append!(resize!(sepval, sepptr[j] - 1), cache)
        end
    end

    sepval
end


# Get the relative indices of every node of a supernodal elimination tree.
function relvals(tree::SupernodeTree{I}, sepptr::Vector{I}, sepval::Vector{I}) where I
    relval = Vector{I}(undef, length(sepval))

    for (j, residual) in enumerate(tree)
        for i in childindices(tree, j)
            p = sepptr[i]
            q = sepptr[j]

            while p < sepptr[i + 1] && sepval[p] in residual
                relval[p] = sepval[p] - residual[begin] + 1
                p += 1
            end

            while p < sepptr[i + 1] && q < sepptr[j + 1]
                if sepval[p] <= sepval[q]
                    relval[p] = length(residual) + q - sepptr[j] + 1
                    p += 1
                end

                q += 1
            end
        end
    end

    relval
end


##########################
# Indexed Tree Interface #
##########################


function AbstractTrees.rootindex(tree::SupernodeTree)
    rootindex(tree.tree)
end


function AbstractTrees.parentindex(tree::SupernodeTree, i::Integer)
    parentindex(tree.tree, i)
end


function firstchildindex(tree::SupernodeTree, i::Integer)
    firstchildindex(tree.tree, i)
end


function AbstractTrees.nextsiblingindex(tree::SupernodeTree, i::Integer)
    nextsiblingindex(tree.tree, i)
end


function rootindices(tree::SupernodeTree)
    rootindices(tree.tree)
end


function AbstractTrees.childindices(tree::SupernodeTree, i::Integer)
    childindices(tree.tree, i)
end


function ancestorindices(tree::SupernodeTree, i::Integer)
    ancestorindices(tree.tree, i)
end


function setrootindex!(tree::SupernodeTree, i::Integer)
    setrootindex!(tree.tree, i) 
end


#############################
# Abstract Vector Interface #
#############################


function Base.getindex(tree::SupernodeTree{I}, i::Integer) where I
    tree.sndptr[i]:tree.sndptr[i + 1] - one(I)
end


function Base.IndexStyle(::Type{SupernodeTree})
    IndexLinear()
end


function Base.size(tree::SupernodeTree)
    size(tree.tree)
end
