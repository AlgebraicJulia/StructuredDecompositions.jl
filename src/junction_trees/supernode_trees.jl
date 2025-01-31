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
    label, etree, upper = eliminationtree(graph, alg)
    lower = copy(transpose(upper))
    rowcount, colcount = supcnt(lower, etree)
    new, ancestor, tree = stree(etree, colcount, snd)
    index = postorder(tree)

    I = indtype(lower)
    eindex = Vector{I}(undef, size(lower, 2))
    sepptr = Vector{I}(undef, length(tree) + 1)
    sndptr = Vector{I}(undef, length(tree) + 1)
    sepptr[1] = sndptr[1] = 1

    for (i, j) in enumerate(invperm(index))
        u = new[j]
        p = eindex[u] = sndptr[i]

        for v in takewhile(v -> v != ancestor[j], ancestorindices(etree, u))
            eindex[v] = p += 1
        end

        sepptr[i + 1] = sndptr[i] + sepptr[i] + colcount[u] - p - 1
        sndptr[i + 1] = p + 1
    end

    invpermute!(label, eindex), SupernodeTree(invpermute!(tree, index), sndptr), eindex, sepptr, lower, upper
end


function supernodetree(graph, alg::PermutationOrAlgorithm, snd::Nodal)
    label, etree, upper = eliminationtree(graph, alg)
    lower = copy(transpose(upper))
    rowcount, colcount = supcnt(lower, etree)
    eindex = postorder(etree)

    I = indtype(lower)
    sepptr = Vector{I}(undef, length(etree) + 1)
    sndptr = Vector{I}(undef, length(etree) + 1)
    sepptr[1] = sndptr[1] = 1

    for (i, j) in enumerate(invperm(eindex))
        sepptr[i + 1] = sepptr[i] + colcount[j] - 1
        sndptr[i + 1] = sndptr[i] + 1
    end

    invpermute!(label, eindex), SupernodeTree(invpermute!(etree, eindex), sndptr), eindex, sepptr, lower, upper
end


function sepdiff(column::AbstractVector{I}, residual::AbstractVector{I}) where I
    i = 1

    for v in takewhile(v -> v in residual, column)
        i += 1
    end

    # i = searchsortedlast(column, residual[begin]; by=v -> v ∉ residual) + 1
    @view column[i:end]
end


# Get the separators of every node of a supernodal elimination tree.
function sepval(lower::SparseMatrixCSC{Bool, I}, tree::SupernodeTree{I}, sepptr::Vector{I}) where I
    function seprange(j)
        sepptr[j]:sepptr[j + 1] - one(I)
    end

    function neighbors(j)
        @view rowvals(lower)[nzrange(lower, j)]
    end

    stack = Vector{I}(undef, maximum(length ∘ seprange, eachindex(tree); init=0))
    sepval = Vector{I}(undef, sepptr[end] - 1)

    function separator(j)
        @view sepval[seprange(j)]
    end

    for (j, residual) in enumerate(tree)
        column = sepdiff(neighbors(residual[begin]), residual)
        state = @view separator(j)[eachindex(column)]
        copy!(state, column)

        for i in childindices(tree, j)
            child = mergesorted!(stack, state, sepdiff(separator(i), residual))
            state = @view separator(j)[eachindex(child)]
            copy!(state, child)
        end
    end

    sepval
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
    tree
end


#############################
# Abstract Vector Interface #
#############################


function Base.getindex(tree::SupernodeTree{I}, i::Integer) where I
    tree.sndptr[i]:tree.sndptr[i + 1] - one(I)
end


function Base.IndexStyle(::Type{<:SupernodeTree})
    IndexLinear()
end


function Base.size(tree::SupernodeTree)
    size(tree.tree)
end
