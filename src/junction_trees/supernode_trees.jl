"""
    SupernodeTree{V} <: AbstractVector{UnitRange{V}}

A supernodal elimination tree.
This type implements the [indexed tree interface](https://juliacollections.github.io/AbstractTrees.jl/stable/#The-Indexed-Tree-Interface).
"""
struct SupernodeTree{V} <: AbstractVector{UnitRange{V}}
    tree::Tree{V}
    sndptr::Vector{V}

    function SupernodeTree{V}(tree::Tree, sndptr::AbstractVector) where {V}
        # validate parameters
        tree != eachindex(sndptr)[1:end - 1] && throw(ArgumentError("tree != eachindex(sndptr)[1:end - 1]"))
        sndptr[1] != 1 && throw(ArgumentError("sndptr[1] != 1"))

        # construct tree
        new{V}(tree, sndptr)
    end
end


function SupernodeTree(tree::Tree{V}, sndptr::AbstractVector{V}) where V
    SupernodeTree{V}(tree, sndptr)
end


function Tree(tree::SupernodeTree)
    Tree(tree.tree)
end


function Tree{V}(tree::SupernodeTree) where V
    Tree{V}(tree.tree)
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
    lower = reverse(upper)
    rowcount, colcount = supcnt(lower, etree)
    new, ancestor, tree = stree(etree, colcount, snd)
    index = postorder(tree)

    V = vtype(lower)
    E = etype(lower)
    eindex = Vector{V}(undef, length(etree))
    sndptr = Vector{V}(undef, length(tree) + 1)
    sepptr = Vector{E}(undef, length(tree) + 1)
    sndptr[1] = sepptr[1] = 1

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
    lower = reverse(upper)
    rowcount, colcount = supcnt(lower, etree)
    eindex = postorder(etree)

    V = vtype(lower)
    E = etype(lower)
    sndptr = Vector{V}(undef, length(etree) + 1)
    sepptr = Vector{E}(undef, length(etree) + 1)
    sndptr[1] = sepptr[1] = 1

    for (i, j) in enumerate(invperm(eindex))
        sepptr[i + 1] = sepptr[i] + colcount[j] - 1
        sndptr[i + 1] = sndptr[i] + 1
    end

    invpermute!(label, eindex), SupernodeTree(invpermute!(etree, eindex), sndptr), eindex, sepptr, lower, upper
end


function sepdiff(column::AbstractVector{V}, residual::UnitRange{V}) where V
    i = 1

    for v in takewhile(v -> v in residual, column)
        i += 1
    end

    # i = searchsortedlast(column, residual[begin]; by=v -> v ∉ residual) + 1
    @view column[i:end]
end


# Get the separators of every node of a supernodal elimination tree.
function sepval(lower::Graph{V, E}, tree::SupernodeTree{V}, sepptr::Vector{E}) where {V, E}
    function seprange(j)
        sepptr[j]:sepptr[j + 1] - one(E)
    end

    stack = Vector{V}(undef, maximum(length ∘ seprange, eachindex(tree); init=0))
    sepval = Vector{V}(undef, sepptr[end] - 1)

    function separator(j)
        @view sepval[seprange(j)]
    end

    for (j, residual) in enumerate(tree)
        vertex = residual[begin]
        column = sepdiff(neighbors(lower, vertex), residual)
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


#############################
# Abstract Vector Interface #
#############################


function Base.getindex(tree::SupernodeTree{V}, i::Integer) where V
    tree.sndptr[i]:tree.sndptr[i + 1] - one(V)
end


function Base.IndexStyle(::Type{<:SupernodeTree})
    IndexLinear()
end


function Base.size(tree::SupernodeTree)
    size(tree.tree)
end
