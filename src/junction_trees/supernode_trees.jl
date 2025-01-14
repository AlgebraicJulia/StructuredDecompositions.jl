struct SupernodeTree <: AbstractVector{UnitRange{Int}}
    tree::Tree
    sndptr::Vector{Int}

    function SupernodeTree(tree::Tree, sndptr::AbstractVector)
        # validate parameters
        tree != eachindex(sndptr)[1:end - 1] && throw(ArgumentError("tree != eachindex(sndptr)[1:end - 1]"))
        sndptr[1] != 1 && throw(ArgumentError("sndptr[1] != 1"))

        # construct tree
        new(tree, sndptr)
    end
end


function SupernodeTree(tree::SupernodeTree)
    SupernodeTree(tree.tree, tree.sndptr)
end


function Tree(tree::SupernodeTree)
    Tree(tree.tree)
end


function Tree(tree::SupernodeTree, root::Integer)
    Tree(tree.tree, root)
end


function supernodetree(matrix::AbstractMatrix; alg::PermutationOrAlgorithm=DEFAULT_ELIMINATION_ALGORITHM, snd::SupernodeType=DEFAULT_SUPERNODE_TYPE)
    supernodetree!(sparse(matrix); alg, snd)
end


function supernodetree!(matrix::SparseMatrixCSC; alg::PermutationOrAlgorithm=DEFAULT_ELIMINATION_ALGORITHM, snd::SupernodeType=DEFAULT_SUPERNODE_TYPE)
    label, tree, index, sepptr, lower, cache = supernodetree!(matrix, alg, snd)
    label, tree
end


# Construct a postordered supernodal elimination tree.
function supernodetree!(matrix::SparseMatrixCSC, alg::PermutationOrAlgorithm, snd::SupernodeType)
    label, etree, upper, cache = eliminationtree!(matrix, alg)
    lower = transpose!(cache, upper)
    rowcount, colcount = supcnt(lower, etree)
    new, ancestor, tree = stree(etree, colcount, snd)
    index = dfs(tree)

    eindex = Vector{Int}(undef, size(lower, 1))
    sepptr = sizehint!(Int[], length(tree) + 1)
    sndptr = sizehint!(Int[], length(tree) + 1)
    push!(sndptr, 1)
    push!(sepptr, 1)

    for j in invperm(index)
        u = v = new[j]
        p = last(sndptr)

        while !isnothing(v) && v != ancestor[j]
            eindex[v] = p
            v = parentindex(etree, v)
            p += 1
        end

        push!(sepptr, last(sndptr) + last(sepptr) + colcount[u] - p)
        push!(sndptr, p)
    end

    invpermute!(label, eindex), SupernodeTree(invpermute!(tree, index), sndptr), eindex, sepptr, lower, upper
end


# Get the separators of every node of a supernodal elimination tree.
function sepvals(lower::SparseMatrixCSC, tree::SupernodeTree, sepptr::AbstractVector)
    cache = sizehint!(Int[], maximum(i -> sepptr[i + 1] - sepptr[i], eachindex(tree)))
    sepval = sizehint!(Int[], last(sepptr) - 1)

    for (j, residual) in enumerate(tree)
        append!(sepval, ifilter(v -> v ∉ residual, view(rowvals(lower), nzrange(lower, first(residual)))))

        for i in childindices(tree, j)
            state = @view sepval[sepptr[j]:end]
            mergesorted!(empty!(cache), state, ifilter(v -> v ∉ residual, view(sepval, sepptr[i]:sepptr[i + 1] - 1)))
            append!(resize!(sepval, sepptr[j] - 1), cache)
        end
    end

    sepval
end


# Get the relative indices of every node of a supernodal elimination tree.
function relvals(tree::SupernodeTree, sepptr::AbstractVector, sepval::AbstractVector)
    relval = Vector{Int}(undef, length(sepval))

    for (j, residual) in enumerate(tree)
        for i in childindices(tree, j)
            p = sepptr[i]
            q = sepptr[j]

            while p < sepptr[i + 1] && sepval[p] in residual
                relval[p] = sepval[p] - first(residual) + 1
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


function Base.show(io::IO, ::MIME"text/plain", tree::SupernodeTree)
    println(io, "$(length(tree))-element SupernodeTree:")
    print_tree(io, IndexNode(tree))
end


###########################
# Abstract Tree Interface #
###########################


function firstchildindex(tree::SupernodeTree, i::Integer)
    firstchildindex(tree.tree, i)
end


function AbstractTrees.rootindex(tree::SupernodeTree)
    rootindex(tree.tree)
end


function AbstractTrees.parentindex(tree::SupernodeTree, i::Integer)
    parentindex(tree.tree, i)
end


function AbstractTrees.nextsiblingindex(tree::SupernodeTree, i::Integer)
    nextsiblingindex(tree.tree, i)
end


function AbstractTrees.childindices(tree::SupernodeTree, i::Integer)
    childindices(tree.tree, i)
end


function AbstractTrees.ParentLinks(::Type{IndexNode{SupernodeTree, Int}})
    StoredParents()
end


function AbstractTrees.SiblingLinks(::Type{IndexNode{SupernodeTree, Int}})
    StoredSiblings()
end


function AbstractTrees.NodeType(::Type{IndexNode{SupernodeTree, Int}})
    HasNodeType()
end


function AbstractTrees.nodetype(::Type{IndexNode{SupernodeTree, Int}})
    IndexNode{SupernodeTree, Int}
end


#############################
# Abstract Vector Interface #
#############################


function Base.getindex(tree::SupernodeTree, i::Integer)
    tree.sndptr[i]:tree.sndptr[i + 1] - 1
end


function Base.IndexStyle(::Type{SupernodeTree})
    IndexLinear()
end


function Base.size(tree::SupernodeTree)
    size(tree.tree)
end
