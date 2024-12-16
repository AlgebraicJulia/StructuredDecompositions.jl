"""
    SupernodeType

A type of supernode. The options are
- [`Node`](@ref)
- [`Maximal`](@ref)
- [`Fundamental`](@ref)
"""
abstract type SupernodeType end


"""
    Node <: Supernode

A single-vertex supernode.
"""
struct Node <: SupernodeType end


"""
    Maximal <: Supernode

A maximal supernode.
"""
struct Maximal <: SupernodeType end


"""
    Fundamental <: Supernode

A fundamental supernode.
"""
struct Fundamental <: SupernodeType end


# Compact Clique Tree Data Structures in Sparse Matrix Factorizations
# Pothen and Sun
# Figure 4: The Clique Tree Algorithm 2
function pothensun(etree::Tree, colcount::AbstractVector, stype::SupernodeType)
    n = treesize(etree)
    new_in_clique = Vector{Int}(undef, n)
    new = Vector{Int}[]
    parent = Int[]
    first_anc = Int[]

    i = 0

    for v in 1:n
        u = child_in_supernode(etree, colcount, stype, v)
 
        if !isnothing(u)
            new_in_clique[v] = new_in_clique[u]
            push!(new[new_in_clique[v]], v)
        else
            new_in_clique[v] = i += 1
            push!(new, [v])
            push!(parent, 0)
            push!(first_anc, 0)
        end

        for s in childindices(etree, v)
            if s !== u
                parent[new_in_clique[s]] = new_in_clique[v]
                first_anc[new_in_clique[s]] = v
            end
        end
    end

    new, Tree(parent)
end


function stree!(order::Order, graph::OrderedGraph, stype::SupernodeType, etree::Tree=etree!(order, graph))
    rowcount, colcount = supcnt(graph, etree)
    supernode, stree = pothensun(etree, colcount, stype)
    
    sndptr = Vector{Int}(undef, treesize(stree) + 1)
    sepptr = Vector{Int}(undef, treesize(stree) + 1)
    postorder = Order(undef, nv(graph))
    partition = Vector{Int}(undef, nv(graph))
    sndptr[1] = 1
    sepptr[1] = 1

    for (i, j) in enumerate(postorder!(stree))
        residual = supernode[j]
        sndptr[i + 1] = sndptr[i] + length(residual)
        sepptr[i + 1] = sepptr[i] + colcount[last(residual)] - 1
        partition[sndptr[i]:sndptr[i + 1] - 1] .= i
        postorder[sndptr[i]:sndptr[i + 1] - 1] = residual
    end

    permute!(order, postorder)
    permute!(graph, postorder)
    stree, sndptr, sepptr, partition
end


# Find a child w of v such that v ∈ supernode(w).
# If no such child exists, return nothing.
function child_in_supernode(tree::Tree, colcount::AbstractVector, stype::Node, v::Integer) end


# Find a child w of v such that v ∈ supernode(w).
# If no such child exists, return nothing.
function child_in_supernode(tree::Tree, colcount::AbstractVector, stype::Maximal, v::Integer)
    u = nothing

    for w in childindices(tree, v)
        if colcount[w] == colcount[v] + 1
            u = w
            break
        end
    end

    u
end


# Find a child w of v such that v ∈ supernode(w).
# If no such child exists, return nothing.
function child_in_supernode(tree::Tree, colcount::AbstractVector, stype::Fundamental, v::Integer)
    u = nothing

    for w in childindices(tree, v)
        if isnothing(u) && colcount[w] == colcount[v] + 1
            u = w
        else
            u = nothing
            break
        end
    end

    u
end


const DEFAULT_SUPERNODE_TYPE = Maximal()
