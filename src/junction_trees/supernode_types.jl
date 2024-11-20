"""
    SupernodeType

A type of supernode. The options are
- [`Node`](@ref)
- [`MaximalSupernode`](@ref)
- [`FundamentalSupernode`](@ref)
"""
abstract type SupernodeType end


"""
    Node <: Supernode

A single-vertex supernode.
"""
struct Node <: SupernodeType end


"""
    MaximalSupernode <: Supernode

A maximal supernode.
"""
struct Maximal <: SupernodeType end


"""
    FundamentalSupernode <: Supernode

A fundamental supernode.
"""
struct Fundamental <: SupernodeType end


# Compact Clique Tree Data Structures in Sparse Matrix Factorizations
# Pothen and Sun
# Figure 4: The Clique Tree Algorithm 2
function stree(etree::EliminationTree, degree::AbstractVector, stype::SupernodeType)
    n = length(etree.tree)
    new_in_clique = Vector{Int}(undef, n)
    new = Vector{Int}[]
    parent = Int[]
    first_anc = Int[]

    i = 0

    for v in 1:n
        u = findchild(etree, degree, stype, v)
 
        if !isnothing(u)
            new_in_clique[v] = new_in_clique[u]
            push!(new[new_in_clique[v]], v)
        else
            new_in_clique[v] = i += 1
            push!(new, [v])
            push!(parent, i)
            push!(first_anc, n)
        end

        for s in childindices(etree.tree, v)
            if s !== u
                parent[new_in_clique[s]] = new_in_clique[v]
                first_anc[new_in_clique[s]] = v
            end
        end
    end

    new, parent, first_anc
end


# Find a child w of v such that
# v ∈ snd(w).
# If no such child exists, return nothing.
function findchild(etree::EliminationTree, degree::AbstractVector, stype::Node, v::Integer) end


# Find a child w of v such that
# v ∈ snd(w).
# If no such child exists, return nothing.
function findchild(etree::EliminationTree, degree::AbstractVector, stype::Maximal, v::Integer)
    for w in childindices(etree.tree, v)
        if degree[w] == degree[v] + 1
            return w
        end
    end
end


# Find a child w of v such that
# v ∈ snd(w).
# If no such child exists, return nothing.
function findchild(etree::EliminationTree, degree::AbstractVector, stype::Fundamental, v::Integer)
    ws = childindices(etree.tree, v)

    if length(ws) == 1
        w = only(ws)

        if degree[w] == degree[v] + 1
            return w
        end
    end
end


const DEFAULT_SUPERNODE_TYPE = Maximal()
