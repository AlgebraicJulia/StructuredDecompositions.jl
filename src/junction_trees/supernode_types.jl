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
function stree(etree::EliminationTree, degree::AbstractVector, stype::SupernodeType)
    n = treesize(etree.tree)
    new_in_clique = Vector{Int}(undef, n)
    new = Vector{Int}[]
    parent = Int[]
    first_anc = Int[]

    i = 0

    for v in 1:n
        u = child_in_supernode(etree, degree, stype, v)
 
        if !isnothing(u)
            new_in_clique[v] = new_in_clique[u]
            push!(new[new_in_clique[v]], v)
        else
            new_in_clique[v] = i += 1
            push!(new, [v])
            push!(parent, 0)
            push!(first_anc, n)
        end

        for s in childindices(etree.tree, v)
            if s !== u
                parent[new_in_clique[s]] = new_in_clique[v]
                first_anc[new_in_clique[s]] = v
            end
        end
    end

    new_in_clique, new, parent, first_anc
end


# Find a child w of v such that v ∈ supernode(w).
# If no such child exists, return nothing.
function child_in_supernode(etree::EliminationTree, degree::AbstractVector, stype::Node, v::Integer) end


# Find a child w of v such that v ∈ supernode(w).
# If no such child exists, return nothing.
function child_in_supernode(etree::EliminationTree, degree::AbstractVector, stype::Maximal, v::Integer)
    u = nothing

    for w in childindices(etree.tree, v)
        if degree[w] == degree[v] + 1
            u = w
            break
        end
    end

    u
end


# Find a child w of v such that v ∈ supernode(w).
# If no such child exists, return nothing.
function child_in_supernode(etree::EliminationTree, degree::AbstractVector, stype::Fundamental, v::Integer)
    u = nothing

    for w in childindices(etree.tree, v)
        if isnothing(u) && degree[w] == degree[v] + 1
            u = w
        else
            u = nothing
            break
        end
    end

    u
end


const DEFAULT_SUPERNODE_TYPE = Maximal()
