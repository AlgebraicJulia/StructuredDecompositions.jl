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


# Find a child w of v such that
# v ∈ snd(w).
# If no such child exists, return nothing.
function findchild(etree::EliminationTree, degree::AbstractVector, stype::Node, v::Integer) end


# Find a child w of v such that
# v ∈ snd(w).
# If no such child exists, return nothing.
function findchild(etree::EliminationTree, degree::AbstractVector, stype::Maximal, v::Integer)
    for w in childindices(etree, v)
        if degree[w] == degree[v] + 1
            return w
        end
    end
end


# Find a child w of v such that
# v ∈ snd(w).
# If no such child exists, return nothing.
function findchild(etree::EliminationTree, degree::AbstractVector, stype::Fundamental, v::Integer)
    ws = childindices(etree, v)

    if length(ws) == 1
        w = only(ws)

        if degree[w] == degree[v] + 1
            return w
        end
    end
end


const DEFAULT_SUPERNODE_TYPE = Maximal()
