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


# Chordal Graphs and Semidefinite Optimization
# Vanderberghe and Andersen
# Algorithm 4.1: Maximal supernodes and supernodal elimination tree.
function stree(etree::EliminationTree, degree::AbstractVector, stype::SupernodeType)
    m = 0
    n = length(etree.tree)
    index = zeros(Int, n)
    snd = Vector{Int}[]
    q = Int[]
    a = Int[]

    for v in 1:n
        ww = findchild(etree, degree, stype, v)
        
        if isnothing(ww)
            index[v] = i = m += 1
            push!(snd, [v])
            push!(q, m)
            push!(a, n)
        else
            index[v] = i = index[ww]
            push!(snd[i], v)
        end

        for w in childindices(etree.tree, v)
            if w !== ww
                j = index[w]
                q[j] = i
                a[j] = v
            end
        end
    end

    snd, q, a
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
