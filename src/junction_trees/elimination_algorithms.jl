"""
    EliminationAlgorithm

A graph elimination algorithm. The options are
- [`CuthillMcKeeJL_RCM`](@ref)
- [`AMDJL_AMD`](@ref)
- [`MetisJL_ND`](@ref)
- [`Order`](@ref)
"""
abstract type EliminationAlgorithm end


"""
    CuthillMcKeeJL_RCM <: EliminationAlgorithm

The reverse Cuthill-McKee algorithm. Uses CuthillMckee.jl.
"""
struct CuthillMcKeeJL_RCM <: EliminationAlgorithm end


"""
    AMDJL_AMD <: EliminationAlgorithm

The approximate minimum degree algorithm. Uses AMD.jl.
"""
struct AMDJL_AMD <: EliminationAlgorithm end


"""
    MetisJL_ND <: EliminationAlgorithm

The nested dissection heuristic. Uses Metis.jl.
"""
struct MetisJL_ND <: EliminationAlgorithm end


# Construct an order using the default graph elimination algorithm.
function Order(graph::AbstractSymmetricGraph)
    Order(graph, DEFAULT_ELIMINATION_ALGORITHM)
end


# Construct an order using the reverse Cuthill-McKee algorithm. Uses
# CuthillMcKee.jl.
function Order(graph::AbstractSymmetricGraph, ealg::CuthillMcKeeJL_RCM)
    order = CuthillMcKee.symrcm(adjacencymatrix(graph))
    Order(order)
end


# Construct an order using the approximate minimum degree algorithm. Uses
# AMD.jl.
function Order(graph::AbstractSymmetricGraph, ealg::AMDJL_AMD)
    order = AMD.symamd(adjacencymatrix(graph))
    Order(order)
end


# Construct an order using the nested dissection heuristic. Uses Metis.jl.
function Order(graph::AbstractSymmetricGraph, ealg::MetisJL_ND)
    order, index = Metis.permutation(adjacencymatrix(graph))
    Order(order, index)
end


# Construct the adjacency matrix of a graph.
function adjacencymatrix(graph::AbstractSymmetricGraph)
    m = ne(graph)
    n = nv(graph)

    colptr = ones(Int, n + 1)
    rowval = sizehint!(Vector{Int}(), 2m)

    for j in 1:n
        ns = collect(neighbors(graph, j))
        sort!(ns)
        colptr[j + 1] = colptr[j] + length(ns)
        append!(rowval, ns)
    end

    nzval = ones(Int, length(rowval))
    SparseMatrixCSC(n, n, colptr, rowval, nzval)
end


const DEFAULT_ELIMINATION_ALGORITHM = AMDJL_AMD()
