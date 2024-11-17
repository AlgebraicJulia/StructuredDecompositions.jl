"""
    EliminationAlgorithm

A graph elimination algorithm. The options are
- [`CuthillMcKeeJL_RCM`](@ref)
- [`AMDJL_AMD`](@ref)
- [`MetisJL_ND`](@ref)
- [`MCS`](@ref)
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


"""
    MCS <: EliminationAlgorithm

The maximum cardinality search algorithm.
"""
struct MCS <: EliminationAlgorithm end


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


# Concstruct an order using the maximum cardinality search algorithm.
function Order(graph::AbstractSymmetricGraph, ealg::MCS)
    order, index = mcs(graph)
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


# Simple Linear-Time Algorithms to Test Chordality of Graphs, Test Acyclicity
# of Hypergraphs, and Selectively Reduce Acyclic Hypergraphs
#
# Tarjan and Yannakakis
#
# Maximum Cardinality Search
function mcs(graph::AbstractSymmetricGraph)
    n = nv(graph)
    α = Vector{Int}(undef, n)
    α⁻¹ = Vector{Int}(undef, n)
    size = Vector{Int}(undef, n)
    set = Vector{Set{Int}}(undef, n)

    for i in 1:n
        size[i] = 1
        set[i] = Set()
    end

    union!(set[1], 1:n)

    i = n
    j = 1

    while i >= 1
        v = pop!(set[j])
        α[v] = i
        α⁻¹[i] = v
        size[v] = 0

        for w in neighbors(graph, v)
            if size[w] >= 1
                delete!(set[size[w]], w)
                size[w] += 1
                push!(set[size[w]], w)
            end
        end

        i -= 1
        j += 1

        while j >= 1 && isempty(set[j])
            j -= 1
        end
    end

    α⁻¹, α
end

const DEFAULT_ELIMINATION_ALGORITHM = AMDJL_AMD()
