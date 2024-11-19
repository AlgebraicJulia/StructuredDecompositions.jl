"""
    EliminationAlgorithm

A graph elimination algorithm. The options are
- [`CuthillMcKeeJL_RCM`](@ref)
- [`AMDJL_AMD`](@ref)
- [`MetisJL_ND`](@ref)
- [`TreeWidthSolverJL_BT`](@ref)
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
    TreeWidthSolverJL_BT <: EliminationAlgorithm

The Bouchitte-Todinca algorithm. Uses TreeWidthSolver.jl.
"""
struct TreeWidthSolverJL_BT <: EliminationAlgorithm end


"""
    MCS <: EliminationAlgorithm

The maximum cardinality search algorithm.
"""
struct MCS <: EliminationAlgorithm end


# Construct an order using the default graph elimination algorithm.
function Order(graph::AbstractSymmetricGraph)
    Order(graph, DEFAULT_ELIMINATION_ALGORITHM)
end


# Construct an order using the reverse Cuthill-McKee algorithm. Uses CuthillMcKee.jl.
function Order(graph::AbstractSymmetricGraph, ealg::CuthillMcKeeJL_RCM)
    order = CuthillMcKee.symrcm(adjacencymatrix(graph))
    Order(order)
end


# Construct an order using the approximate minimum degree algorithm. Uses AMD.jl.
function Order(graph::AbstractSymmetricGraph, ealg::AMDJL_AMD)
    order = AMD.symamd(adjacencymatrix(graph))
    Order(order)
end


# Construct an order using the nested dissection heuristic. Uses Metis.jl.
function Order(graph::AbstractSymmetricGraph, ealg::MetisJL_ND)
    order, index = Metis.permutation(adjacencymatrix(graph))
    Order(order, index)
end


# Construct an order using the Bouchitte-Todinca algorithm. Uses TreeWidthSolver.jl.
function Order(graph::AbstractSymmetricGraph, ealg::TreeWidthSolverJL_BT)
    n = nv(graph)
    T = TreeWidthSolver.LongLongUInt{n ÷ 64 + 1}
    fadjlist = Vector{Vector{Int}}(undef, n)
    bitfadjlist = Vector{T}(undef, n)
    
    for i in 1:n
        fadjlist[i] = sort(collect(outneighbors(graph, i)))
        bitfadjlist[i] = TreeWidthSolver.bmask(T, fadjlist[i])
    end

    bitgraph = TreeWidthSolver.MaskedBitGraph(bitfadjlist, fadjlist, TreeWidthSolver.bmask(T, 1:n))
    decomposition = TreeWidthSolver.bt_algorithm(bitgraph, TreeWidthSolver.all_pmc_enmu(bitgraph, false), ones(n), false, true)
    order = reverse(vcat(TreeWidthSolver.EliminationOrder(decomposition.tree).order...))
    Order(order)
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


# Simple Linear-Time Algorithms to Test Chordality of Graphs, Test Acyclicity of Hypergraphs, and Selectively Reduce Acyclic Hypergraphs
# Tarjan and Yannakakis
# Maximum Cardinality Search
function mcs(graph::AbstractSymmetricGraph)
    n = nv(graph)
    α = Vector{Int}(undef, n)
    β = Vector{Int}(undef, n)
    size = Vector{Int}(undef, n)
    set = Vector{LinkedLists.LinkedList{Int}}(undef, n)
    pointer = Vector{LinkedLists.ListNode{Int}}(undef, n)

    for i in 1:n
        size[i] = 1
        set[i] = LinkedLists.LinkedList{Int}()
        pointer[i] = push!(set[1], i)
    end

    i = n
    j = 1

    while i >= 1
        v = first(set[j])
        deleteat!(set[j], pointer[v])        
        α[v] = i
        β[i] = v
        size[v] = 0

        for w in neighbors(graph, v)
            if size[w] >= 1
                deleteat!(set[size[w]], pointer[w])
                size[w] += 1
                pointer[w] = push!(set[size[w]], w)
            end
        end

        i -= 1
        j += 1

        while j >= 1 && isempty(set[j])
            j -= 1
        end
    end

    β, α
end


const DEFAULT_ELIMINATION_ALGORITHM = AMDJL_AMD()
