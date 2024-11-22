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


"""
    Order(graph[, ealg::EliminationAlgorithm])

Construct an elimination order for a simple graph, optionally specifying an elimination algorithm.
"""
function Order(graph, ealg::EliminationAlgorithm=DEFAULT_ELIMINATION_ALGORITHM)
    Order(adjacencymatrix(graph), ealg)
end


# Construct an elimination order.
function Order(graph::AbstractMatrix)
    Order(graph, DEFAULT_ELIMINATION_ALGORITHM)
end


# Construct an order using the reverse Cuthill-McKee algorithm. Uses CuthillMcKee.jl.
function Order(graph::AbstractMatrix, ealg::CuthillMcKeeJL_RCM)
    order = CuthillMcKee.symrcm(graph)
    Order(order)
end


# Construct an order using the approximate minimum degree algorithm. Uses AMD.jl.
function Order(graph::AbstractMatrix, ealg::AMDJL_AMD)
    order = AMD.symamd(graph)
    Order(order)
end


# Construct an order using the nested dissection heuristic. Uses Metis.jl.
function Order(graph::AbstractMatrix, ealg::MetisJL_ND)
    order, index = Metis.permutation(graph)
    Order(order, index)
end


# Construct an order using the Bouchitte-Todinca algorithm. Uses TreeWidthSolver.jl.
function Order(graph::AbstractSparseMatrixCSC, ealg::TreeWidthSolverJL_BT)
    n = size(graph, 1)
    T = TreeWidthSolver.LongLongUInt{n ÷ 64 + 1}
    fadjlist = Vector{Vector{Int}}(undef, n)
    bitfadjlist = Vector{T}(undef, n)
    
    for i in 1:n
        fadjlist[i] = rowvals(graph)[nzrange(graph, i)]
        bitfadjlist[i] = TreeWidthSolver.bmask(T, fadjlist[i])
    end

    bitgraph = TreeWidthSolver.MaskedBitGraph(bitfadjlist, fadjlist, TreeWidthSolver.bmask(T, 1:n))
    decomposition = TreeWidthSolver.bt_algorithm(bitgraph, TreeWidthSolver.all_pmc_enmu(bitgraph, false), ones(n), false, true)
    order = reverse(vcat(TreeWidthSolver.EliminationOrder(decomposition.tree).order...))
    Order(order)
end


# Construct an order using the maximum cardinality search algorithm.
function Order(graph::AbstractMatrix, ealg::MCS)
    order, index = mcs(graph)
    Order(order, index)
end


# Construct the adjacency matrix of a graph.
function adjacencymatrix(graph::BasicGraphs.AbstractSymmetricGraph)
    m = BasicGraphs.ne(graph)
    n = BasicGraphs.nv(graph)
    colptr = Vector{Int}(undef, n + 1)
    rowval = Vector{Int}(undef, m)
    count = 1

    for i in 1:n
        colptr[i] = count
        neighbor = collect(BasicGraphs.all_neighbors(graph, i))
        sort!(neighbor)

        for j in neighbor
            rowval[count] = j
            count += 1
        end
    end

    colptr[n + 1] = m + 1
    nzval = ones(Bool, m)
    SparseMatrixCSC(n, n, colptr, rowval, nzval)
end


# Construct the adjacency matrix of a graph.
function adjacencymatrix(graph::AbstractGraph)
    adjacency_matrix(graph; dir=:both)
end


# Simple Linear-Time Algorithms to Test Chordality of Graphs, Test Acyclicity of Hypergraphs, and Selectively Reduce Acyclic Hypergraphs
# Tarjan and Yannakakis
# Maximum Cardinality Search
function mcs(graph::AbstractSparseMatrixCSC)
    n = size(graph, 1)
    α = Vector{Int}(undef, n)
    β = Vector{Int}(undef, n)
    len = Vector{Int}(undef, n)
    set = Vector{LinkedLists.LinkedList{Int}}(undef, n)
    pointer = Vector{LinkedLists.ListNode{Int}}(undef, n)

    for i in 1:n
        len[i] = 1
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
        len[v] = 0

        for w in view(rowvals(graph), nzrange(graph, v))
            if len[w] >= 1
                deleteat!(set[len[w]], pointer[w])
                len[w] += 1
                pointer[w] = push!(set[len[w]], w)
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
