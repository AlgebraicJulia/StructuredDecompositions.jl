"""
    EliminationAlgorithm

A graph elimination algorithm. The options are
- [`CuthillMcKeeJL_RCM`](@ref)
- [`SymRCMJL_RCM`](@ref)
- [`AMDJL_AMD`](@ref)
- [`AMDJL_SYMAMD`](@ref)
- [`MetisJL_ND`](@ref)
- [`TreeWidthSolverJL_BT`](@ref)
- [`MCS`](@ref)
"""
abstract type EliminationAlgorithm end


"""
    PermutationOrAlgorithm = Union{AbstractVector, EliminationAlgorithm}

Either a permutation or a graph elimination algorithm.
"""
const PermutationOrAlgorithm = Union{AbstractVector, EliminationAlgorithm}


"""
    CuthillMcKeeJL_RCM <: EliminationAlgorithm

The reverse Cuthill-McKee algorithm. Uses CuthillMckee.jl.
"""
struct CuthillMcKeeJL_RCM <: EliminationAlgorithm end


"""
    SymRCMJL_RCM <: EliminationAlgorithm

The reverse Cuthill-McKee algorithm. Uses SymRCM.jl.
"""
struct SymRCMJL_RCM <: EliminationAlgorithm end


"""
    AMDJL_AMD <: EliminationAlgorithm

The approximate minimum degree algorithm. Uses AMD.jl.
"""
struct AMDJL_AMD <: EliminationAlgorithm
    meta::AMD.Amd
end

"""
    AMDJL_SYMAMD{T} <: EliminationAlgorithm

The SYMAMD algorithm. Uses AMD.jl.
"""
struct AMDJL_SYMAMD{T} <: EliminationAlgorithm
    meta::AMD.Colamd{T}
end


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


function AMDJL_AMD()
    AMDJL_AMD(AMD.Amd())
end


function AMDJL_SYMAMD{T}() where T
    AMDJL_SYMAMD(AMD.Colamd{T}())
end


function AMDJL_SYMAMD()
    AMDJL_SYMAMD{Int}()
end


"""
    permutation(matrix::SparseMatrixCSC, alg::PermutationOrAlgorithm)

Construct a [fill-reducing permutation](https://www.mathworks.com/help/matlab/math/sparse-matrix-reordering.html)
``\\sigma: V \\to V`` of the vertices of a [simple graph](https://mathworld.wolfram.com/SimpleGraph.html) ``G = (V, E)``, represented by its adjacency matrix.
Returns a tuple ``(\\sigma, \\sigma^{-1})``.
"""
permutation(matrix::SparseMatrixCSC, alg::PermutationOrAlgorithm)


function permutation(matrix::SparseMatrixCSC, alg::AbstractVector)
    order = collect(alg)
    order, invperm(order)
end


# Construct an order using the reverse Cuthill-McKee algorithm. Uses CuthillMcKee.jl.
function permutation(matrix::SparseMatrixCSC, alg::CuthillMcKeeJL_RCM)
    order = CuthillMcKee.symrcm(matrix)
    order, invperm(order)
end


# Construct an order using the reverse Cuthill-McKee algorithm. Uses SYMRCMJL_RCM.jl.
function permutation(matrix::SparseMatrixCSC, alg::SymRCMJL_RCM)
    order = SymRCM.symrcm(matrix)
    order, invperm(order)
end


# Construct an order using the approximate minimum degree algorithm. Uses AMD.jl.
function permutation(matrix::SparseMatrixCSC, alg::AMDJL_AMD)
    order = AMD.amd(matrix, alg.meta)
    order, invperm(order)
end


# Construct an order using the SYMAMD algorithm. Uses AMD.jl.
function permutation(matrix::SparseMatrixCSC, alg::AMDJL_SYMAMD)
    order = AMD.symamd(matrix, alg.meta)
    order, invperm(order)
end


# Construct an order using the nested dissection heuristic. Uses Metis.jl.
function permutation(matrix::SparseMatrixCSC, alg::MetisJL_ND)
    order, index = Metis.permutation(matrix)
    order, index
end


# Construct an order using the Bouchitte-Todinca algorithm. Uses TreeWidthSolver.jl.
function permutation(matrix::SparseMatrixCSC, alg::TreeWidthSolverJL_BT)
    T = TreeWidthSolver.LongLongUInt{size(matrix, 2) ÷ 64 + 1}
    fadjlist = Vector{Vector{Int}}(undef, size(matrix, 2))
    bitgraph = Vector{T}(undef, size(matrix, 2))

    for j in axes(matrix, 2)
        neighbors = getindex(rowvals(matrix), nzrange(matrix, j))
        fadjlist[j] = neighbors
        bitgraph[j] = TreeWidthSolver.bmask(T, neighbors)
    end

    graph = TreeWidthSolver.MaskedBitGraph(bitgraph, fadjlist, TreeWidthSolver.bmask(T, axes(matrix, 2)))
    cliques = TreeWidthSolver.all_pmc_enmu(graph, false)
    decomposition = TreeWidthSolver.bt_algorithm(graph, cliques, ones(size(matrix, 2)), false, true)
    order = reverse!(reduce(vcat, TreeWidthSolver.EliminationOrder(decomposition.tree).order))
    order, invperm(order)
end


# Construct an order using the maximum cardinality search algorithm.
function permutation(matrix::SparseMatrixCSC, alg::MCS)
    order = mcs(matrix)
    order, invperm(order)
end


const DEFAULT_ELIMINATION_ALGORITHM = AMDJL_AMD()
