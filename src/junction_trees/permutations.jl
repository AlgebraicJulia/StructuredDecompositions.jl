"""
    Permutation <: AbstractVector{Int}

A [permutation](https://en.wikipedia.org/wiki/Permutation) ``\\sigma`` of the set ``\\{1, \\dots, n\\}``.
This type implements the [abstract array interface](https://docs.julialang.org/en/v1/manual/interfaces/#man-interface-array).
"""
struct Permutation <: AbstractVector{Int}
    inner::Vector{Int} # permutation
    index::Vector{Int} # inverse permutation
end


"""
    Permutation(permutation::AbstractVector)

Construct a permutation ``\\sigma`` from a sequence ``(\\sigma(1), \\dots, \\sigma(n))``.
"""
function Permutation(permutation::AbstractVector)
    Permutation(permutation, invperm(permutation))
end


"""
    Permutation(matrix::SparseMatrixCSC, alg::Union{AbstractVector, EliminationAlgorithm})

Construct a fill-reducing permutation of the vertices of a graph.
"""
Permutation(matrix::SparseMatrixCSC, alg::Union{AbstractVector, EliminationAlgorithm})


function Permutation(matrix::SparseMatrixCSC, alg::AbstractVector)
    Permutation(deepcopy(alg))
end


# Construct an order using the reverse Cuthill-McKee algorithm. Uses CuthillMcKee.jl.
function Permutation(matrix::SparseMatrixCSC, alg::CuthillMcKeeJL_RCM)
    order = CuthillMcKee.symrcm(matrix)
    Permutation(order)
end


# Construct an order using the reverse Cuthill-McKee algorithm. Uses SYMRCMJL_RCM.jl.
function Permutation(matrix::SparseMatrixCSC, alg::SymRCMJL_RCM)
    order = SymRCM.symrcm(matrix)
    Permutation(order)
end


# Construct an order using the approximate minimum degree algorithm. Uses AMD.jl.
function Permutation(matrix::SparseMatrixCSC, alg::AMDJL_AMD)
    order = AMD.amd(matrix, alg.meta)
    Permutation(order)
end


# Construct an order using the SYMAMD algorithm. Uses AMD.jl.
function Permutation(matrix::SparseMatrixCSC, alg::AMDJL_SYMAMD)
    order = AMD.symamd(matrix, alg.meta)
    Permutation(order)
end


# Construct an order using the nested dissection heuristic. Uses Metis.jl.
function Permutation(matrix::SparseMatrixCSC, alg::MetisJL_ND)
    order, index = Metis.permutation(matrix)
    Permutation(order, index)
end


# Construct an order using the Bouchitte-Todinca algorithm. Uses TreeWidthSolver.jl.
function Permutation(matrix::SparseMatrixCSC, alg::TreeWidthSolverJL_BT)
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
    Permutation(order)
end


# Construct an order using the maximum cardinality search algorithm.
function Permutation(matrix::SparseMatrixCSC, alg::MCS)
    order = mcs(matrix)
    Permutation(order)
end


function permutation(matrix::SparseMatrixCSC, alg::Union{AbstractVector, EliminationAlgorithm})
    order = Permutation(matrix, alg)
    order.inner, order.index
end


function Base.invperm(permutation::Permutation)
    Permutation(permutation.index, permutation.inner)
end


function Base.permute!(permutation::Permutation, other::AbstractVector)
    permute!(permutation.inner, other)
    permutation.index[permutation.inner] = eachindex(permutation.inner)
    permutation
end


function Base.invpermute!(permutation::Permutation, other::AbstractVector)
    invpermute!(permutation.inner, other)
    permutation.index[permutation.inner] = eachindex(permutation.inner)
    permutation
end


function Base.convert(::Type{Vector{Int}}, permutation::Permutation)
    permutation.inner
end


function Base.copy(permutation::Permutation)
    Permutation(permutation.inner, permutation.index)
end


function Base.deepcopy(permutation::Permutation)
    Permutation(copy(permutation.inner), copy(permutation.index))
end


#############################
# Abstract Vector Interface #
#############################


function Base.getindex(permutation::Permutation, i)
    permutation.inner[i]
end


function Base.IndexStyle(::Type{Permutation})
    IndexLinear()
end


function Base.size(permutation::Permutation)
    (length(permutation.inner),)
end
