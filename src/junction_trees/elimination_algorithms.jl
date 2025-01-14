"""
    EliminationAlgorithm

A graph elimination algorithm. The options are

| type                 | name                              | complexity  |
| :------------------- | :-------------------------------- | :---------- |
| [`MCS`](@ref)        | maximum cardinality search        | O(m + n)    |
| [`RCM`](@ref)        | reverse Cuthill-Mckee             | O(mΔ)       |
| [`AMD`](@ref)        | approximate minimum degree        | O(mn)       |
| [`SymAMD`](@ref)     | column approximate minimum degree | O(mn)       |
| [`MMD`](@ref)        | multiple minimum degree           | O(mn²)      |
| [`NodeND`](@ref)     | nested dissection                 |             |
| [`FlowCutter`](@ref) | FlowCutter                        |             |
| [`Spectral`](@ref)   | spectral ordering                 |             |
| [`BT`](@ref)         | Bouchitte-Todinca                 | O(2.6183ⁿ)  |

for a graph with m edges, n vertices, and maximum degree Δ.
"""
abstract type EliminationAlgorithm end


"""
    PermutationOrAlgorithm = Union{AbstractVector, EliminationAlgorithm}

Either a permutation or a graph elimination algorithm.
"""
const PermutationOrAlgorithm = Union{AbstractVector, EliminationAlgorithm}


"""
    MCS <: EliminationAlgorithm

    MCS()

The maximum cardinality search algorithm.
"""
mutable struct MCS <: EliminationAlgorithm end


"""
    RCM <: EliminationAlgorithm

    RCM(; sortbydeg=true)

The [reverse Cuthill-McKee algorithm](https://en.wikipedia.org/wiki/Cuthill–McKee_algorithm). Uses [SymRCM.jl](https://github.com/PetrKryslUCSD/SymRCM.jl).
- `sortbydeg`: whether to sort neighbor lists by degree
"""
mutable struct RCM <: EliminationAlgorithm
    sortbydeg::Bool

    function RCM(; sortbydeg=true)
        new(sortbydeg)
    end
end


"""
    AMD <: EliminationAlgorithm

    AMD(; dense=nothing, aggressive=nothing)

The approximate minimum degree algorithm. Uses [AMD.jl](https://github.com/JuliaSmoothOptimizers/AMD.jl).
- `dense`: dense row parameter
- `aggressive`: aggressive absorbtion
"""
mutable struct AMD <: EliminationAlgorithm
    meta::AMDPkg.Amd

    function AMD(; dense=nothing, aggressive=nothing)
        meta = AMDPkg.Amd()

        if !isnothing(dense)
            meta.control[AMDPkg.AMD_DENSE] = dense
        end

        if !isnothing(aggressive)
            meta.control[AMDPkg.AMD_AGGRESSIVE] = aggressive
        end

        new(meta)
    end
end


"""
    SymAMD{Index} <: EliminationAlgorithm

    SymAMD{Index}(; dense_row=nothing, dense_col=nothing, aggressive=nothing) where Index

    SymAMD(; dense_row=nothing, dense_col=nothing, aggressive=nothing)

The column approximate minimum degree algorithm. Uses [AMD.jl](https://github.com/JuliaSmoothOptimizers/AMD.jl).
    - `Index`: either `Int` or `Cint`
    - `dense_row`: dense row parameter
    - `dense_column`: dense column parameter
    - `aggressive`: aggressive absorbtion
"""
mutable struct SymAMD{Index} <: EliminationAlgorithm
    meta::AMDPkg.Colamd{Index}

    function SymAMD(; dense_row=nothing, dense_col=nothing, aggressive=nothing)
        SymAMD{Int}(; dense_row, dense_col, aggressive)
    end

    function SymAMD{Index}(; dense_row=nothing, dense_col=nothing, aggressive=nothing) where Index
        meta = AMDPkg.Colamd{Index}()

        if !isnothing(dense_row)
            meta.knobs[AMDPkg.COLAMD_DENSE_ROW] = dense_row
        end

        if !isnothing(dense_col)
            meta.knobs[AMDPkg.COLAMD_DENSE_COL] = dense_col
        end

        if !isnothing(aggressive)
            meta.knobs[AMDPkg.COLAMD_AGGRESSIVE] = aggressive
        end

        new{Index}(meta)
    end
end


"""
    MMD <: EliminationAlgorithm

    MMD()

The multiple minimum degree algorithm. Uses [Sparspak.jl](https://github.com/PetrKryslUCSD/Sparspak.jl/tree/main).
"""
mutable struct MMD <: EliminationAlgorithm end



"""
    NodeND <: EliminationAlgorithm

    NodeND()

The [nested dissection heuristic](https://en.wikipedia.org/wiki/Nested_dissection). Uses [Metis.jl](https://github.com/JuliaSparse/Metis.jl).
"""
mutable struct NodeND <: EliminationAlgorithm end


"""
    FlowCutter <: EliminationAlgorithm

    FlowCutter(; time=10, seed=0)

The FlowCutter algorithm. Uses [FlowCutterPACE17_jll.jl](https://github.com/JuliaBinaryWrappers/FlowCutterPACE17_jll.jl). 
    - `time`: running time
    - `seed`: random seed
"""
mutable struct FlowCutter <: EliminationAlgorithm
    time::Int
    seed::Int
    history::Vector{String}

    function FlowCutter(; time=10, seed=0)
        new(time, seed, String[])
    end
end


"""
    Spectral <: EliminationAlgorithm

    Spectral(; tol=sqrt(eps(Float64)), restarts=200, mindim=nothing, maxdim=nothing)

Spectral ordering. Uses [ArnoldiMethod.jl](https://github.com/JuliaLinearAlgebra/ArnoldiMethod.jl).
    - `tol`: tolerance for convergence
    - `restarts`: maximum number of restarts
    - `mindim`: minimum Krylov dimension (≥ 2)
    - `maxdim`: maximum Krylov dimension (≥ mindim)
"""
mutable struct Spectral <: EliminationAlgorithm
    tol::Float64
    restarts::Int
    mindim::Union{Int, Nothing}
    maxdim::Union{Int, Nothing}
    history::Union{ArnoldiMethod.History, Nothing}

    function Spectral(; tol=sqrt(eps(Float64)), restarts=200, mindim=nothing, maxdim=nothing)
        new(tol, restarts, mindim, maxdim, nothing)
    end    
end


"""
    BT <: EliminationAlgorithm

    BT()

The Bouchitte-Todinca algorithm. Uses [TreeWidthSolver.jl](https://github.com/ArrogantGao/TreeWidthSolver.jl).
"""
mutable struct BT <: EliminationAlgorithm end



"""
    permutation(matrix::AbstractMatrix, alg::PermutationOrAlgorithm)

Construct a fill-reducing permutation of the vertices of a [simple graph](https://mathworld.wolfram.com/SimpleGraph.html).
- `matrix`: adjacency matrix
- `alg`: ordering algortihm
```julia
julia> using StructuredDecompositions.JunctionTrees

julia> graph = [
    0 1 1 0 0 0 0 0
    1 0 1 0 0 1 0 0
    1 1 0 1 1 0 0 0
    0 0 1 0 1 0 0 0
    0 0 1 1 0 0 1 1
    0 1 0 0 0 0 1 0
    0 0 0 0 1 1 0 1
    0 0 0 0 1 0 1 0
];

julia> order, index = permutation(graph, MCS());

julia> order
8-element Vector{Int64}:
 1
 6
 2
 3
 4
 5
 7
 8

julia> index == invperm(order)
true
```
"""
permutation(matrix::AbstractMatrix, alg::PermutationOrAlgorithm)


function permutation(matrix::AbstractMatrix, alg::EliminationAlgorithm)
    permutation(sparse(matrix), alg)
end


function permutation(matrix::AbstractMatrix, alg::AbstractVector)
    order = collect(alg)
    order, invperm(order)
end


function permutation(matrix::SparseMatrixCSC, alg::MCS)
    order = mcs(matrix)
    order, invperm(order)
end


function permutation(matrix::SparseMatrixCSC, alg::RCM)
    order = SymRCM.symrcm(matrix)
    order, invperm(order)
end


function permutation(matrix::SparseMatrixCSC, alg::AMD)
    order = AMDPkg.amd(matrix, alg.meta)
    order, invperm(order)
end


function permutation(matrix::SparseMatrixCSC, alg::SymAMD)
    order = AMDPkg.symamd(matrix, alg.meta)
    order, invperm(order)
end


function permutation(matrix::SparseMatrixCSC, alg::MMD)
    order = collect(axes(matrix, 2))
    index = collect(axes(matrix, 2))
    SpkMmd._generalmmd(size(matrix, 2), getcolptr(matrix), rowvals(matrix), order, index)
    order, index
end


function permutation(matrix::SparseMatrixCSC, alg::NodeND)
    order, index = Metis.permutation(matrix)
    order, index
end


function permutation(matrix::SparseMatrixCSC, alg::FlowCutter)
    nb, tw, nv, bagptr, bagval, treerow, treecol, history = flowcutter(matrix, alg.time, alg.seed)
    tree = Tree(sparse(treerow, treecol, ones(Bool, 2nb - 2), nb, nb), nb)
    order = sizehint!(Int[], nv)

    for i in tree
        bag = @view bagval[bagptr[i]:bagptr[i + 1] - 1]
        sort!(bag)
    end

    for i in invperm(dfs(tree))
        j = parentindex(tree, i)

        if isnothing(j)
            for v in @view bagval[bagptr[i]:bagptr[i + 1] - 1]
                push!(order, v)
            end
        else
            left = @view bagval[bagptr[i]:bagptr[i + 1] - 1]
            right = @view bagval[bagptr[j]:bagptr[j + 1] - 1]
            diffsorted!(order, left, right)
        end
    end

    append!(alg.history, history)
    order, invperm(order)
end


function permutation(matrix::SparseMatrixCSC, alg::Spectral)
    kwargs = Dict{Symbol, Real}(:tol => alg.tol, :restarts => alg.restarts)
       
    if !isnothing(alg.mindim)
        kwargs[:mindim] = alg.mindim        
    end

    if !isnothing(alg.maxdim)
        kwargs[:maxdim] = alg.maxdim
    end

    order, history = spectralorder(matrix; kwargs...)
    alg.history = history
    order, invperm(order)
end


function permutation(matrix::SparseMatrixCSC, alg::BT)
    T = TreeWidthSolver.LongLongUInt{size(matrix, 2) ÷ 64 + 1}
    fadjlist = Vector{Vector{Int}}(undef, size(matrix, 2))
    bitgraph = Vector{T}(undef, size(matrix, 2))

    for j in axes(matrix, 2)
        neighbors = rowvals(matrix)[nzrange(matrix, j)]
        fadjlist[j] = neighbors
        bitgraph[j] = TreeWidthSolver.bmask(T, neighbors)
    end

    graph = TreeWidthSolver.MaskedBitGraph(bitgraph, fadjlist, TreeWidthSolver.bmask(T, axes(matrix, 2)))
    cliques = TreeWidthSolver.all_pmc_enmu(graph, false)
    decomposition = TreeWidthSolver.bt_algorithm(graph, cliques, ones(size(matrix, 2)), false, true)
    order = reverse!(reduce(vcat, TreeWidthSolver.EliminationOrder(decomposition.tree).order))
    order, invperm(order)
end


function Base.show(io::IO, alg::RCM)
    output = "RCM:\n"
    output *= "    parameters:\n"
    output *= "        sortbydeg: $(alg.sortbydeg)\n"
    print(io, output)
end


function Base.show(io::IO, alg::AMD)
    meta = alg.meta
    output = "AMD:\n"
    output *= "    parameters:\n"
    output *= "        dense: $(meta.control[AMDPkg.AMD_DENSE])\n"
    output *= "        aggressive: $(meta.control[AMDPkg.AMD_AGGRESSIVE])\n"
    output *= "    information:\n"
    output *= "        status: $(AMDPkg.amd_statuses[meta.info[AMDPkg.AMD_STATUS]])\n"
    output *= "        matrix size: $(meta.info[AMDPkg.AMD_N])\n"
    output *= "        number of nonzeros: $(meta.info[AMDPkg.AMD_NZ])\n"
    output *= "        pattern symmetry: $(meta.info[AMDPkg.AMD_SYMMETRY])\n"
    output *= "        number of nonzeros on diagonal: $(meta.info[AMDPkg.AMD_NZDIAG])\n"
    output *= "        number of nonzeros in A + Aᵀ: $(meta.info[AMDPkg.AMD_NZ_A_PLUS_AT])\n"
    output *= "        number of dense columns: $(meta.info[AMDPkg.AMD_NDENSE])\n"
    output *= "        memory used: $(meta.info[AMDPkg.AMD_MEMORY])\n"
    output *= "        number of garbage collections: $(meta.info[AMDPkg.AMD_NCMPA])\n"
    output *= "        approximate number of nonzeros in factor: $(meta.info[AMDPkg.AMD_LNZ])\n"
    output *= "        number of float divides: $(meta.info[AMDPkg.AMD_NDIV])\n"
    output *= "        number of float * or - for LDL: $(meta.info[AMDPkg.AMD_NMULTSUBS_LDL])\n"
    output *= "        number of float * or - for LU: $(meta.info[AMDPkg.AMD_NMULTSUBS_LU])\n"
    output *= "        max nonzeros in any column of factor: $(meta.info[AMDPkg.AMD_DMAX])\n"
    print(io, output)
end


function Base.show(io::IO, alg::SymAMD)
    meta = alg.meta
    output = "SymAMD:\n"
    output *= "    parameters:\n"
    output *= "        dense row: $(meta.knobs[AMDPkg.COLAMD_DENSE_ROW])\n"
    output *= "        dense col: $(meta.knobs[AMDPkg.COLAMD_DENSE_COL])\n"
    output *= "        aggressive: $(meta.knobs[AMDPkg.COLAMD_AGGRESSIVE])\n"
    output *= "    information:\n"
    output *= "        status: $(AMDPkg.colamd_statuses[meta.stats[AMDPkg.COLAMD_STATUS]])\n"
    output *= "        memory defragmentation: $(meta.stats[AMDPkg.COLAMD_DEFRAG_COUNT])\n"
    print(io, output)
end


function Base.show(io::IO, alg::FlowCutter)
    output = "FlowCutter:\n"
    output *= "    parameters:\n"
    output *= "        time: $(alg.time)\n"
    output *= "        seed: $(alg.seed)\n"
    
    if !isempty(alg.history)
        output *= "    information:\n"
    
        for line in alg.history
            if startswith(line, "status")
                break
            end

            output *= "        $line\n"
        end
    end

    print(io, output)
end


function Base.show(io::IO, alg::Spectral)
    output = "Spectral:\n"
    output *= "    parameters:\n"
    output *= "        tol: $(alg.tol)\n"
    output *= "        restarts: $(alg.restarts)\n"

    if !isnothing(alg.mindim)
        output *= "        mindim: $(alg.mindim)\n"
    end

    if !isnothing(alg.maxdim)
        output *= "        maxdim: $(alg.maxdim)\n"
    end

    if !isnothing(alg.history)
        history = lowercase(string(alg.history))
        output *= "    information:\n"
        output *= "        $history\n"
    end

    print(io, output)
end


const DEFAULT_ELIMINATION_ALGORITHM = AMD()
