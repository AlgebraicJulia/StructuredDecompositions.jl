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

Either a permutation or an algorithm.
"""
const PermutationOrAlgorithm = Union{AbstractVector, EliminationAlgorithm}


"""
    MCS <: EliminationAlgorithm

    MCS()

The maximum cardinality search algorithm.
"""
struct MCS <: EliminationAlgorithm end


"""
    RCM <: EliminationAlgorithm

    RCM(; sortbydeg=true)

The [reverse Cuthill-McKee algorithm](https://en.wikipedia.org/wiki/Cuthill–McKee_algorithm). Uses [SymRCM.jl](https://github.com/PetrKryslUCSD/SymRCM.jl).
- `sortbydeg`: whether to sort neighbor lists by degree
"""
struct RCM <: EliminationAlgorithm
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
struct AMD <: EliminationAlgorithm
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
struct SymAMD{Index} <: EliminationAlgorithm
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

The [multiple minimum degree algorithm](https://en.wikipedia.org/wiki/Minimum_degree_algorithm).
Uses [Sparspak.jl](https://github.com/PetrKryslUCSD/Sparspak.jl/tree/main).
"""
struct MMD <: EliminationAlgorithm end



"""
    NodeND <: EliminationAlgorithm

    NodeND()

The [nested dissection algorithm](https://en.wikipedia.org/wiki/Nested_dissection). Uses [Metis.jl](https://github.com/JuliaSparse/Metis.jl).
"""
struct NodeND <: EliminationAlgorithm end


"""
    FlowCutter <: EliminationAlgorithm

    FlowCutter(; time=10, seed=0)

The FlowCutter algorithm. Uses [FlowCutterPACE17_jll.jl](https://github.com/JuliaBinaryWrappers/FlowCutterPACE17_jll.jl). 
- `time`: running time
- `seed`: random seed
"""
struct FlowCutter <: EliminationAlgorithm
    time::Int
    seed::Int
    history::Vector{String}

    function FlowCutter(; time=10, seed=0)
        new(time, seed, String[])
    end
end


"""
    Spectral <: EliminationAlgorithm

    Spectral(; tol=0.0)

The spectral ordering algorithm. Uses [Laplacians.jl](https://github.com/danspielman/Laplacians.jl).
- `tol`: tolerance for convergence
"""
struct Spectral <: EliminationAlgorithm
    tol::Float64

    function Spectral(; tol=0.0)
        new(tol)
    end    
end


"""
    BT <: EliminationAlgorithm

    BT()

The Bouchitte-Todinca algorithm. Uses [TreeWidthSolver.jl](https://github.com/ArrogantGao/TreeWidthSolver.jl).
"""
struct BT <: EliminationAlgorithm end



"""
    permutation(matrix::AbstractMatrix, alg::PermutationOrAlgorithm)

Construct a fill-reducing permutation of the vertices of a simple graph.
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
    order = spectralorder(matrix; tol=alg.tol)
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
    println(io, "RCM:")
    println(io, "    parameters:")
    println(io, "        sortbydeg: $(alg.sortbydeg)")
end


function Base.show(io::IO, alg::AMD)
    meta = alg.meta
    println(io, "AMD:")
    println(io, "    parameters:")
    println(io, "        dense: $(meta.control[AMDPkg.AMD_DENSE])")
    println(io, "        aggressive: $(meta.control[AMDPkg.AMD_AGGRESSIVE])")
    println(io, "    information:")
    println(io, "        status: $(AMDPkg.amd_statuses[meta.info[AMDPkg.AMD_STATUS]])")
    println(io, "        matrix size: $(meta.info[AMDPkg.AMD_N])")
    println(io, "        number of nonzeros: $(meta.info[AMDPkg.AMD_NZ])")
    println(io, "        pattern symmetry: $(meta.info[AMDPkg.AMD_SYMMETRY])")
    println(io, "        number of nonzeros on diagonal: $(meta.info[AMDPkg.AMD_NZDIAG])")
    println(io, "        number of nonzeros in A + Aᵀ: $(meta.info[AMDPkg.AMD_NZ_A_PLUS_AT])")
    println(io, "        number of dense columns: $(meta.info[AMDPkg.AMD_NDENSE])")
    println(io, "        memory used: $(meta.info[AMDPkg.AMD_MEMORY])")
    println(io, "        number of garbage collections: $(meta.info[AMDPkg.AMD_NCMPA])")
    println(io, "        approximate number of nonzeros in factor: $(meta.info[AMDPkg.AMD_LNZ])")
    println(io, "        number of float divides: $(meta.info[AMDPkg.AMD_NDIV])")
    println(io, "        number of float * or - for LDL: $(meta.info[AMDPkg.AMD_NMULTSUBS_LDL])")
    println(io, "        number of float * or - for LU: $(meta.info[AMDPkg.AMD_NMULTSUBS_LU])")
    println(io, "        max nonzeros in any column of factor: $(meta.info[AMDPkg.AMD_DMAX])")
end


function Base.show(io::IO, alg::SymAMD)
    meta = alg.meta
    println(io, "SymAMD:")
    println(io, "    parameters:")
    println(io, "        dense row: $(meta.knobs[AMDPkg.COLAMD_DENSE_ROW])")
    println(io, "        dense col: $(meta.knobs[AMDPkg.COLAMD_DENSE_COL])")
    println(io, "        aggressive: $(meta.knobs[AMDPkg.COLAMD_AGGRESSIVE])")
    println(io, "    information:")
    println(io, "        status: $(AMDPkg.colamd_statuses[meta.stats[AMDPkg.COLAMD_STATUS]])")
    println(io, "        memory defragmentation: $(meta.stats[AMDPkg.COLAMD_DEFRAG_COUNT])")
end


function Base.show(io::IO, alg::FlowCutter)
    println(io, "FlowCutter:")
    println(io, "    parameters:")
    println(io, "        time: $(alg.time)")
    println(io, "        seed: $(alg.seed)")
    
    if !isempty(alg.history)
        println(io, "    information:")
    
        for line in alg.history
            if startswith(line, "status")
                break
            end

            println(io, "        $line")
        end
    end
end


function Base.show(io::IO, alg::Spectral)
    println(io, "Spectral:")
    println(io, "    parameters:")
    println(io, "        tol: $(alg.tol)")
end


const DEFAULT_ELIMINATION_ALGORITHM = AMD()
