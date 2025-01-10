"""
    EliminationAlgorithm

A graph elimination algorithm. The options are

| type                 | name                              | complexity |
| :------------------- | :-------------------------------- | :--------- |
| [`RCM`](@ref)        | reverse Cuthill-Mckee             |            |
| [`AMD`](@ref)        | approximate minimum degree        |            |
| [`SymAMD`](@ref)     | column approximate minimum degree |            |
| [`NodeND`](@ref)     | nested dissection                 |            |
| [`BT`](@ref)         | Bouchitte-Todinca                 |            |
| [`MCS`](@ref)        | maximum cardinality search        | O(m + n)   |
| [`FlowCutter`](@ref) | FlowCutter                        |            |

for a graph with m vertices and n edges.
"""
abstract type EliminationAlgorithm end


"""
    PermutationOrAlgorithm = Union{AbstractVector, EliminationAlgorithm}

Either a permutation or a graph elimination algorithm.
"""
const PermutationOrAlgorithm = Union{AbstractVector, EliminationAlgorithm}


"""
    RCM <: EliminationAlgorithm

The [reverse Cuthill-McKee algorithm](https://en.wikipedia.org/wiki/Cuthill–McKee_algorithm). Uses [SymRCM.jl](https://github.com/PetrKryslUCSD/SymRCM.jl).
"""
mutable struct RCM <: EliminationAlgorithm
    sortbydeg::Bool
end


"""
    AMD <: EliminationAlgorithm

The approximate minimum degree algorithm. Uses [AMD.jl](https://github.com/JuliaSmoothOptimizers/AMD.jl).
"""
mutable struct AMD <: EliminationAlgorithm
    meta::AMDPkg.Amd
end

"""
    SymAMD{T} <: EliminationAlgorithm

The column approximate minimum degree algorithm. Uses [AMD.jl](https://github.com/JuliaSmoothOptimizers/AMD.jl).
"""
mutable struct SymAMD{T} <: EliminationAlgorithm
    meta::AMDPkg.Colamd{T}
end


"""
    NodeND <: EliminationAlgorithm

The [nested dissection heuristic](https://en.wikipedia.org/wiki/Nested_dissection). Uses [Metis.jl](https://github.com/JuliaSparse/Metis.jl).
"""
mutable struct NodeND <: EliminationAlgorithm end


"""
    BT <: EliminationAlgorithm

The Bouchitte-Todinca algorithm. Uses [TreeWidthSolver.jl](https://github.com/ArrogantGao/TreeWidthSolver.jl).
"""
mutable struct BT <: EliminationAlgorithm end


"""
    MCS <: EliminationAlgorithm

The maximum cardinality search algorithm.
"""
mutable struct MCS <: EliminationAlgorithm end


"""
    FlowCutter <: EliminationAlgorithm

The FlowCutter algorithm. Uses [FlowCutterPACE17_jll.jl](https://github.com/JuliaBinaryWrappers/FlowCutterPACE17_jll.jl). 
"""
mutable struct FlowCutter <: EliminationAlgorithm
    time::Int
    seed::Int
    comment::Vector{String}

    function FlowCutter(time, seed, comment)
        time < 0 && throw(ArgumentError("time < 0"))
        seed < 0 && throw(ArgumentError("seed < 0"))
        new(time, seed, comment)
    end
end


"""
    RCM(; sortbydeg=true)

Parameters for the RCM algorithm.
- `sortbydeg`: whether to sort neighbor lists by degree
"""
function RCM(; sortbydeg=true)
    RCM(sortbydeg)
end


"""
    AMD(; dense=nothing, aggressive=nothing)

Parameters for the AMD algorithm.
- `dense`: dense row parameter
- `aggressive`: aggressive absorbtion
"""
function AMD(; dense=nothing, aggressive=nothing)
    meta = AMDPkg.Amd()

    if !isnothing(dense)
        meta.control[AMDPkg.AMD_DENSE] = dense
    end

    if !isnothing(aggressive)
        meta.control[AMDPkg.AMD_AGGRESSIVE] = aggressive
    end

    AMD(meta)
end


function SymAMD{T}(; dense_row=nothing, dense_col=nothing, aggressive=nothing) where T
    meta = AMDPkg.Colamd{T}()

    if !isnothing(dense_row)
        meta.knobs[AMDPkg.COLAMD_DENSE_ROW] = dense_row
    end

    if !isnothing(dense_col)
        meta.knobs[AMDPkg.COLAMD_DENSE_COL] = dense_col
    end

    if !isnothing(aggressive)
        meta.knobs[AMDPkg.COLAMD_AGGRESSIVE] = aggressive
    end

    SymAMD(meta)
end


"""
    SymAMD(; dense_row=nothing, dense_col=nothing, aggressive=nothing)

Parameters for the SYMAMD algorithm.
- `dense_row`: dense row parameter
- `dense_column`: dense column parameter
- `aggressive`: aggressive absorbtion
"""
function SymAMD(; dense_row=nothing, dense_col=nothing, aggressive=nothing)
    SymAMD{Int}(; dense_row, dense_col, aggressive)
end


"""
    FlowCutter(; time=10, seed=0, verbose=false)

Parameters for the FlowCutter algorithm.
- `time`: running time
- `seed`: random seed
"""
function FlowCutter(; time=10, seed=0)
    FlowCutter(time, seed, String[])
end


"""
    ischordal(matrix::AbstractMatrix)

Determine whether a [simple graph](https://mathworld.wolfram.com/SimpleGraph.html), represented by its adjacency matrix `matrix`,
is [chordal](https://en.wikipedia.org/wiki/Chordal_graph).
"""
function ischordal(matrix::AbstractMatrix)
    ischordal(sparse(matrix))
end


function ischordal(matrix::SparseMatrixCSC)
    isperfect(matrix, permutation(matrix, MCS())...)
end


"""
    isperfect(matrix::AbstractMatrix, order::AbstractVector[, index::AbstractVector])

Determine whether an fill-reducing permutation is perfect.
"""
function isperfect(matrix::AbstractMatrix, order::AbstractVector, index::AbstractVector=invperm(order))
    isperfect(sparse(matrix))
end


# Simple Linear-Time Algorithms to Test Chordality of Graphs, Test Acyclicity of Hypergraphs, and Selectively Reduce Acyclic Hypergraphs
# Tarjan and Yannakakis
# Test for Zero Fill-In.
#
# Determine whether a fill-reducing permutation is perfect.
function isperfect(matrix::SparseMatrixCSC, order::AbstractVector, index::AbstractVector=invperm(order))
    f = Vector{Int}(undef, size(matrix, 2))
    findex = similar(f)
    
    for (i, w) in enumerate(order)
        f[w] = w
        findex[w] = i
        
        for v in @view rowvals(matrix)[nzrange(matrix, w)]
            if index[v] < i
                findex[v] = i
                
                if f[v] == v
                    f[v] = w
                end
            end
        end
        
        for v in @view rowvals(matrix)[nzrange(matrix, w)]
            if index[v] < i && findex[f[v]] < i
                return false
            end
        end
    end
    
    true
end


"""
    permutation(matrix::AbstractMatrix, alg::PermutationOrAlgorithm)

Construct a fill-reducing permutation of the vertices of a [simple graph](https://mathworld.wolfram.com/SimpleGraph.html),
represented by its adjacency matrix `matrix`.

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


function permutation(matrix::SparseMatrixCSC, alg::NodeND)
    order, index = Metis.permutation(matrix)
    order, index
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


function permutation(matrix::SparseMatrixCSC, alg::MCS)
    order = mcs(matrix)
    order, invperm(order)
end


function permutation(matrix::SparseMatrixCSC, alg::FlowCutter)
    comment, nb, tw, nv, bagptr, bagval, tree = flowcutter(matrix, alg.time, alg.seed)
    alg.comment = comment
    tree = Tree(tree, nb)
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

    order, invperm(order)
end


function Base.show(io::IO, alg::EliminationAlgorithm)
    string = "Parameters:\n"
    print(io, string)
end


function Base.show(io::IO, alg::RCM)
    string = "Parameters:\n"
    string *= "  sortbydeg: $(alg.sortbydeg)\n"
    print(io, string)
end


function Base.show(io::IO, alg::AMD)
    print(io, alg.meta)
end


function Base.show(io::IO, alg::SymAMD)
    print(io, alg.meta)
end


function Base.show(io::IO, alg::FlowCutter)
    string = "Parameters:\n"
    string *= "  time: $(alg.time)\n"
    string *= "  seed: $(alg.seed)\n"
    string *= "Information:\n"

    for c in alg.comment
        string *= "  $c\n"
    end

    print(io, string)
end


const DEFAULT_ELIMINATION_ALGORITHM = AMD()
