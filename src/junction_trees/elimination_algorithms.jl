"""
    EliminationAlgorithm

A graph elimination algorithm. The options are

| type                 | name                              | complexity  | connected |
| :------------------- | :-------------------------------- | :---------- | :-------- |
| [`MCS`](@ref)        | maximum cardinality search        | O(m + n)    | no        |
| [`RCM`](@ref)        | reverse Cuthill-Mckee             | O(mΔ)       | yes       |
| [`AMD`](@ref)        | approximate minimum degree        | O(mn)       | no        |
| [`SymAMD`](@ref)     | column approximate minimum degree |             | no        |
| [`MMD`](@ref)        | multiple minimum degree           | O(mn²)      | no        |
| [`NodeND`](@ref)     | nested dissection                 |             | no        |
| [`Spectral`](@ref)   | spectral ordering                 |             | yes       |
| [`BT`](@ref)         | Bouchitte-Todinca                 | O(2.6183ⁿ)  | no        |

for a graph with m edges, n vertices, and maximum degree Δ.
The algorithms [`RCM`](@ref) and [`Spectral`](@ref) only work on connected graphs.
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

The [reverse Cuthill-McKee algorithm](https://en.wikipedia.org/wiki/Cuthill–McKee_algorithm)
only works on connected graphs.
"""
struct RCM <: EliminationAlgorithm end


"""
    AMD <: EliminationAlgorithm

    AMD(; dense=10.0, aggressive=1.0)

The approximate minimum degree algorithm.
- `dense`: dense row parameter
- `aggressive`: aggressive absorbtion
"""
struct AMD <: EliminationAlgorithm
    dense::Float64
    aggressive::Float64

    function AMD(; dense=10.0, aggressive=1.0)
        new(dense, aggressive)
    end
end


"""
    SymAMD <: EliminationAlgorithm

    SymAMD(; dense_row=10.0, dense_col=10.0, aggressive=1.0)

The column approximate minimum degree algorithm.
- `dense_row`: dense row parameter
- `dense_column`: dense column parameter
- `aggressive`: aggressive absorbtion
"""
struct SymAMD <: EliminationAlgorithm
    dense_row::Float64
    dense_col::Float64
    aggressive::Float64

    function SymAMD(; dense_row=10.0, dense_col=10.0, aggressive=1.0)
        new(dense_row, dense_col, aggressive)
    end
end


"""
    MMD <: EliminationAlgorithm

    MMD()

The [multiple minimum degree algorithm](https://en.wikipedia.org/wiki/Minimum_degree_algorithm).
"""
struct MMD <: EliminationAlgorithm end



"""
    NodeND <: EliminationAlgorithm

    NodeND()

The [nested dissection algorithm](https://en.wikipedia.org/wiki/Nested_dissection).
In order to use it, import the package [Metis](https://github.com/JuliaSparse/Metis.jl).
"""
struct NodeND <: EliminationAlgorithm end


"""
    Spectral <: EliminationAlgorithm

    Spectral(; tol=0.0)

The spectral ordering algorithm only works on connected graphs.
In order to use it, import the package [Laplacians](https://github.com/danspielman/Laplacians.jl).
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

The Bouchitte-Todinca algorithm.
In order to use it, import the package [TreeWidthSolver](https://github.com/ArrogantGao/TreeWidthSolver.jl).
"""
struct BT <: EliminationAlgorithm end



"""
    permutation(graph;
        alg::PermutationOrAlgorithm=DEFAULT_ELIMINATION_ALGORITHM)

Construct a fill-reducing permutation of the vertices of a simple graph.
```julia
julia> using SparseArrays, StructuredDecompositions

julia> graph = sparse([
           0 1 1 0 0 0 0 0
           1 0 1 0 0 1 0 0
           1 1 0 1 1 0 0 0
           0 0 1 0 1 0 0 0
           0 0 1 1 0 0 1 1
           0 1 0 0 0 0 1 0
           0 0 0 0 1 1 0 1
           0 0 0 0 1 0 1 0
       ]);

julia> order, index = permutation(graph);

julia> order
8-element Vector{Int64}:
 4
 8
 7
 6
 5
 1
 3
 2

julia> index == invperm(order)
true
```
"""
function permutation(graph; alg::PermutationOrAlgorithm=DEFAULT_ELIMINATION_ALGORITHM)
    permutation(graph, alg)
end


function permutation(graph, alg::PermutationOrAlgorithm)
    permutation(Graph(graph), alg)
end


function permutation(graph::Graph{V}, alg::AbstractVector) where V
    order::Vector{V} = alg
    order, invperm(order)
end


function permutation(graph::SparseMatrixCSC{<:Any, I}, alg::AbstractVector) where I
    order::Vector{I} = alg
    order, invperm(order)
end


function permutation(graph::Graph, alg::MCS)
    index = mcs(graph)
    invperm(index), index
end


function permutation(graph::Graph, alg::RCM)
    order = rcm(graph)
    order, invperm(order)
end


function permutation(graph::Graph, alg::Union{AMD, SymAMD, MMD, NodeND, Spectral, BT})
    matrix = SparseMatrixCSC{Bool}(graph)
    permutation(matrix, alg)
end


function permutation(matrix::SparseMatrixCSC{<:Any, I}, alg::AMD) where I
    meta = AMDPkg.Amd()
    meta.control[AMDPkg.AMD_DENSE] = alg.dense
    meta.control[AMDPkg.AMD_AGGRESSIVE] = alg.aggressive
    order::Vector{I} = AMDPkg.amd(matrix, meta)
    order, invperm(order)
end


function permutation(matrix::SparseMatrixCSC{<:Any, I}, alg::SymAMD) where I
    meta = AMDPkg.Colamd{I}()
    meta.knobs[AMDPkg.COLAMD_DENSE_ROW] = alg.dense_row
    meta.knobs[AMDPkg.COLAMD_DENSE_COL] = alg.dense_col
    meta.knobs[AMDPkg.COLAMD_AGGRESSIVE] = alg.aggressive
    order::Vector{I} = AMDPkg.symamd(matrix, meta)
    order, invperm(order)
end


function permutation(matrix::SparseMatrixCSC{<:Any, I}, alg::MMD) where I
    order::Vector{I} = axes(matrix, 2)
    index::Vector{I} = axes(matrix, 2)
    SpkMmd._generalmmd(size(matrix, 2), getcolptr(matrix), rowvals(matrix), order, index)
    order, index
end


function Base.show(io::IO, ::MIME"text/plain", alg::AMD)
    println(io, "AMD:")
    println(io, "   dense: $(alg.dense)")
    println(io, "   aggressive: $(alg.aggressive)")
end


function Base.show(io::IO, ::MIME"text/plain", alg::SymAMD)
    println(io, "SymAMD:")
    println(io, "   dense row: $(alg.dense_row)")
    println(io, "   dense col: $(alg.dense_col)")
    println(io, "   aggressive: $(alg.aggressive)")
end


function Base.show(io::IO, ::MIME"text/plain", alg::Spectral)
    println(io, "Spectral:")
    println(io, "   tol: $(alg.tol)")
end


"""
    DEFAULT_ELIMINATION_ALGORITHM = AMD()

The default algorithm.
"""
const DEFAULT_ELIMINATION_ALGORITHM = AMD()
