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

The [reverse Cuthill-McKee algorithm](https://en.wikipedia.org/wiki/Cuthill–McKee_algorithm).
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
In order to use this algorithm, import the package [Metis](https://github.com/JuliaSparse/Metis.jl).
"""
struct NodeND <: EliminationAlgorithm end


"""
    Spectral <: EliminationAlgorithm

    Spectral(; tol=0.0)

The spectral ordering algorithm.
In order to use this algorithm, import the package [Laplacians](https://github.com/danspielman/Laplacians.jl).
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
In order to use this algorithm, import the package [TreeWidthSolver](https://github.com/ArrogantGao/TreeWidthSolver.jl).
"""
struct BT <: EliminationAlgorithm end



"""
    permutation(matrix::AbstractMatrix;
        alg::PermutationOrAlgorithm=DEFAULT_ELIMINATION_ALGORITHM)

Construct a fill-reducing permutation of the vertices of a simple graph.
```julia
julia> using StructuredDecompositions

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

julia> order, index = permutation(graph; alg=MCS());

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
function permutation(matrix::AbstractMatrix; alg::PermutationOrAlgorithm=DEFAULT_ELIMINATION_ALGORITHM)
    permutation(matrix, alg)
end


function permutation(matrix::AbstractMatrix, alg::EliminationAlgorithm)
    permutation(sparse(matrix), alg)
end


function permutation(matrix::SparseMatrixCSC, alg::T) where T <: Union{NodeND, Spectral, BT}
    throwextension(T)
end


function permutation(matrix::AbstractMatrix, alg::AbstractVector)
    order = Vector{Int}(alg)
    order, invperm(order)
end


function permutation(matrix::SparseMatrixCSC, alg::MCS)
    index, size = mcs(matrix)
    invperm(index), index
end


function permutation(matrix::SparseMatrixCSC, alg::RCM)
    order = SymRCM.symrcm(matrix)
    order, invperm(order)
end


function permutation(matrix::SparseMatrixCSC, alg::AMD)
    meta = AMDPkg.Amd()
    meta.control[AMDPkg.AMD_DENSE] = alg.dense
    meta.control[AMDPkg.AMD_AGGRESSIVE] = alg.aggressive
    order = AMDPkg.amd(matrix, meta)
    order, invperm(order)
end


function permutation(matrix::SparseMatrixCSC{T, I}, alg::SymAMD) where {T, I}
    meta = AMDPkg.Colamd{I}()
    meta.knobs[AMDPkg.COLAMD_DENSE_ROW] = alg.dense_row
    meta.knobs[AMDPkg.COLAMD_DENSE_COL] = alg.dense_col
    meta.knobs[AMDPkg.COLAMD_AGGRESSIVE] = alg.aggressive
    order = AMDPkg.symamd(matrix, meta)
    order, invperm(order)
end


function permutation(matrix::SparseMatrixCSC, alg::MMD)
    order = collect(axes(matrix, 2))
    index = collect(axes(matrix, 2))
    SpkMmd._generalmmd(size(matrix, 2), getcolptr(matrix), rowvals(matrix), order, index)
    order, index
end


function throwextension(::Type{NodeND})
    throw(ArgumentError("In order to use the algorithm `NodeND`, import the package `Metis`."))
end


function throwextension(::Type{Spectral})
    throw(ArgumentError("In order to use the algorithm `Spectral`, import the package `Laplacians`."))
end


function throwextension(::Type{BT})
    throw(ArgumentError("In order to use the algorithm `BT`, import the package `TreeWidthSolver`."))
end


function Base.show(io::IO, alg::RCM)
    println(io, "RCM:")
    println(io, "   sortbydeg: $(alg.sortbydeg)")
end


function Base.show(io::IO, alg::AMD)
    println(io, "AMD:")
    println(io, "   dense: $(alg.dense)")
    println(io, "   aggressive: $(alg.aggressive)")
end


function Base.show(io::IO, alg::SymAMD)
    println(io, "SymAMD:")
    println(io, "   dense row: $(alg.dense_row)")
    println(io, "   dense col: $(alg.dense_col)")
    println(io, "   aggressive: $(alg.aggressive)")
end


function Base.show(io::IO, alg::Spectral)
    println(io, "Spectral:")
    println(io, "   tol: $(alg.tol)")
end


"""
    DEFAULT_ELIMINATION_ALGORITHM = AMD()

The default algorithm.
"""
const DEFAULT_ELIMINATION_ALGORITHM = AMD()
