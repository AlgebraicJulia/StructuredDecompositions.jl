# Construct the subtree graph of a junction tree.
function eliminationgraph(tree::JunctionTree)
    eliminationgraph(true, tree)
end


"""
    eliminationgraph([element=true,] tree::JunctionTree)

See [`eliminationgraph!`]. The function returns a sparse matrix whose structural nonzeros are filled with `element`.
"""
function eliminationgraph(element::T, tree::JunctionTree) where T
    matrix = eliminationgraph(T, tree)
    fill!(nonzeros(matrix), element)
    matrix
end


"""
    eliminationgraph(T::Type, tree::JunctionTree)

See [`eliminationgraph!`]. The function returns a sparse matrix with elements of type `T`.
"""
function eliminationgraph(T::Type, tree::JunctionTree)
    n = last(residual(last(tree)))
    eliminationgraph!(spzeros(T, n, n), tree)
end


"""
    eliminationgraph!(target::SparseMatrixCSC, tree::JunctionTree)

Construct the [subtree graph](https://en.wikipedia.org/wiki/Chordal_graph) of 
a junction tree. The result is stored in `target`.
"""
function eliminationgraph!(target::SparseMatrixCSC, tree::JunctionTree)
    sizehint!(empty!(rowvals(target)), nnz(tree))
    push!(empty!(getcolptr(target)), 1)

    for bag in tree
        res = residual(bag)
        sep = separator(bag)

        for i in eachindex(res)
            append!(rowvals(target), res[i + 1:end])
            append!(rowvals(target), sep)
            push!(getcolptr(target), length(rowvals(target)) + 1)
        end
    end

    resize!(nonzeros(target), nnz(tree))
    target
end


"""
    eliminationgraph(matrix::AbstractMatrix;
        alg::PermutationOrAlgorithm=DEFAULT_ELIMINATION_ALGORITHM)

Construct the elimination graph of a simple graph.
```julia
julia> using StructuredDecompositions, SparseArrays

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

julia> label, filled = label, filled = eliminationgraph(graph);

julia> filled
8×8 SparseMatrixCSC{Int64, Int64} with 13 stored entries:
 ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅
 ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅
 ⋅  1  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅
 ⋅  ⋅  1  ⋅  ⋅  ⋅  ⋅  ⋅
 ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅
 1  1  1  1  ⋅  ⋅  ⋅  ⋅
 1  ⋅  ⋅  ⋅  1  1  ⋅  ⋅
 ⋅  ⋅  ⋅  1  1  1  1  ⋅

julia> isfilled(filled)
true

julia> ischordal(filled + filled')
true
"""
function eliminationgraph(matrix::AbstractMatrix; alg::PermutationOrAlgorithm=DEFAULT_ELIMINATION_ALGORITHM)
    eliminationgraph!(sparse(matrix); alg)
end


"""
    eliminationgraph!(matrix::SparseMatrixCSC;
        alg::PermutationOrAlgorithm=DEFAULT_ELIMINATION_ALGORITHM)

A mutating version of [`eliminationgraph`](@ref).
"""
function eliminationgraph!(matrix::SparseMatrixCSC; alg::PermutationOrAlgorithm=DEFAULT_ELIMINATION_ALGORITHM)
    label, tree = junctiontree!(matrix; alg)
    label, eliminationgraph!(matrix, tree)
end


"""
    nnz(tree::JunctionTree)

Compute the number of edges in the [intersection graph](https://en.wikipedia.org/wiki/Intersection_graph) of a junction tree.
"""
function SparseArrays.nnz(tree::JunctionTree)
    sum(tree; init=0) do bag
        m = length(residual(bag))
        n = length(separator(bag))
        (n)m + (m - 1)m ÷ 2
    end 
end


"""
    ischordal(matrix::AbstractMatrix)

Determine whether a simple graph is [chordal](https://en.wikipedia.org/wiki/Chordal_graph).
"""
function ischordal(matrix::AbstractMatrix)
    ischordal(sparse(matrix))
end


# Determine whether a graph is chordal.
function ischordal(matrix::SparseMatrixCSC)
    # validate arguments
    size(matrix, 1) != size(matrix, 2) && throw(ArgumentError("size(matrix, 1) != size(matrix, 2)"))

    # run algorithm
    ischordal(size(matrix, 2)) do j
        @view rowvals(matrix)[nzrange(matrix, j)]
    end
end


function ischordal(neighbors::Function, n::Integer)
    # validate arguments
    n < 0 && throw(ArgumentError("n < 0"))

    # run algorithm
    index, size = mcs(neighbors, n)
    isperfect(neighbors, invperm(index), index)
end


"""
    isfilled(matrix::AbstractMatrix)

Determine whether a directed graph is filled.
"""
function isfilled(matrix::AbstractMatrix)
    isfilled(sparse(matrix))
end


function isfilled(matrix::SparseMatrixCSC)
    # validate arguments
    size(matrix, 1) != size(matrix, 2) && throw(ArgumentError("size(matrix, 1) != size(matrix, 2)"))

    # run algorithm
    isfilled(size(matrix, 2)) do j
        @view rowvals(matrix)[nzrange(matrix, j)]
    end
end


function isfilled(neighbors::Function, n::Integer)
    # validate arguments
    n < 0 && throw(ArgumentError("n < 0"))

    # run algorithm
    isperfect(neighbors, 1:n, 1:n)
end


"""
    isperfect(matrix::AbstractMatrix, order::AbstractVector[, index::AbstractVector])

Determine whether an fill-reducing permutation is perfect.
"""
function isperfect(matrix::AbstractMatrix, order::AbstractVector, index::AbstractVector=invperm(order))
    isperfect(sparse(matrix))
end


function isperfect(matrix::SparseMatrixCSC, order::AbstractVector, index::AbstractVector=invperm(order))
    # validate arguments
    size(matrix, 1) != size(matrix, 2) && throw(ArgumentError("size(matrix, 1) != size(matrix, 2)"))
    axes(matrix, 2) != eachindex(order) && throw(ArgumentError("axes(matrix, 2) != eachindex(order)"))

    # run algorithm
    isperfect(order, index) do j
        @view rowvals(matrix)[nzrange(matrix, j)]
    end
end


# Simple Linear-Time Algorithms to Test Chordality of Graphs, Test Acyclicity of Hypergraphs, and Selectively Reduce Acyclic Hypergraphs
# Tarjan and Yannakakis
# Test for Zero Fill-In.
#
# Determine whether a fill-reducing permutation is perfect.
function isperfect(neighbors::Function, order::AbstractVector, index::AbstractVector=invperm(order))
    # validate arguments
    eachindex(order) != eachindex(index) && throw(ArgumentError("eachindex(order) != eachindex(index)"))

    # run algorithm
    f = Vector{Int}(undef, length(order))
    findex = Vector{Int}(undef, length(order))

    for (i, w) in enumerate(order)
        f[w] = w
        findex[w] = i

        for v in neighbors(w)
            if index[v] < i
                findex[v] = i

                if f[v] == v
                    f[v] = w
                end
            end
        end

        for v in neighbors(w)
            if index[v] < i && findex[f[v]] < i
                return false
            end
        end
    end

    true
end
