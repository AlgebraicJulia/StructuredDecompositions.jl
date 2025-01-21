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
    isperfect(matrix, permutation(matrix, MCS())...)
end


"""
    isfilled(matrix::AbstractMatrix)

Determine whether a directed graph is filled.
"""
function isfilled(matrix::AbstractMatrix)
    isfilled(sparse(matrix))
end


# Determine whether a directed graph is filled.
function isfilled(matrix::SparseMatrixCSC)
    isperfect(matrix, axes(matrix)...)
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


# Construct the intersection graph of a junction tree.
function filledgraph(tree::JunctionTree)
    filledgraph(true, tree)
end


"""
    filledgraph([element=true,] tree::JunctionTree)

See below. The function returns a sparse matrix whose structural nonzeros are filled with `element`.
"""
function filledgraph(element::T, tree::JunctionTree) where T
    matrix = filledgraph(T, tree)
    fill!(nonzeros(matrix), element)
    matrix
end


"""
    filledgraph(T::Type, tree::JunctionTree)

See below. The function returns a sparse matrix with elements of type `T`.
"""
function filledgraph(T::Type, tree::JunctionTree)
    n = last(residual(last(tree)))
    filledgraph!(spzeros(T, n, n), tree)
end


"""
    filledgraph!(target::SparseMatrixCSC, tree::JunctionTree)

Construct the [intersection graph](https://en.wikipedia.org/wiki/Intersection_graph) of the subtrees of
a junction tree. The result is stored in `target`.
"""
function filledgraph!(target::SparseMatrixCSC, tree::JunctionTree)
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
