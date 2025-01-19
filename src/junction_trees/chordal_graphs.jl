"""
    nnz(tree::JunctionTree)

Compute the number of edges in the [intersection graph](https://en.wikipedia.org/wiki/Intersection_graph) of a junction tree.
"""
function SparseArrays.nnz(tree::JunctionTree)
    sum(tree) do bag
        r = length(residual(bag))
        s = length(separator(bag))
        div((2s + r - 1)r, 2)
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


# Construct the intersection graph of a junction tree.
function chordalgraph(tree::JunctionTree)
    chordalgraph(true, tree)
end


# Construct a chordal completion of a connected simple graph.
function chordalgraph(matrix::AbstractMatrix; alg::PermutationOrAlgorithm=DEFAULT_ELIMINATION_ALGORITHM)
    label, tree = junctiontree(matrix; alg)
    chordalgraph(tree)
end


"""
    chordalgraph([element=true,] tree::JunctionTree)

See below. The function returns a sparse matrix whose structural nonzeros are filled with `element`.
"""
function chordalgraph(element::T, tree::JunctionTree) where T
    matrix = chordalgraph(T, tree)
    fill!(nonzeros(matrix), element)
    matrix
end


"""
    chordalgraph([element=true,] matrix::AbstractMatrix)

See below. The function returns a sparse matrix whose structural nonzeros are filled with `element`.
"""
function chordalgraph(element, matrix::AbstractMatrix; alg::PermutationOrAlgorithm=DEFAULT_ELIMINATION_ALGORITHM)
    label, tree = junctiontree(matrix; alg)
    label, chordalgraph(element, tree)
end


"""
    chordalgraph(T::Type, tree::JunctionTree)

Construct the [intersection graph](https://en.wikipedia.org/wiki/Intersection_graph) of the subtrees of
a junction tree. The function returns a sparse matrix with elements of type `T`
and the same sparsity structure as the lower triangular part of the graph's adjacency matrix.
"""
function chordalgraph(T::Type, tree::JunctionTree)
    n = last(residual(last(tree)))
    colptr = sizehint!(Int[], n + 1)
    rowval = sizehint!(Int[], nnz(tree))
    push!(colptr, 1) 

    for bag in tree
        res = residual(bag)
        sep = separator(bag)

        for i in eachindex(res)
            append!(rowval, res[i + 1:end])
            append!(rowval, sep)
            push!(colptr, length(rowval) + 1)
        end
    end

    nzval = Vector{T}(undef, nnz(tree))
    SparseMatrixCSC(n, n, colptr, rowval, nzval)
end


"""
    chordalgraph(T::Type, matrix::AbstractMatrix)

Construct the chordal completion of a connected simple graph. 
The function returns a sparse matrix with elements of type `T`
and the same sparsity structure as the lower triangular part of the graph's adjacency matrix.
"""
chordalgraph(T::Type, matrix::AbstractMatrix)
