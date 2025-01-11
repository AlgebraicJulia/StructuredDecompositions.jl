"""
    nnz(tree::JunctionTree)

Compute the number of edges in the [intersection graph](https://en.wikipedia.org/wiki/Intersection_graph) of a junction tree.
"""
function SparseArrays.nnz(tree::JunctionTree)
    sum(tree) do bag
        r = length(residual(bag))
        s = length(separator(bag))
        (2s + r - 1) * r ÷ 2
    end 
end


"""
    ischordal(matrix::AbstractMatrix)

Determine whether a [simple graph](https://mathworld.wolfram.com/SimpleGraph.html), represented by its adjacency matrix `matrix`,
is [chordal](https://en.wikipedia.org/wiki/Chordal_graph).
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
    chordalgraph(Element::Type, tree::JunctionTree)

Construct the [intersection graph](https://en.wikipedia.org/wiki/Intersection_graph) of the subtrees of
a junction tree. The function returns a sparse matrix with elements of type `Element`
and the same sparsity structure as the lower triangular part of the graph's adjacency matrix.
"""
function chordalgraph(Element::Type, tree::JunctionTree)
    colptr = Vector{Int}(undef, last(tree.sndptr))
    rowval = Vector{Int}(undef, nnz(tree))
    nzval = Vector{Element}(undef, nnz(tree))
    j = p = colptr[1] = 1

    for bag in tree
        for i in eachindex(residual(bag))
            for v in @view bag[i + 1:end]
                rowval[p] = v
                p += 1
            end

            colptr[j + 1] = p
            j += 1
        end
    end

    SparseMatrixCSC(last(tree.sndptr) - 1, last(tree.sndptr) - 1, colptr, rowval, nzval)
end
