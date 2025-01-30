"""
    eliminationgraph([element=true,] tree::JunctionTree)

See [`eliminationgraph!`](@ref). The function returns a sparse matrix whose structural nonzeros are filled with `element`.
"""
function eliminationgraph(tree::JunctionTree)
    eliminationgraph(true, tree)
end


function eliminationgraph(element::T, tree::JunctionTree) where T
    matrix = eliminationgraph(T, tree)
    fill!(nonzeros(matrix), element)
    matrix
end


function eliminationgraph(T::Type, tree::JunctionTree{I}) where I
    n = residual(tree[end])[end]
    eliminationgraph!(spzeros(T, I, n, n), tree)
end


"""
    eliminationgraph!(graph, tree::JunctionTree)

Construct the [subtree graph](https://en.wikipedia.org/wiki/Chordal_graph) of 
a junction tree. The result is stored in `graph`.
"""
function eliminationgraph!(matrix::SparseMatrixCSC, tree::JunctionTree)
    sizehint!(empty!(rowvals(matrix)), nnz(tree))
    push!(empty!(getcolptr(matrix)), 1)

    for bag in tree
        res = residual(bag)
        sep = separator(bag)

        for i in eachindex(res)
            append!(rowvals(matrix), res[i + 1:end])
            append!(rowvals(matrix), sep)
            push!(getcolptr(matrix), length(rowvals(matrix)) + 1)
        end
    end

    resize!(nonzeros(matrix), nnz(tree))
    matrix
end


"""
    eliminationgraph([element=true,] graph;
        alg::PermutationOrAlgorithm=DEFAULT_ELIMINATION_ALGORITHM,
        snd::SupernodeType=DEFAULT_SUPERNODE_TYPE)

Construct the elimination graph of a simple graph.
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

julia> label, filled = eliminationgraph(graph);

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
```
"""
function eliminationgraph(graph; alg::PermutationOrAlgorithm=DEFAULT_ELIMINATION_ALGORITHM, snd::SupernodeType=DEFAULT_SUPERNODE_TYPE)
    eliminationgraph(true, graph, alg, snd)
end


function eliminationgraph(element, graph; alg::PermutationOrAlgorithm=DEFAULT_ELIMINATION_ALGORITHM, snd::SupernodeType=DEFAULT_SUPERNODE_TYPE)
    eliminationgraph(element, graph, alg, snd)
end


function eliminationgraph(element::T, graph, alg::PermutationOrAlgorithm, snd::SupernodeType) where T
    label, matrix = eliminationgraph(T, graph, alg, snd)
    fill!(nonzeros(matrix), element)
    label, matrix
end


function eliminationgraph(T::Type, graph, alg::PermutationOrAlgorithm, snd::SupernodeType)
    label, tree = junctiontree(graph, alg, snd)
    label, eliminationgraph(T, tree)
end


function eliminationgraph(T::Type, graph, alg::PermutationOrAlgorithm, snd::Nodal)
    label, tree = junctiontree(graph, alg, snd)
    I = eltype(eltype(tree))
    nzval = Vector{T}(undef, length(tree.sepval))
    label, SparseMatrixCSC{T, I}(length(tree), length(tree), tree.sepptr, tree.sepval, nzval)
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
    ischordal(graph)

Determine whether a simple graph is [chordal](https://en.wikipedia.org/wiki/Chordal_graph).
"""
function ischordal(matrix::SparseMatrixCSC{<:Any, I}) where I
    # validate arguments
    size(matrix, 1) != size(matrix, 2) && throw(ArgumentError("size(matrix, 1) != size(matrix, 2)"))

    # run algorithm
    vertices::OneTo{I} = axes(matrix, 2)

    ischordal(vertices) do j
        @view rowvals(matrix)[nzrange(matrix, j)]
    end
end


function ischordal(neighbors::Function, vertices::AbstractVector)
    index = mcs(neighbors, vertices)
    isperfect(neighbors, invperm(index), index)
end


"""
    isfilled(graph)

Determine whether a directed graph is filled.
"""
function isfilled(matrix::SparseMatrixCSC{<:Any, I}) where I
    # validate arguments
    size(matrix, 1) != size(matrix, 2) && throw(ArgumentError("size(matrix, 1) != size(matrix, 2)"))

    # run algorithm
    vertices::OneTo{I} = axes(matrix, 2)

    isfilled(vertices) do j
        @view rowvals(matrix)[nzrange(matrix, j)]
    end
end


function isfilled(neighbors::Function, vertices::AbstractVector)
    isperfect(neighbors, vertices, vertices)
end


"""
    isperfect(graph, order::AbstractVector[, index::AbstractVector])

Determine whether an fill-reducing permutation is perfect.
"""
function isperfect(graph, order::AbstractVector)
    isperfect(graph, order, invperm(order))
end


function isperfect(matrix::SparseMatrixCSC{<:Any, I}, order::AbstractVector{I}, index::AbstractVector{I}) where I
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
function isperfect(neighbors::Function, order::AbstractVector{I}, index::AbstractVector{I}) where I
    # validate arguments
    eachindex(order) != eachindex(index) && throw(ArgumentError("eachindex(order) != eachindex(index)"))

    # run algorithm
    f = Vector{I}(undef, length(order))
    findex = Vector{I}(undef, length(order))

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
