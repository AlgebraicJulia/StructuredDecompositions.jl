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


function eliminationgraph(T::Type, tree::JunctionTree{<:Any, E}) where E
    n = residual(tree[end])[end]
    eliminationgraph!(spzeros(T, E, n, n), tree)
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
    I = promote_type(vtype(tree), etype(tree))
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
function ischordal(graph)
    ischordal(Graph(graph))
end


function ischordal(graph::Graph)
    index = mcs(graph)
    isperfect(graph, invperm(index), index)
end


"""
    isfilled(graph)

Determine whether a directed graph is filled.
"""
function isfilled(graph)
    isfilled(Graph(graph))
end


function isfilled(graph::Graph)
    isperfect(graph, vertices(graph), vertices(graph))
end


"""
    isperfect(graph, order::AbstractVector[, index::AbstractVector])

Determine whether an fill-reducing permutation is perfect.
"""
function isperfect(graph, order::AbstractVector, index::AbstractVector=invperm(order))
    isperfect(Graph(graph), order, index)
end


# Simple Linear-Time Algorithms to Test Chordality of Graphs, Test Acyclicity of Hypergraphs, and Selectively Reduce Acyclic Hypergraphs
# Tarjan and Yannakakis
# Test for Zero Fill-In.
#
# Determine whether a fill-reducing permutation is perfect.
function isperfect(graph::Graph{V}, order::AbstractVector{V}, index::AbstractVector{V}) where V
    # validate arguments
    vertices(graph) != eachindex(index) && throw(ArgumentError("vertices(graph) != eachindex(index)"))
    eachindex(order) != eachindex(index) && throw(ArgumentError("eachindex(order) != eachindex(index)"))

    # run algorithm
    f = Vector{V}(undef, nv(graph))
    findex = Vector{V}(undef, nv(graph))

    for (i, w) in enumerate(order)
        f[w] = w
        findex[w] = i

        for v in neighbors(graph, w)
            if index[v] < i
                findex[v] = i

                if f[v] == v
                    f[v] = w
                end
            end
        end

        for v in neighbors(graph, w)
            if index[v] < i && findex[f[v]] < i
                return false
            end
        end
    end

    true
end
