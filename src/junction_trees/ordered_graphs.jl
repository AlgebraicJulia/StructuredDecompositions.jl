"""
    OrderedGraph <: AbstractSimpleGraph{Int}

An [ordered graph](https://en.wikipedia.org/wiki/Ordered_graph) ``(G, \\sigma)``.
This type implements the [abstract graph interface](https://juliagraphs.org/Graphs.jl/stable/core_functions/interface/).
"""
struct OrderedGraph <: AbstractSimpleGraph{Int}
    lower::SparseMatrixCSC{Bool, Int} # adjacency matrix (lower triangular)
    upper::SparseMatrixCSC{Bool, Int} # adjacency matrix (upper triangular)
    order::Order                      # permutation
end


"""
    OrderedGraph(graph[, ealg::Union{Order, EliminationAlgorithm}])

Construct an ordered graph, optionally specifying an elimination algorithm.
"""
function OrderedGraph(graph, ealg::Union{Order, EliminationAlgorithm}=DEFAULT_ELIMINATION_ALGORITHM)
    OrderedGraph(adjacencymatrix(graph), ealg)
end


# Given a graph G, construct the ordered graph
#    (G, σ),
# where σ is a permutation computed using an elimination algorithm.
# ----------------------------------------
#    graph     simple connected graph
#    ealg      elimination algorithm
# ----------------------------------------
function OrderedGraph(graph::AbstractMatrix, ealg::EliminationAlgorithm=DEFAULT_ELIMINATION_ALGORITHM)
    OrderedGraph(graph, Order(graph, ealg))
end


# Given a graph H and permutation σ, construct the ordered graph
#    (G, σ)
# ----------------------------------------
#    graph     simple connected graph
#    order     vertex order
# ----------------------------------------
function OrderedGraph(graph::AbstractSparseMatrixCSC, order::Order)
    graph = permute(graph, order, order)
    OrderedGraph(tril(graph), triu(graph), order)
end


# Given an ordered graph (G, σ) and permutation μ, construct the ordered graph
#    (G, σ ∘ μ).
# ----------------------------------------
#    graph           ordered graph
#    permutation     permutation
# ----------------------------------------
function OrderedGraph(graph::OrderedGraph, permutation::Order)
    order = graph.order
    graph = OrderedGraph(adjacencymatrix(graph), permutation)
    OrderedGraph(graph.lower, graph.upper, compose(permutation, order))
end


"""
    Order(graph::OrderedGraph)

Construct the permutation ``\\sigma``.
"""
function Order(graph::OrderedGraph)
    copy(graph.order)
end


# Construct the adjacency matrix of an ordered graph.
function adjacencymatrix(graph::OrderedGraph)
    graph.lower .|| graph.upper
end


# A Compact Row Storage Scheme for Cholesky Factors Using Elimination Trees
# Liu
# Algorithm 4.2: Elimination Tree by Path Compression.
function etree(graph::OrderedGraph)
    n = nv(graph)
    parent = Vector{Int}(undef, n)
    ancestor = Vector{Int}(undef, n)

    for i in 1:n
        parent[i] = 0
        ancestor[i] = 0

        for k in inneighbors(graph, i)
            r = k

            while !iszero(ancestor[r]) && ancestor[r] != i
                t = ancestor[r]
                ancestor[r] = i
                r = t
            end

            if iszero(ancestor[r])
                ancestor[r] = i
                parent[r] = i
            end
        end
    end

    parent[n] = n
    parent
end


# Construct a copy of an ordered graph.
function Base.copy(graph::OrderedGraph)
    OrderedGraph(graph.lower, graph.upper, graph.order)
end


# Multiline printing.
function Base.show(io::IO, ::MIME"text/plain", graph::OrderedGraph)
    print(io, "ordered graph:\n")
    SparseArrays._show_with_braille_patterns(io, graph.lower)
end 


############################
# Abstract Graph Interface #
############################


function SimpleGraphs.is_directed(::Type{OrderedGraph})
    true
end


function SimpleGraphs.edgetype(graph::OrderedGraph)
    SimpleEdge{Int}
end


function SimpleGraphs.ne(graph::OrderedGraph)
    last(graph.lower.colptr) - 1
end


function SimpleGraphs.nv(graph::OrderedGraph)
    size(graph.lower, 1)
end


function SimpleGraphs.badj(graph::OrderedGraph, i::Integer)
    view(rowvals(graph.upper), nzrange(graph.upper, i))
end


function SimpleGraphs.fadj(graph::OrderedGraph, i::Integer)
    view(rowvals(graph.lower), nzrange(graph.lower, i))
end


function SimpleGraphs.all_neighbors(graph::OrderedGraph, i::Integer)
    [inneighbors(graph, i); outneighbors(graph, i)]
end


function SimpleGraphs.has_edge(graph::OrderedGraph, edge::SimpleEdge{Int})
    i = src(edge)
    j = dst(edge)
    i < j && insorted(j, outneighbors(graph, i))
end
