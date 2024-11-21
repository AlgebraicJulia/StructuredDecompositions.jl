# An ordered graph (G, σ).
struct OrderedGraph
    graph::Graph # graph
    order::Order # permutation
end


# Given a graph G, construct the ordered graph
#    (G, σ),
# where σ is a permutation computed using an elimination algorithm.
# ----------------------------------------
#    graph     simple connected graph
#    ealg      elimination algorithm
# ----------------------------------------
function OrderedGraph(graph, ealg::Union{Order, EliminationAlgorithm}=DEFAULT_ELIMINATION_ALGORITHM)
    OrderedGraph(adjacencymatrix(graph), ealg)
end


# Given a matrix M, construct the ordered graph
#    (G, σ),
# where G is the sparsity graph if M and σ is a permutation computed using an elimination
# algorithm.
# ----------------------------------------
#    matrix    symmetric matrix
#    ealg      elimination algorithm
# ----------------------------------------
function OrderedGraph(matrix::AbstractMatrix, ealg::EliminationAlgorithm=DEFAULT_ELIMINATION_ALGORITHM)
    OrderedGraph(matrix, Order(matrix, ealg))
end


# Given a matrix M and permutation σ, construct the ordered graph
#    (G, σ)
# where G is the sparsity graph of M.
# ----------------------------------------
#    graph     symmetric matrix
#    order     vertex order
# ----------------------------------------
function OrderedGraph(matrix::SparseMatrixCSC, order::Order)
    n = size(matrix, 1)
    graph = Graph(n)

    for u in 1:n
        for v in matrix.rowval[matrix.colptr[u]:matrix.colptr[u + 1] - 1]
            i = inverse(order, u)
            j = inverse(order, v)

            if i < j
                add_edge!(graph, i, j)
            end
        end
    end

    OrderedGraph(graph, order)
end


# Given an ordered graph (G, σ) and permutation μ, construct the ordered graph
#    (G, σ ∘ μ).
# ----------------------------------------
#    graph     ordered graph
#    order     permutation
# ----------------------------------------
function OrderedGraph(graph::OrderedGraph, order::Order)
    digraph = Graph(nv(graph))

    for edge in edges(graph)
        i = inverse(order, src(graph, edge))
        j = inverse(order, tgt(graph, edge))

        if i < j
            add_edge!(digraph, i, j)
        else
            add_edge!(digraph, j, i)
        end
    end
    
    OrderedGraph(digraph, compose(order, graph.order))
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

            while ancestor[r] != 0 && ancestor[r] != i
                t = ancestor[r]
                ancestor[r] = i
                r = t
            end

            if ancestor[r] == 0
                ancestor[r] = i
                parent[r] = i
            end
        end
    end

    parent[n] = n
    parent
end

    
function Base.deepcopy(graph::OrderedGraph)
    OrderedGraph(deepcopy(graph.order), deepcopy(graph.graph))
end


# Get the vertex σ(i).
function permutation(graph::OrderedGraph, i)
    graph.order[i]
end


# Get the index σ⁻¹(v).
function inverse(graph::OrderedGraph, v)
    inverse(graph.order, v)
end


############################
# Abstract Graph Interface #
############################


function BasicGraphs.ne(graph::OrderedGraph)
    ne(graph.graph)
end


function BasicGraphs.nv(graph::OrderedGraph)
    nv(graph.graph)
end


function BasicGraphs.inneighbors(graph::OrderedGraph, i)
    inneighbors(graph.graph, i)
end


function BasicGraphs.outneighbors(graph::OrderedGraph, i)
    outneighbors(graph.graph, i)
end


function BasicGraphs.edges(graph::OrderedGraph)
    edges(graph.graph)
end


function BasicGraphs.vertices(graph::OrderedGraph)
    vertices(graph.graph)
end
    

function BasicGraphs.src(graph::OrderedGraph, i)
    src(graph.graph, i)
end


function BasicGraphs.tgt(graph::OrderedGraph, i)
    tgt(graph.graph, i)
end
