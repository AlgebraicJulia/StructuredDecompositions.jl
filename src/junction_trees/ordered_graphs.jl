# An ordered graph (G, σ).
struct OrderedGraph
    lower::SparseMatrixCSC{Int64, Int64} # adjacency matrix (lower triangular)
    upper::SparseMatrixCSC{Int64, Int64} # adjacency matrix (upper triangular)
    order::Order                         # permutation
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
    m = last(nzrange(graph, size(graph, 1)))
    n = size(graph, 1)

    colptr_lower = Vector{Int}(undef, n + 1)
    colptr_upper = Vector{Int}(undef, n + 1)

    rowval_lower = Vector{Int}(undef, m ÷ 2)
    rowval_upper = Vector{Int}(undef, m ÷ 2)

    colptr_lower[1] = 1
    colptr_upper[1] = 1

    count_lower = 1
    count_upper = 1

    for i in 1:n
        colptr_lower[i] = count_lower
        colptr_upper[i] = count_upper
        neighbors = inverse(order, rowvals(graph)[nzrange(graph, order[i])])
        sort!(neighbors)

        for j in neighbors
            if i < j
                rowval_lower[count_lower] = j                
                count_lower += 1
            else
                rowval_upper[count_upper] = j
                count_upper += 1
            end
        end
    end

    colptr_lower[n + 1] = m ÷ 2 + 1
    colptr_upper[n + 1] = m ÷ 2 + 1  

    nzval_lower = ones(Int, m ÷ 2)
    nzval_upper = ones(Int, m ÷ 2) 

    lower = SparseMatrixCSC(n, n, colptr_lower, rowval_lower, nzval_lower)
    upper = SparseMatrixCSC(n, n, colptr_upper, rowval_upper, nzval_upper)

    OrderedGraph(lower, upper, order)
end


# Given an ordered graph (G, σ) and permutation μ, construct the ordered graph
#    (G, σ ∘ μ).
# ----------------------------------------
#    graph     ordered graph
#    order     permutation
# ----------------------------------------
function OrderedGraph(graph::OrderedGraph, order::Order)
    newgraph = OrderedGraph(adjacencymatrix(graph), order)
    OrderedGraph(newgraph.lower, newgraph.upper, compose(order, graph.order))
end


# Construct the permutation σ.
function Order(graph::OrderedGraph)
    Order(graph.order)
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
    last(graph.lower.colptr) - 1
end


function BasicGraphs.nv(graph::OrderedGraph)
    size(graph.lower, 1)
end


function BasicGraphs.inneighbors(graph::OrderedGraph, i)
    rowvals(graph.upper)[nzrange(graph.upper, i)]
end


function BasicGraphs.outneighbors(graph::OrderedGraph, i)
    rowvals(graph.lower)[nzrange(graph.lower, i)]
end


function BasicGraphs.all_neighbors(graph::OrderedGraph, i)
    Base.Iterators.flatten((inneighbors(graph, i), outneighbors(graph, i)))
end


#function BasicGraphs.edges(graph::OrderedGraph)
#    1:ne(graph)
#end


function BasicGraphs.vertices(graph::OrderedGraph)
    1:nv(graph)
end


#function BasicGraphs.src(graph::OrderedGraph, i)
#    src(graph.graph, i)
#end


#function BasicGraphs.tgt(graph::OrderedGraph, i)
#    tgt(graph.graph, i)
#end
