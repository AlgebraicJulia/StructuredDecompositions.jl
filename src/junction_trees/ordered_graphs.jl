# An ordered graph (G, σ).
struct OrderedGraph <: AbstractSimpleGraph{Int}
    lower::SparseMatrixCSC{Bool, Int} # adjacency matrix (lower triangular)
    upper::SparseMatrixCSC{Bool, Int} # adjacency matrix (upper triangular)
    order::Order                      # permutation
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

    count_lower = 1
    count_upper = 1

    for i in 1:n
        colptr_lower[i] = count_lower
        colptr_upper[i] = count_upper
        neighbor = sort(inverse(order, rowvals(graph)[nzrange(graph, order[i])]))

        for j in neighbor
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

    nzval_lower = ones(Bool, m ÷ 2)
    nzval_upper = ones(Bool, m ÷ 2) 

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


# Construct the adjacency matrix of an ordered graph.
function adjacencymatrix(graph::OrderedGraph)
    m = ne(graph)
    n = nv(graph)
    colptr = Vector{Int}(undef, n + 1)
    rowval = Vector{Int}(undef, 2m)
    count = 1

    for i in 1:n
        colptr[i] = count

        for j in sort(all_neighbors(graph, i))
            rowval[count] = j
            count += 1
        end
    end

    colptr[n + 1] = 2m + 1
    nzval = ones(Bool, 2m)
    SparseMatrixCSC(n, n, colptr, rowval, nzval)
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
    permutation(graph.order, i)
end


# Get the index σ⁻¹(v).
function inverse(graph::OrderedGraph, v)
    inverse(graph.order, v)
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
