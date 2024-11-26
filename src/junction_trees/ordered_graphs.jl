"""
    OrderedGraph <: AbstractSimpleGraph{Int}

A directed simple graph whose edges i → j are oriented from lower to higher vertices.
This type implements the [abstract graph interface](https://juliagraphs.org/Graphs.jl/stable/core_functions/interface/).
"""
struct OrderedGraph <: AbstractSimpleGraph{Int}
    symmetric::SparseMatrixCSC{Bool, Int} # adjacency matrix (symmetric)
    lower::SparseMatrixCSC{Bool, Int}     # adjacency matrix (lower triangular)
    upper::SparseMatrixCSC{Bool, Int}     # adjacency matrix (upper triangular)
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
    OrderedGraph(graph, tril(graph), triu(graph))
end


function adjacencymatrix(graph::OrderedGraph)
    graph.symmetric
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

    parent
end


# An Efficient Algorithm to Compute Row and Column Counts for Sparse Cholesky Factorization
# Gilbert, Ng, and Peyton
# Figure 3: Implementation of algorithm to compute row and column counts.
function supcnt(graph::OrderedGraph, tree::Tree)
    order = postorder(tree)
    rc, cc = supcnt(OrderedGraph(graph, order), PostorderTree(tree, order))
    view(rc, inv(order)), view(cc, inv(order))
end


# An Efficient Algorithm to Compute Row and Column Counts for Sparse Cholesky Factorization
# Gilbert, Ng, and Peyton
# Figure 3: Implementation of algorithm to compute row and column counts.
function supcnt(graph::OrderedGraph, tree::PostorderTree)
    n = treesize(tree)
    
    #### Disjoint Set Union ####
    
    rvert = collect(1:n)
    index = collect(1:n)
    forest = IntDisjointSets(n)
    
    function find(u)
        index[find_root!(forest, u)]
    end
    
    function union(u, v)
        w = max(u, v)
        rvert[w] = root_union!(forest, rvert[u], rvert[v])
        index[rvert[w]] = w
    end
    
    ############################
    
    prev_p = zeros(Int, n)
    prev_nbr = zeros(Int, n)
    rc = ones(Int, n)
    wt = ones(Int, n)

    for u in 1:n - 1
        wt[parentindex(tree, u)] = 0
    end
    
    for p in 1:n - 1
        wt[parentindex(tree, p)] -= 1

        for u in outneighbors(graph, p)
            if fdesc(tree, p) > prev_nbr[u]
                wt[p] += 1
                pp = prev_p[u]
                
                if iszero(pp)
                    rc[u] += level(tree, p) - level(tree, u)
                else
                    q = find(pp)
                    rc[u] += level(tree, p) - level(tree, q)
                    wt[q] -= 1
                end
    
                prev_p[u] = p
            end

            prev_nbr[u] = p
        end

        union(p, parentindex(tree, p))
    end

    cc = wt

    for v in 1:n - 1
        cc[parentindex(tree, v)] += cc[v]
    end

    rc, cc
end


# Compute higher degree of every vertex in the elimination graph of
#    (G, σ).
function outdegrees(graph::OrderedGraph, tree::Tree)
    rc, cc = supcnt(graph, tree)
    cc .- 1
end


# Construct a copy of an ordered graph.
function Base.copy(graph::OrderedGraph)
    OrderedGraph(graph.symmetric, graph.lower, graph.upper, graph.order)
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
    view(rowvals(graph.symmetric), nzrange(graph.symmetric, i))
end


function SimpleGraphs.has_edge(graph::OrderedGraph, edge::SimpleEdge{Int})
    i = src(edge)
    j = dst(edge)
    i < j && insorted(j, outneighbors(graph, i))
end
