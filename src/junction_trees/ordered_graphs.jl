"""
    OrderedGraph <: AbstractOrderedGraph

A directed simple graph whose edges ``(i, j)`` satisfy the inequality ``i < j``.
This type implements the [abstract graph interface](https://juliagraphs.org/Graphs.jl/stable/core_functions/interface/).
"""
struct OrderedGraph <: AbstractSimpleGraph{Int}
    matrix::SparseMatrixCSC{Bool, Int}
    colptr::Vector{Int}
end


"""
    OrderedGraph(graph[, ealg::Union{Order, EliminationAlgorithm}])

Construct an ordered graph by permuting the vertices of a simple graph and directing them from lower to higher.
"""
function OrderedGraph(graph, ealg::Union{Order, EliminationAlgorithm}=DEFAULT_ELIMINATION_ALGORITHM)
    OrderedGraph(adjacencymatrix(graph), ealg)
end


# Construct an ordered graph by permuting the vertices of a simple graph and directing them from lower to higher.
# ----------------------------------------
#    graph     simple graph
#    ealg      elimination algorithm
# ----------------------------------------
function OrderedGraph(graph::AbstractMatrix, ealg::EliminationAlgorithm=DEFAULT_ELIMINATION_ALGORITHM)
    OrderedGraph(graph, Order(graph, ealg))
end


# Construct an ordered graph by permuting the vertices of a simple graph and directing them from lower to higher.
# ----------------------------------------
#    graph     simple graph
#    order     vertex order
# ----------------------------------------
function OrderedGraph(graph::AbstractSparseMatrixCSC, order::Order)
    OrderedGraph(permute(graph, order, order))
end


function OrderedGraph(graph::AbstractSparseMatrixCSC)
    n = size(graph, 1)
    colptr = Vector{Int}(undef, n)

    for i in 1:n
        colptr[i] = graph.colptr[i] + searchsortedfirst(view(rowvals(graph), nzrange(graph, i)), i) - 1
    end

    OrderedGraph(graph, colptr)
end


function Base.permute!(graph::OrderedGraph, permutation::AbstractVector)
    permute!(graph.matrix, permutation, permutation)

    map!(graph.colptr, vertices(graph)) do i
        first(nzrange(graph.matrix, i)) + searchsortedfirst(all_neighbors(graph, i), i) - 1
    end

    graph 
end


# Construct the adjacency matrix of an ordered graph.
function adjacencymatrix(graph::OrderedGraph)
    graph.matrix
end


# A Compact Row Storage Scheme for Cholesky Factors Using Elimination Trees
# Liu
# Algorithm 4.2: Elimination Tree by Path Compression.
# ----------------------------------------
#    graph     simple connected graph
# ----------------------------------------
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

    Tree(parent)
end


function etree!(order::Order, graph::OrderedGraph)
    tree = etree(graph)
    postorder = postorder!(tree)
    permute!(order, postorder)
    permute!(graph, postorder)
    tree
end


# An Efficient Algorithm to Compute Row and Column Counts for Sparse Cholesky Factorization
# Gilbert, Ng, and Peyton
# Figure 3: Implementation of algorithm to compute row and column counts.
# ----------------------------------------
#    graph     simple connected graph
#    tree      elimination tree
# ----------------------------------------
function supcnt(graph::OrderedGraph, tree::Tree, level::AbstractVector=levels(tree), fdesc::AbstractVector=firstdescendants(tree))
    n = treesize(tree)
    sets = DisjointSets(n)
    
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
            if fdesc[p] > prev_nbr[u]
                wt[p] += 1
                pp = prev_p[u]
                
                if iszero(pp)
                    rc[u] += level[p] - level[u]
                else
                    q = find!(sets, pp)
                    rc[u] += level[p] - level[q]
                    wt[q] -= 1
                end
    
                prev_p[u] = p
            end

            prev_nbr[u] = p
        end

        union!(sets, p, parentindex(tree, p))
    end

    cc = wt

    for v in 1:n - 1
        cc[parentindex(tree, v)] += cc[v]
    end

    rc, cc
end


# Multiline printing.
function Base.show(io::IO, ::MIME"text/plain", graph::OrderedGraph)
    print(io, "ordered graph:\n")
    SparseArrays._show_with_braille_patterns(io, adjacencymatrix(graph))
end 


############################
# Abstract Graph Interface #
############################


function SimpleGraphs.ne(graph::OrderedGraph)
    nnz(graph.matrix) ÷ 2
end


function SimpleGraphs.nv(graph::OrderedGraph)
    size(graph.matrix, 1)
end


function SimpleGraphs.badj(graph::OrderedGraph, i::Integer)
    view(rowvals(graph.matrix), graph.matrix.colptr[i]:graph.colptr[i] - 1)
end


function SimpleGraphs.fadj(graph::OrderedGraph, i::Integer)
    view(rowvals(graph.matrix), graph.colptr[i]:graph.matrix.colptr[i + 1] - 1)
end


function SimpleGraphs.all_neighbors(graph::OrderedGraph, i::Integer)
    view(rowvals(graph.matrix), nzrange(graph.matrix, i))
end


function SimpleGraphs.is_directed(::Type{OrderedGraph})
    true
end


function SimpleGraphs.edgetype(graph::OrderedGraph)
    SimpleEdge{Int}
end


function SimpleGraphs.has_edge(graph::OrderedGraph, edge::SimpleEdge{Int})
    i = src(edge)
    j = dst(edge)
    i < j && insorted(j, outneighbors(graph, i))
end
