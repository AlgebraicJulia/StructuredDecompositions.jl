"""
    OrderedGraph <: AbstractOrderedGraph

A directed simple graph whose edges ``(i, j)`` satisfy the inequality ``i < j``.
This type implements the [abstract graph interface](https://juliagraphs.org/Graphs.jl/stable/core_functions/interface/).
"""
struct OrderedGraph <: AbstractOrderedGraph
    colptr::Vector{Int}
    adjptr::Vector{Int}
    rowval::Vector{Int}
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

    OrderedGraph(colptr, graph.colptr, graph.rowval)
end



# Construct the adjacency matrix of an ordered graph.
function adjacencymatrix(graph::OrderedGraph)
    SparseMatrixCSC(nv(graph), nv(graph), graph.adjptr, graph.rowval, ones(Bool, length(graph.rowval)))
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

    parent
end


# An Efficient Algorithm to Compute Row and Column Counts for Sparse Cholesky Factorization
# Gilbert, Ng, and Peyton
# Figure 3: Implementation of algorithm to compute row and column counts.
# ----------------------------------------
#    graph     simple connected graph
#    tree      elimination tree
# ----------------------------------------
function supcnt(graph::OrderedGraph, tree::Tree)
    order = postorder(tree)
    rc, cc = supcnt(OrderedGraph(graph, order), PostorderTree(tree, order))
    view(rc, inv(order)), view(cc, inv(order))
end


# An Efficient Algorithm to Compute Row and Column Counts for Sparse Cholesky Factorization
# Gilbert, Ng, and Peyton
# Figure 3: Implementation of algorithm to compute row and column counts.
# ----------------------------------------
#    graph     simple connected graph
#    tree      elimination tree
# ----------------------------------------
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
            if first(descendantindices(tree, p)) > prev_nbr[u]
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


# Multiline printing.
function Base.show(io::IO, ::MIME"text/plain", graph::OrderedGraph)
    print(io, "ordered graph:\n")
    SparseArrays._show_with_braille_patterns(io, adjacencymatrix(graph))
end 


############################
# Abstract Graph Interface #
############################


function SimpleGraphs.ne(graph::OrderedGraph)
    (last(graph.adjptr) - 1) ÷ 2
end


function SimpleGraphs.nv(graph::OrderedGraph)
    length(graph.colptr)
end


function SimpleGraphs.badj(graph::OrderedGraph, i::Integer)
    view(graph.rowval, graph.adjptr[i]:graph.colptr[i] - 1)
end


function SimpleGraphs.fadj(graph::OrderedGraph, i::Integer)
    view(graph.rowval, graph.colptr[i]:graph.adjptr[i + 1] - 1)
end


function SimpleGraphs.all_neighbors(graph::OrderedGraph, i::Integer)
    view(graph.rowval, graph.adjptr[i]:graph.adjptr[i + 1] - 1)
end
