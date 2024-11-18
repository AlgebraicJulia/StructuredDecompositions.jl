# An ordered graph (G, σ).
struct OrderedGraph
    graph::Graph # graph
    order::Order # permutation
end


# Given a graph G, construct the ordered graph
#    (G, σ),
# where the permutation σ is computed using an elimination algorithm.
function OrderedGraph(sgraph::AbstractSymmetricGraph, ealg::EliminationAlgorithm=DEFAULT_ELIMINATION_ALGORITHM)
    OrderedGraph(sgraph, Order(sgraph, ealg))
end


# Given a graph G and permutation σ, construct the ordered graph
#    (G, σ).
function OrderedGraph(sgraph::AbstractSymmetricGraph, order::Order)
    n = nv(sgraph)
    graph = Graph(n)

    for e in edges(sgraph)
        u = src(sgraph, e)
        v = tgt(sgraph, e)
        
        if order(u, v)
            add_edge!(graph, inverse(order, u), inverse(order, v))
        end
    end

    OrderedGraph(graph, order)
end


# Given an ordered graph (G, σ) and permutation μ, construct the ordered graph
#    (G, σ ∘ μ).
function OrderedGraph(ograph::OrderedGraph, order::Order)
    n = nv(ograph)
    graph = Graph(n)

    for e in edges(ograph)
        u = src(ograph, e)
        v = tgt(ograph, e)

        if order(u, v)
            add_edge!(graph, inverse(order, u), inverse(order, v))
        else
            add_edge!(graph, inverse(order, v), inverse(order, u))
        end
    end
    
    OrderedGraph(graph, compose(order, ograph.order))
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

    
function Base.deepcopy(ograph::OrderedGraph)
    order = deepcopy(ograph.order)
    graph = deepcopy(ograph.graph)
    OrderedGraph(graph, order)
end


# Get the vertex σ(i).
function permutation(ograph::OrderedGraph, i)
    ograph.order[i]
end


# Get the index σ⁻¹(v),
function inverse(ograph::OrderedGraph, v)
    inverse(ograph.order, v)
end


############################
# Abstract Graph Interface #
############################


function BasicGraphs.nv(ograph::OrderedGraph)
    nv(ograph.graph)
end


function BasicGraphs.inneighbors(ograph::OrderedGraph, i)
    inneighbors(ograph.graph, i)
end


function BasicGraphs.outneighbors(ograph::OrderedGraph, i)
    outneighbors(ograph.graph, i)
end


function BasicGraphs.edges(ograph::OrderedGraph)
    edges(ograph.graph)
end


function BasicGraphs.vertices(ograph::OrderedGraph)
    vertices(ograph.graph)
end
    

function BasicGraphs.src(ograph::OrderedGraph, i)
    src(ograph.graph, i)
end


function BasicGraphs.tgt(ograph::OrderedGraph, i)
    tgt(ograph.graph, i)
end
