# An ordered graph (G, σ).
struct OrderedGraph
    graph::Graph
    order::Order
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


# Construct the elimination graph of an ordered graph.
function eliminationgraph(ograph::OrderedGraph)
    ograph = deepcopy(ograph)
    
    for i in vertices(ograph)
        ns = collect(outneighbors(ograph, i))
        n = length(ns)
        
        for j in 1:n
            for k in j + 1:n
                add_edge!(ograph, ns[j], ns[k])
            end
        end    
    end
    
    ograph
end

    
function Base.deepcopy(ograph::OrderedGraph)
    order = deepcopy(ograph.order)
    graph = deepcopy(ograph.graph)
    OrderedGraph(graph, order)
end


# Get the vertex σ(i).
function order(ograph::OrderedGraph, i)
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


function BasicGraphs.has_edge(ograph::OrderedGraph, i, j)
    u = min(i, j)
    v = max(i, j)
    has_edge(ograph.graph, u, v)
end


function BasicGraphs.add_edge!(ograph::OrderedGraph, i, j)
    u = min(i, j)
    v = max(i, j)
    
    if !has_edge(ograph.graph, u, v)
        add_edge!(ograph.graph, u, v)
        true
    else
        false
    end
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
