# The edge set of an ordered graoh.
struct OrderedGraphEdges <: Graphs.AbstractEdgeIter
    graph::OrderedGraph
end


function Base.show(io::IO, edges::OrderedGraphEdges)
    n = length(edges)
    print(io, "OrderedGraphEdges $n")
end 


######################
# Iterator Interface #
######################


function Base.iterate(edges::OrderedGraphEdges, (i, p)::Tuple{Integer, Integer}=(1, 1))
    if p <= ne(edges.graph)
        while p ∉ nzrange(edges.graph.upper, i)
            i += 1
        end

        Graphs.SimpleEdge{Int}(i, edges.graph.upper.rowval[p]), (i, p + 1)
    else
        nothing
    end
end


function Base.iterate(reversed::Iterators.Reverse{OrderedGraphEdges}, (i, p)::Tuple{Integer, Integer}=(nv(reversed.itr.graph), ne(reversed.itr.graph)))
    if p >= 1
        while p ∉ nzrange(reversed.itr.graph.upper, i)
            i -= 1
        end

        Graphs.SimpleEdge{Int}(i, reversed.itr.graph.upper.rowval[p]), (i, p - 1)
    else
        nothing
    end
end


function Base.eltype(::Type{OrderedGraphEdges})
    Graphs.SimpleEdge{Int}
end


function Base.hasfastin(::Type{OrderedGraphEdges})
    true
end


function Base.length(edges::OrderedGraphEdges)
    ne(edges.graph)
end


function Base.in(edge, edges::OrderedGraphEdges)
    has_edge(edges.graph, edge)
end


############################
# Abstract Graph Interface #
############################


function Graphs.edges(graph::OrderedGraph)
    OrderedGraphEdges(graph)
end
