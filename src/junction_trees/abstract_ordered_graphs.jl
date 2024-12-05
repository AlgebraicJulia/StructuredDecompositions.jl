abstract type AbstractOrderedGraph <: AbstractSimpleGraph{Int} end


############################
# Abstract Graph Interface #
############################


function SimpleGraphs.is_directed(::Type{AbstractOrderedGraph})
    true
end


function SimpleGraphs.edgetype(graph::AbstractOrderedGraph)
    SimpleEdge{Int}
end


function SimpleGraphs.has_edge(graph::AbstractOrderedGraph, edge::SimpleEdge{Int})
    i = src(edge)
    j = dst(edge)
    i < j && insorted(j, outneighbors(graph, i))
end
