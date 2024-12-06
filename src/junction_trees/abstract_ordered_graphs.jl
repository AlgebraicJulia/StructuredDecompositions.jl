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


# Multiline printing.
function Base.show(io::IO, ::MIME"text/plain", graph::AbstractOrderedGraph)
    print(io, "ordered graph:\n")
    SparseArrays._show_with_braille_patterns(io, adjacencymatrix(graph))
end 

