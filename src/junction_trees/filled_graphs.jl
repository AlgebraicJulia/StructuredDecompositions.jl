struct FilledGraph <: AbstractSimpleGraph{Int}
    lower::SparseMatrixCSC{Bool, Int}
    tree::PostorderTree
end


############################
# Abstract Graph Interface #
############################


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
