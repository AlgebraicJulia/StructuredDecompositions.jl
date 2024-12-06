struct FilledGraph <: AbstractOrderedGraph
    colptr::Vector{Int}
    rowval::Vector{Int}
    tree::PostorderTree
end


function adjacencymatrix(graph::FilledGraph)
    lower = SparseMatrixCSC(nv(graph), nv(graph), graph.colptr, graph.rowval, ones(Bool, ne(graph)))
    lower .|| lower'
end


############################
# Abstract Graph Interface #
############################


function SimpleGraphs.ne(graph::FilledGraph)
    last(graph.colptr) - 1
end


function SimpleGraphs.nv(graph::FilledGraph)
    length(graph.colptr)
end


function SimpleGraphs.badj(graph::FilledGraph, i::Integer)
    descendantindices(graph.tree, i)
end


function SimpleGraphs.fadj(graph::FilledGraph, i::Integer)
    view(graph.rowval, graph.colptr[i]:graph.colptr[i + 1] - 1)
end


#=
function SimpleGraphs.all_neighbors(graph::OrderedGraph, i::Integer)
    view(graph.rowval, graph.adjptr[i]:graph.adjptr[i + 1] - 1)
end
=#
