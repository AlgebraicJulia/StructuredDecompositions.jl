module TreeWidthSolverExt


using TreeWidthSolver
using TreeWidthSolver.Graphs
using SparseArrays
using StructuredDecompositions.JunctionTrees


function JunctionTrees.permutation(graph::Graph{I}, alg::BT) where I
    order::Vector{I} = reverse!(reduce(vcat, elimination_order(graph); init=Int[]))
    order, invperm(order)
end


function JunctionTrees.permutation(matrix::SparseMatrixCSC{<:Any, I}, alg::BT) where I
    order::Vector{I}, index::Vector{I} = permutation(simple_graph(matrix), alg)
    order, index
end


end
