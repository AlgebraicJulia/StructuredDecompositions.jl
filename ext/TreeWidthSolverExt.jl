module TreeWidthSolverExt


import TreeWidthSolver


using SparseArrays
using StructuredDecompositions.JunctionTrees


function JunctionTrees.permutation(matrix::SparseMatrixCSC, alg::BT)
    graph = TreeWidthSolver.simple_graph(matrix)
    order = reverse!(reduce(vcat, TreeWidthSolver.elimination_order(graph); init=Int[]))
    order, invperm(order)
end


end
