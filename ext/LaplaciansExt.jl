module LaplaciansExt


import Laplacians


using SparseArrays
using StructuredDecompositions.JunctionTrees


function JunctionTrees.permutation(matrix::SparseMatrixCSC, alg::Spectral)
    order = spectralorder(matrix; tol=alg.tol)
    order, invperm(order)
end


# A Spectral Algorithm for Envelope Reduction of Sparse Matrices
# Barnard, Pothen, and Simon
# Algorithm 1: Spectral Algorithm
#
# Compute the spectral ordering of a graph.
function spectralorder(matrix::SparseMatrixCSC; tol=0.0)
    matrix = SparseMatrixCSC{Float64, Int}(matrix)
    fill!(nonzeros(matrix), 1)
    fkeep!((i, j, v) -> i != j, matrix)
    value, vector = Laplacians.fiedler(matrix; tol)
    sortperm(reshape(vector, size(matrix, 2)))
end


end
