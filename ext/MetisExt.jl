module MetisExt


import Metis


using SparseArrays
using StructuredDecompositions.JunctionTrees


function JunctionTrees.permutation(matrix::SparseMatrixCSC, alg::NodeND)
    order::Vector{Int}, index::Vector{Int} = Metis.permutation(matrix)
    order, index
end


end
