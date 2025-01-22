module MetisExt


import Metis


using SparseArrays
using StructuredDecompositions.JunctionTrees


function JunctionTrees.permutation(matrix::SparseMatrixCSC, alg::NodeND)
    order, index = Metis.permutation(matrix)
    order = convert(Vector{Int}, order)
    index = convert(Vector{Int}, index)
    order, index
end


end
