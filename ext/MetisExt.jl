module MetisExt


using Metis
using SparseArrays
using StructuredDecompositions.JunctionTrees


function JunctionTrees.permutation(matrix::SparseMatrixCSC{T, I}, alg::NodeND) where {T, I}
    order::Vector{I}, index::Vector{I} = Metis.permutation(matrix)
    order, index
end


end
