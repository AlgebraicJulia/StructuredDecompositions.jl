using BenchmarkTools
using Catlab
using LinearAlgebra
using MatrixMarket
using StructuredDecompositions.Decompositions
using StructuredDecompositions.Decompositions: adjacency_matrix
using StructuredDecompositions.DecidingSheaves
using StructuredDecompositions.FunctorUtils
using SuiteSparseMatrixCollection
using SparseArrays


import SimpleGraphs
import SimpleGraphAlgorithms
import QDLDL
import GenericTensorNetworks
import Graphs as GraphsPkg


const SUITE = BenchmarkGroup()


# fixing bug upstream
function Catlab.WiringDiagramAlgebras.make_homomorphism(row, X::StructACSet{S}, Y::StructACSet{S}) where S
  components = let i = 0
    NamedTuple{ob(S)}(Int[row[i+=1] for _ in parts(X,c)] for c in ob(S))
  end
  ACSetTransformation(components, X, Y)
end


include("graph_coloring_fixed.jl")
