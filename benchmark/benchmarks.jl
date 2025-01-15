using BenchmarkTools
using Catlab.Graphs
using Catlab.CategoricalAlgebra
using LinearAlgebra
using MatrixMarket
using StructuredDecompositions.Decompositions
using StructuredDecompositions.Decompositions: adjacency_matrix
using StructuredDecompositions.DecidingSheaves
using StructuredDecompositions.FunctorUtils
using StructuredDecompositions.JunctionTrees
using SuiteSparseMatrixCollection
using SimpleGraphs
using SimpleGraphAlgorithms
using SparseArrays
using QDLDL

const SUITE = BenchmarkGroup()

include("junction_trees.jl")
include("graph_coloring_fixed.jl")
