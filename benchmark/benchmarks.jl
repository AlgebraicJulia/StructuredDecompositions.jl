import Graphs as GraphsPkg

using BenchmarkTools
using Catlab.Graphs
using Catlab.CategoricalAlgebra
using MatrixMarket
using StructuredDecompositions.Decompositions
using StructuredDecompositions.DecidingSheaves
using StructuredDecompositions.FunctorUtils
using StructuredDecompositions.JunctionTrees
using SuiteSparseMatrixCollection

const SUITE = BenchmarkGroup()

include("JunctionTrees.jl")
include("GraphColoring.jl")
