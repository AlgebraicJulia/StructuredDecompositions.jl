module Decompositions

export StructuredDecomposition, StrDecomp, 
      DecompType, Decomposition, CoDecomposition, 
      𝐃, bags, adhesions, adhesionSpans, 
      ∫

using ..JunctionTrees
using ..JunctionTrees: PermutationOrAlgorithm, EliminationAlgorithm, SupernodeType, DEFAULT_ELIMINATION_ALGORITHM, DEFAULT_SUPERNODE_TYPE

using PartialFunctions
using MLStyle

using Base: DEFAULT_STABLE, ForwardOrdering
using Base.Threads
using Catlab
using Catlab.CategoricalAlgebra
using Catlab.Graphs
using Catlab.ACSetInterface
using Catlab.CategoricalAlgebra.Diagrams
using SparseArrays

import Catlab.CategoricalAlgebra.Diagrams: ob_map, hom_map, colimit, limit

include("decompositions/str_decomps.jl")
include("decompositions/junction_trees.jl")

end
