module Decompositions

export StructuredDecomposition, StrDecomp, 
      DecompType, Decomposition, CoDecomposition, 
      𝐃, bags, adhesions, adhesionSpans, 
      ∫

using AbstractTrees
using CliqueTrees
using CliqueTrees: PermutationOrAlgorithm, EliminationAlgorithm, SupernodeType, DEFAULT_ELIMINATION_ALGORITHM, DEFAULT_SUPERNODE_TYPE

using PartialFunctions
using MLStyle

using Catlab
using Catlab.CategoricalAlgebra
using Catlab.Graphs
using Catlab.Graphs: add_edges!
using Catlab.ACSetInterface
using Catlab.CategoricalAlgebra.Diagrams

import Catlab.CategoricalAlgebra.Diagrams: ob_map, hom_map, colimit, limit

include("str_decomps.jl")
include("clique_trees.jl")

end
