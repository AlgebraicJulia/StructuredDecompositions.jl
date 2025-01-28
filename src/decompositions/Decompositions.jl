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
using Catlab.Graphs: add_edges!
using Catlab.ACSetInterface
using Catlab.CategoricalAlgebra.Diagrams
using SparseArrays: getcolptr
using SparseArrays

import Catlab.CategoricalAlgebra.Diagrams: ob_map, hom_map, colimit, limit

include("str_decomps.jl")
include("junction_trees.jl")

end
