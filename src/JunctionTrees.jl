module JunctionTrees


using AbstractTrees
using Base: DEFAULT_STABLE, OneTo
using Base.Order
using DataStructures: IntDisjointSets, find_root!, root_union!
using LinearAlgebra
using SparseArrays


import AMD
import Base.Iterators
import CuthillMcKee
import Metis
import SymRCM
import TreeWidthSolver


# Permutations
export Permutation


# Elimination Algorithms
export AMDJL_AMD, AMDJL_SYMAMD, CuthillMcKeeJL_RCM, SymRCMJL_RCM, MetisJL_ND, TreeWidthSolverJL_BT, MCS, permutation


# Ordered Graphs
export OrderedGraph


# Supernode Types
export Node, Maximal, Fundamental 


# Junction Trees
export JunctionTree, junctiontree, junctiontree!, treewidth, separator, residual, relative


include("junction_trees/stacks.jl")
include("junction_trees/indices.jl")
include("junction_trees/disjoint_trees.jl")
include("junction_trees/disjoint_lists.jl")
include("junction_trees/sparse.jl")
include("junction_trees/elimination_algorithms.jl")
include("junction_trees/supernode_types.jl")
include("junction_trees/permutations.jl")
include("junction_trees/children.jl")
include("junction_trees/trees.jl")
include("junction_trees/bags.jl")
include("junction_trees/junction_trees.jl")


end
