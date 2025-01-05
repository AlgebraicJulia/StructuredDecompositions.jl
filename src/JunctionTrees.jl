module JunctionTrees


using AbstractTrees
using Base: DEFAULT_STABLE, OneTo
using Base.Order
using DataStructures: IntDisjointSets, find_root!, root_union!
using LinearAlgebra
using LinkedLists: ListNode, LinkedList
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
export AMDJL_AMD, AMDJL_SYMAMD, CuthillMcKeeJL_RCM, SymRCMJL_RCM, MetisJL_ND, TreeWidthSolverJL_BT, MCS


# Ordered Graphs
export OrderedGraph


# Supernode Types
export Node, Maximal, Fundamental 


# Junction Trees
export JunctionTree, junctiontree, treewidth, seperator, residual, clique, relative


include("junction_trees/sparse.jl")
include("junction_trees/elimination_algorithms.jl")
include("junction_trees/supernode_types.jl")
include("junction_trees/stacks.jl")
include("junction_trees/disjoint_sets.jl")
include("junction_trees/permutations.jl")
include("junction_trees/children.jl")
include("junction_trees/trees.jl")
include("junction_trees/bags.jl")
include("junction_trees/junction_trees.jl")


end
