module JunctionTrees


import AMD
import Catlab.BasicGraphs
import CuthillMcKee
import LinkedLists
import Metis
import TreeWidthSolver

using AbstractTrees
using Base.Order: Ordering
using DataStructures
using Graphs
using Graphs.SimpleGraphs
using LinearAlgebra
using SparseArrays
using SparseArrays: AbstractSparseMatrixCSC


# Orders
export Order


# Elimination Algorithms
export AMDJL_AMD, CuthillMcKeeJL_RCM, MetisJL_ND, TreeWidthSolverJL_BT, MCS


# Ordered Graphs
export OrderedGraph


# Supernode Types
export Node, Maximal, Fundamental 


# Junction Trees
export JunctionTree, treewidth, seperator, residual, clique, find_clique, lift_par, lift_sep, lift


include("junction_trees/fixed_stacks.jl")
include("junction_trees/disjoint_sets.jl")
include("junction_trees/orders.jl")
include("junction_trees/elimination_algorithms.jl")
include("junction_trees/trees.jl")
include("junction_trees/child_indices.jl")
include("junction_trees/ordered_graphs.jl")
include("junction_trees/supernode_types.jl")
include("junction_trees/junction_trees.jl")


end
