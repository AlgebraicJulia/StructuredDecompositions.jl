module JunctionTrees


import AMD
import CuthillMcKee
import LinkedLists
import Metis
import TreeWidthSolver

using AbstractTrees
using Catlab.BasicGraphs
using DataStructures
using SparseArrays
using SparseArrays: AbstractSparseMatrixCSC


# Orders
export Order, inverse


# Elimination Algorithms
export AMDJL_AMD, CuthillMcKeeJL_RCM, MetisJL_ND, TreeWidthSolverJL_BT, MCS


# Supernode Types
export Node, Maximal, Fundamental 


# Junction Trees
export JunctionTree, width, height, seperator, residual, clique, seperator_to_parent, seperator_to_self


include("junction_trees/orders.jl")
include("junction_trees/elimination_algorithms.jl")
include("junction_trees/ordered_graphs.jl")
include("junction_trees/trees.jl")
include("junction_trees/postorder_trees.jl")
include("junction_trees/elimination_trees.jl")
include("junction_trees/supernode_types.jl")
include("junction_trees/supernode_trees.jl")
include("junction_trees/junction_trees.jl")


end
