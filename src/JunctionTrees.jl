module JunctionTrees


import AMD
import CuthillMcKee
import Metis

using AbstractTrees
using Catlab.BasicGraphs
using DataStructures
using SparseArrays


# Orders
export Order


# Elimination Algorithms
export AMDJL_AMD, CuthillMcKeeJL_RCM, MetisJL_ND


# Supernode Types
export Node, Maximal, Fundamental 


# Supernode Trees
export SupernodeTree, width, height, seperators, supernode, order, inverse


include("junction_trees/orders.jl")
include("junction_trees/elimination_algorithms.jl")
include("junction_trees/ordered_graphs.jl")
include("junction_trees/trees.jl")
include("junction_trees/postorder_trees.jl")
include("junction_trees/elimination_trees.jl")
include("junction_trees/supernode_types.jl")
include("junction_trees/supernode_trees.jl")


end
