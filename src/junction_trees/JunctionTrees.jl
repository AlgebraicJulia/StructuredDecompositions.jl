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
export SupernodeTree, width, height


include("orders.jl")
include("elimination_algorithms.jl")
include("ordered_graphs.jl")
include("abstract_trees.jl")
include("trees.jl")
include("postorder_trees.jl")
include("elimination_trees.jl")
include("supernode_types.jl")
include("supernode_trees.jl")


end
