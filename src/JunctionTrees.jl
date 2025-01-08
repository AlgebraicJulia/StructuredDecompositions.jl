module JunctionTrees


using AbstractTrees
using Base: DEFAULT_STABLE, OneTo
using Base.Order
using DataStructures: IntDisjointSets, find_root!, root_union!
using LinearAlgebra
using SparseArrays
using SparseArrays: getcolptr


import AMD
import Base.Iterators
import CuthillMcKee
import Metis
import SymRCM
import TreeWidthSolver


const AbstractScalar{T} = AbstractArray{T, 0}
const Scalar{T} = Array{T, 0}


# Elimination Algorithms
export AMDJL_AMD, AMDJL_SymAMD, CuthillMcKeeJL_RCM, SymRCMJL_RCM, MetisJL_ND, TreeWidthSolverJL_BT, MCS, permutation


# Ordered Graphs
export OrderedGraph


# Supernode Types
export Nodal, Maximal, Fundamental 


# Bags
export Bag, separator, residual


# Junction Trees
export JunctionTree, junctiontree, junctiontree!, treewidth, relative


include("junction_trees/disjoint_trees.jl")
include("junction_trees/disjoint_lists.jl")
include("junction_trees/utils.jl")
include("junction_trees/elimination_algorithms.jl")
include("junction_trees/supernode_types.jl")
include("junction_trees/children.jl")
include("junction_trees/trees.jl")
include("junction_trees/bags.jl")
include("junction_trees/junction_trees.jl")


end
