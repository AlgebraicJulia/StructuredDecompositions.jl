module JunctionTrees


using AbstractTrees
using Base.Order: Ordering
using DataStructures
using LinearAlgebra
using SparseArrays
using SparseArrays: AbstractSparseMatrixCSC


import AMD
import Base.Iterators
import Catlab
import CuthillMcKee
import Graphs
import LinkedLists: ListNode, LinkedList
import Metis
import SymRCM
import TreeWidthSolver


# Permutations
export Permutation


# Elimination Algorithms
export AMDJL_AMD, AMDJL_SYMAMD, CuthillMcKeeJL_RCM, SymRCMJL_RCM, MetisJL_ND, TreeWidthSolverJL_BT, MCS, DFS


# Ordered Graphs
export OrderedGraph


# Supernode Types
export Node, Maximal, Fundamental 


# Junction Trees
export JunctionTree, jtree, treewidth, seperator, residual, clique, relative


include("junction_trees/elimination_algorithms.jl")
include("junction_trees/supernode_types.jl")
include("junction_trees/fixed_stacks.jl")
include("junction_trees/disjoint_sets.jl")
include("junction_trees/sum_vectors.jl")
include("junction_trees/permutations.jl")
include("junction_trees/trees/abstract_trees.jl")
include("junction_trees/trees/trees.jl")
include("junction_trees/trees/tree_edges.jl")
include("junction_trees/trees/tree_children.jl")
include("junction_trees/trees/tree_parent.jl")
include("junction_trees/graphs/abstract_graphs.jl")
include("junction_trees/graphs/ordered_graphs.jl")
include("junction_trees/graphs/ordered_graph_edges.jl")
include("junction_trees/junction_trees.jl")
include("junction_trees/utilities.jl")


end
