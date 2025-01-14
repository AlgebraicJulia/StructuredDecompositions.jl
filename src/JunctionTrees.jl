module JunctionTrees


using AbstractTrees
using Base.Order
using Base.Iterators: drop, filter as ifilter, map as imap, reverse as ireverse, take
using DataStructures: IntDisjointSets, find_root!, root_union!
using LinearAlgebra
using SparseArrays
using SparseArrays: getcolptr
using Sparspak: SpkMmd


import AMD as AMDPkg
import FlowCutterPACE17_jll
import Laplacians
import Metis
import SymRCM
import TreeWidthSolver


const AbstractScalar{T} = AbstractArray{T, 0}
const Scalar{T} = Array{T, 0}
const MAX_ITEMS_PRINTED = 5


# Linked Lists
export SinglyLinkedList


# Elimination Algorithms
export MCS, RCM, AMD, SymAMD, MMD, NodeND, FlowCutter, Spectral, BT, permutation


# Trees
export Tree, eliminationtree, eliminationtree!, firstchildindex, rootindex, nextsiblingindex, parentindex, childindices


# Supernode Types
export Nodal, Maximal, Fundamental 


# Bags
export Bag, separator, residual


# Supernode Trees
export SupernodeTree, supernodetree, supernodetree!


# Junction Trees
export JunctionTree, junctiontree, junctiontree!, treewidth, treewidth!, relative


# Chordal Graphs
export ischordal, isperfect, chordalgraph


include("junction_trees/singly_linked_lists.jl")
include("junction_trees/linked_lists.jl")
include("junction_trees/elimination_algorithms.jl")
include("junction_trees/children.jl")
include("junction_trees/trees.jl")
include("junction_trees/supernode_types.jl")
include("junction_trees/bags.jl")
include("junction_trees/supernode_trees.jl")
include("junction_trees/junction_trees.jl")
include("junction_trees/chordal_graphs.jl")
include("junction_trees/utils.jl")

end
