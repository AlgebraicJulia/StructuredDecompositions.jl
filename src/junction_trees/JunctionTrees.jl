module JunctionTrees


using AbstractTrees
using Base.Order
using Base.Iterators: drop, filter as ifilter, map as imap, reverse as ireverse, peel, take
using DataStructures: IntDisjointSets, find_root!, root_union!
using LinearAlgebra
using SparseArrays
using SparseArrays: getcolptr
using Sparspak: SpkMmd


import AMD as AMDPkg
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
export MCS, RCM, AMD, SymAMD, MMD, NodeND, Spectral, BT, permutation


# Trees
export Tree, eliminationtree, eliminationtree!, rootindex, firstchildindex, nextsiblingindex, parentindex, rootindices, childindices


# Supernode Types
export Nodal, Maximal, Fundamental 


# Bags
export Bag, separator, residual


# Supernode Trees
export SupernodeTree, supernodetree, supernodetree!


# Junction Trees
export JunctionTree, junctiontree, junctiontree!, treewidth, treewidth!, relative


# Chordal Graphs
export ischordal, isperfect, filledgraph


include("singly_linked_lists.jl")
include("linked_lists.jl")
include("elimination_algorithms.jl")
include("trees.jl")
include("supernode_types.jl")
include("bags.jl")
include("supernode_trees.jl")
include("junction_trees.jl")
include("chordal_graphs.jl")
include("utils.jl")

end
