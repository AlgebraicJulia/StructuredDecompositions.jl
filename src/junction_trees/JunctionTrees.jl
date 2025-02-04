module JunctionTrees


using AbstractTrees
using Base: OneTo
using Base.Order
using Base.Iterators: take, takewhile
using DataStructures: IntDisjointSets, find_root!, root_union!
using LinearAlgebra
using SparseArrays
using SparseArrays: getcolptr, indtype
using Sparspak: SpkMmd


import AMD as AMDPkg


const AbstractScalar{T} = AbstractArray{T, 0}
const Scalar{T} = Array{T, 0}
const MAX_ITEMS_PRINTED = 5


# Linked Lists
export SinglyLinkedList


# Elimination Algorithms
export MCS, RCM, AMD, SymAMD, MMD, NodeND, Spectral, BT, permutation


# Trees
export Tree, eliminationtree, rootindex, firstchildindex, nextsiblingindex, parentindex, rootindices, childindices, ancestorindices, setrootindex!


# Supernode Types
export Nodal, Maximal, Fundamental 


# Bags
export Bag, separator, residual


# Supernode Trees
export SupernodeTree, supernodetree


# Junction Trees
export JunctionTree, junctiontree, treewidth, relative, relative!


# Chordal Graphs
export eliminationgraph, eliminationgraph!, ischordal, isfilled, isperfect


include("singly_linked_lists.jl")
include("linked_lists.jl")
include("graphs.jl")
include("elimination_algorithms.jl")
include("trees.jl")
include("supernode_types.jl")
include("bags.jl")
include("supernode_trees.jl")
include("junction_trees.jl")
include("abstract_trees.jl")
include("chordal_graphs.jl")
include("utils.jl")


end
