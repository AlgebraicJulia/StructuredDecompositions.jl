module JunctionTrees


using AbstractTrees
using Base: DEFAULT_STABLE, OneTo
using Base.Order
using Base.Iterators: drop, filter as ifilter, map as imap, reverse as ireverse
using DataStructures: IntDisjointSets, find_root!, root_union!
using LinearAlgebra
using SparseArrays
using SparseArrays: getcolptr
using Sparspak: SpkMmd


import AMD as AMDPkg
import ArnoldiMethod
import FlowCutterPACE17_jll
import Metis
import SymRCM
import TreeWidthSolver


const AbstractScalar{T} = AbstractArray{T, 0}
const Scalar{T} = Array{T, 0}


# Elimination Algorithms
export RCM, MMD, AMD, SymAMD, NodeND, BT, MCS, FlowCutter, Spectral, ischordal, isperfect, permutation


# Ordered Graphs
export OrderedGraph


# Supernode Types
export Nodal, Maximal, Fundamental 


# Bags
export Bag, separator, residual


# Junction Trees
export JunctionTree, junctiontree, junctiontree!, treewidth, treewidth!, chordalgraph, relative


include("junction_trees/disjoint_trees.jl")
include("junction_trees/disjoint_lists.jl")
include("junction_trees/utils.jl")
include("junction_trees/elimination_algorithms.jl")
include("junction_trees/supernode_types.jl")
include("junction_trees/children.jl")
include("junction_trees/trees.jl")
include("junction_trees/bags.jl")
include("junction_trees/junction_trees.jl")
include("junction_trees/chordal_graphs.jl")


end
