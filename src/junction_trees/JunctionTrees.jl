module JunctionTrees


import AMD
import CuthillMcKee
import Metis

using AbstractTrees
using Catlab.BasicGraphs
using DataStructures
using SparseArrays

# Elimination Algorithms
export EliminationAlgorithm, AMDJL_AMD, CuthillMcKeeJL_RCM, MetisJL_ND, MCS

# Supernodes
export Supernode, Node, MaximalSupernode, FundamentalSupernode

# Orders
export Order

# Elimination Trees
export EliminationTree
export getwidth, getsupernode, getsubtree, getlevel

# Junction Trees
export JunctionTree
export getseperator, getresidual


# Add an element x to a sorted set v.
# Returns true if x ∉ v.
# Returns false if x ∈ v.
function insertsorted!(v::Vector, x::Integer)
    i = searchsortedfirst(v, x)

    if i > length(v) || v[i] != x
        insert!(v, i, x)
        true
    else
        false
    end
end


# Delete an element x from a sorted set v.
# Returns true if x ∈ v.
# Returns false if x ∉ v.
function deletesorted!(v::Vector, x::Integer)
    i = searchsortedfirst(v, x)

    if i <= length(v) && v[i] == x
        deleteat!(v, i)
        true
    else
        false
    end
end


# Delete the elements xs from a sorted set v.
# Returns true if xs and v intersect.
# Returns false if xs and v are disjoint.
function deletesorted!(v::Vector, xs::AbstractVector)
    isintersecting = true

    for x in xs
        isintersecting = deletesorted!(v, x) || isintersecting
    end

    isintersecting
end


include("elimination_algorithms.jl")
include("supernodes.jl")
include("orders.jl")
include("trees.jl")
include("elimination_trees.jl")
include("junction_trees.jl")


end
