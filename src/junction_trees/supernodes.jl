"""
    Supernode

A type of supernode. The options are
- [`Node`](@ref)
- [`MaximalSupernode`](@ref)
- [`FundamentalSupernode`](@ref)
"""
abstract type Supernode end


"""
    Node <: Supernode

A single-vertex supernode.
"""
struct Node <: Supernode end


"""
    MaximalSupernode <: Supernode

A maximal supernode.
"""
struct MaximalSupernode <: Supernode end


"""
    FundamentalSupernode <: Supernode

A fundamental supernode.
"""
struct FundamentalSupernode <: Supernode end


const DEFAULT_SUPERNODE = MaximalSupernode()
