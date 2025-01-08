"""
    SupernodeType

A type of supernode. The options are
- [`Nodal`](@ref)
- [`Maximal`](@ref)
- [`Fundamental`](@ref)
"""
abstract type SupernodeType end


"""
    Nodal <: Supernode

A single-vertex supernode.
"""
struct Nodal <: SupernodeType end


"""
    Maximal <: Supernode

A maximal supernode.
"""
struct Maximal <: SupernodeType end


"""
    Fundamental <: Supernode

A fundamental supernode.
"""
struct Fundamental <: SupernodeType end


const DEFAULT_SUPERNODE_TYPE = Maximal()
