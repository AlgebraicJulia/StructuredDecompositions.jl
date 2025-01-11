"""
    SupernodeType

A type of supernode partition. The options are

| type                  | name                            |
| :-------------------- | :------------------------------ |
| [`Nodal`](@ref)       | nodal supernode partition       |
| [`Maximal`](@ref)     | maximal supernode partition     |
| [`Fundamental`](@ref) | fundamental supernode partition |
"""
abstract type SupernodeType end


"""
    Nodal <: SupernodeType

A nodal  supernode partition.
"""
struct Nodal <: SupernodeType end


"""
    Maximal <: SupernodeType

A maximal supernode partition.
"""
struct Maximal <: SupernodeType end


"""
    Fundamental <: SupernodeType

A fundamental supernode partition.
"""
struct Fundamental <: SupernodeType end


const DEFAULT_SUPERNODE_TYPE = Maximal()
