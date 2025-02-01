"""
    AbstractTree{I} = Union{Tree{I}, SupernodeTree{I}, JunctionTree{I}}

A rooted forest.
This type implements the [indexed tree interface](https://juliacollections.github.io/AbstractTrees.jl/stable/#The-Indexed-Tree-Interface).
"""
const AbstractTree{I} = Union{Tree{I}, SupernodeTree{I}, JunctionTree{I}}


function Base.show(io::IO, ::MIME"text/plain", tree::T) where T <: AbstractTree
    println(io, "$(length(tree))-element $T:")

    for i in rootindices(tree)
        print_tree(io, IndexNode(tree, i))
    end
end


##########################
# Indexed Tree Interface #
##########################


"""
    rootindex(tree::AbstractTree)

Get the first root of a rooted forest.
"""
AbstractTrees.rootindex(tree::AbstractTree)

"""
    rootindices(tree::AbstractTree)

Construct an iterator over the roots of a rooted forest.
"""
rootindices(tree::AbstractTree)


"""
    parentindex(tree::AbstractTree, i::Integer)

Get the parent index of node `i`. Returns `nothing` if `i` is a root.
"""
AbstractTrees.parentindex(tree::AbstractTree, i::Integer)


"""
    firstchildindex(tree::AbstractTree, i::Integer)

Get the first child of node `i`. Returns `nothing` if `i` is a leaf.
"""
firstchildindex(tree::AbstractTree, i::Integer)


"""
    nextsiblingindex(tree::AbstractTree, i::Integer)

Get the next sibling of node `i`. Returns `nothing` if `i` is the greatest sibling.
"""
AbstractTrees.nextsiblingindex(tree::AbstractTree, i::Integer)


"""
    childindices(tree::AbstractTree, i::Integer)

Construct an iterator over the children of node `i`.
"""
AbstractTrees.childindices(tree::AbstractTree, i::Integer)


"""
    ancestorindices(tree::AbstractTree, i::Integer)

Construct an iterator over the ancestors of node `i`.
"""
ancestorindices(tree::AbstractTree, i::Integer)


function AbstractTrees.ParentLinks(::Type{IndexNode{T, I}}) where {I, T <: AbstractTree{I}}
    StoredParents()
end


function AbstractTrees.SiblingLinks(::Type{IndexNode{T, I}}) where {I, T <: AbstractTree{I}}
    StoredSiblings()
end


function AbstractTrees.NodeType(::Type{IndexNode{T, I}}) where {I, T <: AbstractTree{I}}
    HasNodeType()
end


function AbstractTrees.nodetype(::Type{IndexNode{T, I}}) where {I, T <: AbstractTree{I}}
    IndexNode{T, I}
end
