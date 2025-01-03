# The children of a vertex in a rooted tree.
struct TreeChildren{T <: AbstractTree}
    tree::T
    index::Int
end


function Base.show(io::IO, children::TreeChildren)
    i = children.index
    print(io, "TreeChildren $i")
end


######################
# Iterator Interface #
######################


function Base.iterate(children::TreeChildren, i::Int=firstchildindex(children.tree, children.index))
    i, nextsiblingindex(children.tree, i)
end


function Base.iterate(children::TreeChildren, ::Nothing) end


function Base.IteratorSize(::Type{TreeChildren{T}}) where T <: AbstractTree
    Base.SizeUnknown()
end


function Base.eltype(::Type{TreeChildren{T}}) where T <: AbstractTree
    Int
end


function Base.hasfastin(::Type{TreeChildren{T}}) where T <: AbstractTree
    true
end


function Base.in(i, children::TreeChildren)
    children.index == parentindex(children.tree, i)
end


##########################
# Indexed Tree Interface #
##########################


function AbstractTrees.childindices(tree::AbstractTree, i::Integer)
    TreeChildren(tree, i)
end
