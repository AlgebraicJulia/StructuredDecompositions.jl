# The parent of a vertex in a rooted tree.
struct TreeParent{T <: AbstractTree}
    tree::T
    index::Int
end


function Base.show(io::IO, parent::TreeParent)
    i = parent.index
    print(io, "TreeParent $i")
end


######################
# Iterator Interface #
######################


function Base.iterate(parent::TreeParent, i::Int=parentindex(parent.tree, parent.index))
    i, nothing
end


function Base.iterate(parent::TreeParent, ::Nothing) end


function Base.length(parent::TreeParent)
    parent.index == rootindex(parent.tree) ? 0 : 1
end


function Base.eltype(::Type{TreeParent{T}}) where T <: AbstractTree
    Int
end


function Base.hasfastin(::Type{TreeParent{T}}) where T <: AbstractTree
    true
end


function Base.in(i, parent::TreeParent)
    i == parentindex(parent.tree, parent.index)
end


##########################
# Indexed Tree Interface #
##########################


function Graphs.inneighbors(tree::AbstractTree, i::Integer)
    TreeParent(tree, i)
end
