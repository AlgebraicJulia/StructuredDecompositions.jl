struct ChildIndices
    tree::Tree
    index::Int
end


function Base.iterate(iterator::ChildIndices)
    iterate(iterator, iterator.tree.child[iterator.index])
end


function Base.iterate(iterator::ChildIndices, i::Integer)
    if iszero(i)
        nothing
    else
        i, iterator.tree.brother[i]
    end
end


function Base.IteratorSize(::Type{ChildIndices})
    Base.SizeUnknown()
end


function Base.eltype(::Type{ChildIndices})
    Int
end


function AbstractTrees.childindices(tree::Tree, i::Integer)
    ChildIndices(tree, i)
end
