# The children of a node in a rooted tree.
struct Children
    child::Int
    brother::Vector{Int}
end


function Base.show(io::IO, children::Children)
    i = children.index
    print(io, "Children $i")
end


######################
# Iterator Interface #
######################


function Base.iterate(children::Children, i::Int=children.child)
    iszero(i) ? nothing : (i, children.brother[i])
end


function Base.IteratorSize(::Type{Children})
    Base.SizeUnknown()
end


function Base.eltype(::Type{Children})
    Int
end


