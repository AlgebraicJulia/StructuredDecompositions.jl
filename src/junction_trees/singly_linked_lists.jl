"""
    SinglyLinkedList{Init <: AbstractScalar{Int}, Next <: AbstractVector{Int}}

A [singly linked list](https://en.wikipedia.org/wiki/Linked_list).
This type supports the [iteration interface](https://docs.julialang.org/en/v1/manual/interfaces/).
"""
struct SinglyLinkedList{Init <: AbstractScalar{Int}, Next <: AbstractVector{Int}}
    head::Init
    next::Next
end


function Base.isempty(list::SinglyLinkedList)
    iszero(list.head[])
end


function Base.pushfirst!(list::SinglyLinkedList, v::Integer)
    list.next[v] = list.head[]
    list.head[] = v
    list
end


function Base.popfirst!(list::SinglyLinkedList)
    v = list.head[]
    list.head[] = list.next[v]
    v
end


function Base.show(io::IO, list::T) where T <: SinglyLinkedList
    items = pushfirst!(map(string, take(list, MAX_ITEMS_PRINTED)), "head")

    if MAX_ITEMS_PRINTED < length(items)
        items[end] = "..."
    end

    println(io, T)
    println(io, join(items, " → "))
end


######################
# Iterator Interface #
######################


function Base.iterate(list::SinglyLinkedList, i::Integer=list.head[])
    iszero(i) ? nothing : (i, list.next[i])
end


function Base.IteratorSize(::Type{T}) where T <: SinglyLinkedList
    Base.SizeUnknown()
end


function Base.eltype(::Type{T}) where T <: SinglyLinkedList
    Int
end

