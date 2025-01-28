# A doubly linked list of distinct integers.
struct SinglyLinkedList{I <: Integer, Init <: AbstractScalar{I}, Next <: AbstractVector{I}}
    head::Init
    next::Next
end


# Evaluate whether a linked list is empty.
function Base.isempty(list::SinglyLinkedList)
    iszero(list.head[])
end


# Append an element `v` to the front of a linked list.
# If `v` ∈ `list`, the behavior of this function is undefined.
function Base.pushfirst!(list::SinglyLinkedList, v::Integer)
    list.next[v] = list.head[]
    list.head[] = v
    list
end


# Remove the first element of a linked list.
function Base.popfirst!(list::SinglyLinkedList)
    v = list.head[]
    list.head[] = list.next[v]
    v
end


function Base.show(io::IO, ::MIME"text/plain", list::L) where L <: SinglyLinkedList
    items = pushfirst!(map(string, take(list, MAX_ITEMS_PRINTED)), "head")

    if MAX_ITEMS_PRINTED < length(items)
        items[end] = "..."
    end

    println(io, L)
    println(io, join(items, " → "))
end


#######################
# Iteration Interface #
#######################


function Base.iterate(list::SinglyLinkedList{I}, i::I=list.head[]) where I
    iszero(i) ? nothing : (i, list.next[i])
end


function Base.IteratorSize(::Type{<:SinglyLinkedList})
    Base.SizeUnknown()
end


function Base.eltype(::Type{<:SinglyLinkedList{I}}) where I
    I
end

