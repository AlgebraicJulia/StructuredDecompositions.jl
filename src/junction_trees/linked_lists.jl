# A doubly linked list of distinct integers.
struct LinkedList{I <: Integer, Init <: AbstractScalar{I}, Prev <: AbstractVector{I}, Next <: AbstractVector{I}}
    head::Init
    prev::Prev
    next::Next
end


# Evaluate whether a linked list is empty.
function Base.isempty(list::LinkedList)
    iszero(list.head[])
end


# Append an element `v` to the front of a linked list.
# If `v` ∈ `list`, the behavior of this function is undefined.
function Base.pushfirst!(list::LinkedList, v::Integer)
    n = list.head[]
    list.head[] = v
    list.prev[v] = 0
    list.next[v] = n

    if !iszero(n)
        list.prev[n] = v
    end

    list
end


# Remove the first element of a linked list.
function Base.popfirst!(list::LinkedList)
    v = list.head[]
    n = list.next[v]
    list.head[] = n

    if !iszero(n)
        list.prev[n] = 0
    end

    v
end


# Delete an element `v` from a linked list.
# If `v` ∉ `list`, the behavior of this function is undefined.
function Base.delete!(list::LinkedList, v::Integer)
    p = list.prev[v]
    n = list.next[v]

    if !iszero(p)
        list.next[p] = n
    else
        list.head[] = n
    end

    if !iszero(n)
        list.prev[n] = p
    end

    list
end


function Base.show(io::IO, ::MIME"text/plain", list::L) where L <: LinkedList
    items = pushfirst!(map(string, take(list, MAX_ITEMS_PRINTED)), "head")

    if MAX_ITEMS_PRINTED < length(items)
        items[end] = "..."
    end

    println(io, L)
    println(io, join(items, " ↔ "))
end


#######################
# Iteration Interface #
#######################


function Base.iterate(list::LinkedList{I}, i::I=list.head[]) where I
    iszero(i) ? nothing : (i, list.next[i])
end


function Base.IteratorSize(::Type{<:LinkedList})
    Base.SizeUnknown()
end


function Base.eltype(::Type{<:LinkedList{I}}) where I
    I
end
