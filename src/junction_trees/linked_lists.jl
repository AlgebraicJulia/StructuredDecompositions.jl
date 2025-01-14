# A doubly linked list.
struct LinkedList{Init <: AbstractScalar{Int}, Prev <: AbstractVector{Int}, Next <: AbstractVector{Int}}
    head::Init
    prev::Prev
    next::Next
end


function Base.isempty(list::LinkedList)
    iszero(list.head[])
end


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


function Base.popfirst!(list::LinkedList)
    v = list.head[]
    n = list.next[v]
    list.head[] = n

    if !iszero(n)
        list.prev[n] = 0
    end

    v
end


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


# Simple Linear-Time Algorithms to Test Chordality of Graphs, Test Acyclicity of Hypergraphs, and Selectively Reduce Acyclic Hypergraphs
# Tarjan and Yannakakis
# Maximum Cardinality Search
#
# Construct a fill-reducing permutation of a simple graph, represented by its adjacency matrix.
# The complexity is O(m + n), where m = |E| and n = |V|.
function mcs(matrix::SparseMatrixCSC)
    # validate arguments
    size(matrix, 1) != size(matrix, 2) && throw(ArgumentError("size(matrix, 1) != size(matrix, 2)"))
    
    # construct disjoint sets data structure
    head = zeros(Int, size(matrix, 2) + 1)
    prev = Vector{Int}(undef, size(matrix, 2) + 1)
    next = Vector{Int}(undef, size(matrix, 2) + 1)
    
    function set(i)
        LinkedList(view(head, i), prev, next)
    end
    
    # run algorithm
    order = Vector{Int}(undef, size(matrix, 2))
    number = Vector{Int}(undef, size(matrix, 2))

    for v in axes(matrix, 2)
        number[v] = 1
        pushfirst!(set(1), v)
    end

    i = size(matrix, 2)
    j = 1

    while i >= 1
        v = popfirst!(set(j))
        order[i] = v
        number[v] = 0

        for w in @view rowvals(matrix)[nzrange(matrix, v)]
            if number[w] >= 1
                delete!(set(number[w]), w)
                number[w] += 1
                pushfirst!(set(number[w]), w)
            end
        end

        i -= 1
        j += 1

        while j >= 1 && isempty(set(j))
            j -= 1
        end
    end

    order
end


function Base.show(io::IO, list::T) where T <: LinkedList
    items = pushfirst!(map(string, take(list, MAX_ITEMS_PRINTED)), "head")

    if MAX_ITEMS_PRINTED < length(items)
        items[end] = "..."
    end

    print(io, "$T:\n$(join(items, " ↔ "))")
end


######################
# Iterator Interface #
######################


function Base.iterate(list::LinkedList, i::Integer=list.head[])
    iszero(i) ? nothing : (i, list.next[i])
end


function Base.IteratorSize(::Type{T}) where T <: LinkedList
    Base.SizeUnknown()
end


function Base.eltype(::Type{T}) where T <: LinkedList
    Int
end
