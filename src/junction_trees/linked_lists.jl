# A doubly linked list of distinct integers.
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


function mcs(matrix::SparseMatrixCSC)
    # validate arguments
    size(matrix, 1) != size(matrix, 2) && throw(ArgumentError("size(matrix, 1) != size(matrix, 2)"))

    # run algorithm
    mcs(size(matrix, 2)) do j
        @view rowvals(matrix)[nzrange(matrix, j)]
    end
end


# Simple Linear-Time Algorithms to Test Chordality of Graphs, Test Acyclicity of Hypergraphs, and Selectively Reduce Acyclic Hypergraphs
# Tarjan and Yannakakis
# Maximum Cardinality Search
#
# Construct a fill-reducing permutation of a simple graph.
# The complexity is O(m + n), where m = |E| and n = |V|.
function mcs(neighbors::Function, n::Integer) 
    # validate arguments
    n < 0 && throw(ArgumentError("n < 0"))

    # construct disjoint sets data structure
    head = zeros(Int, n + 1)
    prev = Vector{Int}(undef, n + 1)
    next = Vector{Int}(undef, n + 1)
    
    function set(i)
        LinkedList(view(head, i), prev, next)
    end
    
    # run algorithm
    alpha = Vector{Int}(undef, n)
    size = Vector{Int}(undef, n)

    for v in n:-1:1
        size[v] = 1
        pushfirst!(set(1), v)
    end

    j = 1

    for i in n:-1:1
        v = popfirst!(set(j))
        alpha[v] = i
        size[v] = 1 - size[v]

        for w in neighbors(v)
            if size[w] >= 1 && v != w
                delete!(set(size[w]), w)
                size[w] += 1
                pushfirst!(set(size[w]), w)
            end
        end

        j += 1

        while j >= 1 && isempty(set(j))
            j -= 1
        end
    end

    alpha, size
end


function Base.show(io::IO, ::MIME"text/plain", list::T) where T <: LinkedList
    items = pushfirst!(map(string, take(list, MAX_ITEMS_PRINTED)), "head")

    if MAX_ITEMS_PRINTED < length(items)
        items[end] = "..."
    end

    println(io, T)
    println(io, join(items, " ↔ "))
end


#######################
# Iteration Interface #
#######################


function Base.iterate(list::LinkedList, i::Integer=list.head[])
    iszero(i) ? nothing : (i, list.next[i])
end


function Base.IteratorSize(::Type{T}) where T <: LinkedList
    Base.SizeUnknown()
end


function Base.eltype(::Type{T}) where T <: LinkedList
    Int
end
