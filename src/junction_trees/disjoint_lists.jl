# A collection of disjoint linked lists.
struct DisjointLists
    last::Vector{Int}
    prev::Vector{Int}
    next::Vector{Int}
end


# Construct n empty lists.
function DisjointLists(n::Integer)
    last = zeros(Int, n)
    prev = Vector{Int}(undef, n)
    next = Vector{Int}(undef, n)
    DisjointLists(last, prev, next)
end


# Is the ith list empty?
function Base.isempty(lists::DisjointLists, i::Integer)
    iszero(lists.last[i])
end


# Append an element v to the end of the ith list.
function Base.push!(lists::DisjointLists, i::Integer, v::Integer)
    p = lists.last[i]
    lists.last[i] = v
    lists.prev[v] = p
    lists.next[v] = 0
    
    if !iszero(p)
        lists.next[p] = v
    end
    
    lists
end


# Remove the last element of the ith list.
function Base.pop!(lists::DisjointLists, i::Integer)
    v = lists.last[i]
    p = lists.prev[v]
    lists.last[i] = p
    
    if !iszero(p)
        lists.next[p] = 0
    end
    
    v
end


# Delete an element v from the ith list.
function Base.delete!(lists::DisjointLists, i::Integer, v::Integer)
    p = lists.prev[v]
    n = lists.next[v]
    
    if !iszero(p)
        lists.next[p] = n
    end
    
    if !iszero(n)
        lists.prev[n] = p
    else
        lists.last[i] = p
    end
    
    lists
end


# Simple Linear-Time Algorithms to Test Chordality of Graphs, Test Acyclicity of Hypergraphs, and Selectively Reduce Acyclic Hypergraphs
# Tarjan and Yannakakis
# Maximum Cardinality Search
#
# Construct a fill-reducing permutation of a simple graph, represented by its adjacency matrix.
function mcs(matrix::SparseMatrixCSC)
    order = Vector{Int}(undef, size(matrix, 2))
    number = Vector{Int}(undef, size(matrix, 2))
    set = DisjointLists(size(matrix, 2) + 1)

    for v in axes(matrix, 2)
        number[v] = 1
        push!(set, 1, v)
    end

    i = size(matrix, 2)
    j = 1

    while i >= 1
        v = pop!(set, j)
        order[i] = v
        number[v] = 0

        for w in @view rowvals(matrix)[nzrange(matrix, v)]
            if number[w] >= 1
                delete!(set, number[w], w)
                number[w] += 1
                push!(set, number[w], w)
            end
        end

        i -= 1
        j += 1

        while j >= 1 && isempty(set, j)
            j -= 1
        end
    end

    order
end
