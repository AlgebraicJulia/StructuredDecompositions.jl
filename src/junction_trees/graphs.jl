# A simple directed graph.
struct Graph{V <: Integer, E <: Integer}
    ptr::Vector{E}
    tgt::Vector{V}
end


function Graph{V, E}(matrix::SparseMatrixCSC) where {V, E}
    # validate argument
    size(matrix, 1) != size(matrix, 2) && throw(ArgumentError("size(matrix, 1) != size(matrix, 2)"))

    # construct graph
    Graph{V, E}(getcolptr(matrix), rowvals(matrix))
end


function Graph(matrix::SparseMatrixCSC{<:Any, I}) where I
    Graph{I, I}(getcolptr(matrix), rowvals(matrix))
end


function SparseArrays.SparseMatrixCSC{T, I}(graph::Graph) where {T, I}
    nzval = Vector{T}(undef, ne(graph))
    SparseMatrixCSC{T, I}(nv(graph), nv(graph), graph.ptr, graph.tgt, nzval)
end


function SparseArrays.SparseMatrixCSC{T}(graph::Graph{V, E}) where {T, V, E}
    I = promote_type(V, E)
    SparseMatrixCSC{T, I}(graph)
end


function vtype(graph::Graph{V}) where V
    V
end


function etype(graph::Graph{<:Any, E}) where E
    E
end


function nv(graph::Graph{V}) where V
    convert(V, length(graph.ptr) - 1)
end


function ne(graph::Graph{<:Any, E}) where E
    convert(E, length(graph.tgt))
end


function vertices(graph::Graph{V}) where V
    OneTo{V}(nv(graph))
end


function neighbors(graph::Graph{<:Any, E}, i::Integer) where E
    range = graph.ptr[i]:graph.ptr[i + 1] - one(E)
    @view graph.tgt[range]
end


# Simple Linear-Time Algorithms to Test Chordality of Graphs, Test Acyclicity of Hypergraphs, and Selectively Reduce Acyclic Hypergraphs
# Tarjan and Yannakakis
# Maximum Cardinality Search
#
# Construct a fill-reducing permutation of a simple graph.
# The complexity is O(m + n), where m = |E| and n = |V|.
function mcs(graph::Graph{V}) where V
    # construct disjoint sets data structure
    head = zeros(V, nv(graph) + 1)
    prev = Vector{V}(undef, nv(graph) + 1)
    next = Vector{V}(undef, nv(graph) + 1)
   
    function set(i)
        LinkedList(view(head, i), prev, next)
    end
   
    # run algorithm
    alpha = Vector{V}(undef, nv(graph))
    size = Vector{V}(undef, nv(graph))

    for v in reverse(vertices(graph))
        size[v] = 1
        pushfirst!(set(1), v)
    end

    j = 1

    for i in reverse(vertices(graph))
        v = popfirst!(set(j))
        alpha[v] = i
        size[v] = 0

        for w in neighbors(graph, v)
            if size[w] >= 1
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

    alpha
end


# Algorithms for Sparse Linear Systems
# Scott and Tuma
# Algorithm 8.3: CM and RCM algorithms for band and profile reduction
#
# Compute the reverse Cuthill-Mckee ordering of a connected simple graph.
function rcm(graph::Graph)
    rcm!(copy(graph))
end


function rcm!(graph::Graph{V}) where V
    # sort neighbors
    degree = diff(graph.ptr)

    for j in vertices(graph)
        sort!(neighbors(graph, j); by=i -> degree[i])
    end

    # run algorithm
    root::V = argmin(degree)
    reverse!(bfs(graph, root))
end


# Perform a breadth-first search of a connected graph.
function bfs(graph::Graph{V}, root::V) where V
    label = fill(false, nv(graph))
    order = sizehint!(V[], nv(graph))

    label[root] = true
    push!(order, root)
    
    for j in order
        for i in neighbors(graph, j)
            if !label[i]
                label[i] = true
                push!(order, i)
            end
        end
    end
    
   order
end


# Direct Methods for Sparse Linear Systems §2.11
# Davis
# cs_symperm
#
# Permute an ordered graph.
function sympermute(graph, index::AbstractVector, order::Ordering)
    sympermute(Graph(graph), index, order)
end


function sympermute(graph::G, index::AbstractVector, order::Ordering) where G <: Graph
    sympermute!(zero(G), graph, index, order)
end


function sympermute!(target::Graph{V, E}, graph::Graph{V, E}, index::AbstractVector{V}, order::Ordering) where {V, E}
    # validate arguments
    vertices(graph) != eachindex(index) && throw(ArgumentError("vertices(graph) != eachindex(index)"))

    # compute column counts
    total::E = 0
    count = Vector{E}(undef, nv(graph) + 1)
    count[begin]  = 1
    count[2:end] .= 0

    @inbounds for j in vertices(graph)
        for i in neighbors(graph, j)
            if lt(order, i, j)
                u = index[i]
                v = index[j]

                if lt(order, v, u)
                    u, v = v, u
                end

                total += 1
                count[v + 1] += 1
            end
        end
    end

    # permute graph
    resize!(target.tgt, total)
    copy!(count, cumsum!(resize!(target.ptr, nv(graph) + 1), count))
    
    @inbounds for j in vertices(graph)
        for i in neighbors(graph, j)
            if lt(order, i, j)
                u = index[i]
                v = index[j]

                if lt(order, v, u)
                    u, v = v, u
                end

                target.tgt[count[v]] = u
                count[v] += 1
            end
        end
    end

    target
end


# Reverse the edges of a simple directed graph.
function Base.reverse(graph::G) where G <: Graph
    reverse!(zero(G), graph)
end


function Base.reverse!(target::Graph{V, E}, graph::Graph{V, E}) where {V, E}
    # compute column counts
    count = Vector{E}(undef, nv(graph) + 1)
    count[begin]  = 1
    count[2:end] .= 0

    @inbounds for j in vertices(graph)
        for i in neighbors(graph, j)
            count[i + 1] += 1
        end
    end

    # permute graph
    resize!(target.tgt, ne(graph))
    copy!(count, cumsum!(resize!(target.ptr, nv(graph) + 1), count))

    @inbounds for j in vertices(graph)
        for i in neighbors(graph, j)
            target.tgt[count[i]] = j
            count[i] += 1
        end
    end

    target
end


function Base.copy(graph::Graph)
    Graph(copy(graph.ptr), copy(graph.tgt))
end


function Base.zero(::Type{Graph{V, E}}) where {V, E}
    Graph(E[1], V[])
end
