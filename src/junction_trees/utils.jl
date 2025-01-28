# Direct Methods for Sparse Linear Systems §2.11
# Davis
# cs_symperm
#
# Permute an ordered graph.
function sympermute(graph, index::AbstractVector{I}, order::Ordering) where I
    target = spzeros(Nothing, I, length(index), length(index))
    sympermute!(target, graph, index, order)
end


function sympermute!(target::SparseMatrixCSC, matrix::SparseMatrixCSC, index::AbstractVector, order::Ordering)
    # validate arguments
    size(target) != size(matrix) && throw(ArgumentError("size(target) != size(matrix)"))

    # run algorithm
    sympermute!(target, index, order) do j
        @view rowvals(matrix)[nzrange(matrix, j)]
    end
end


function sympermute!(neighbors::Function, target::SparseMatrixCSC{Nothing, I}, index::AbstractVector{I}, order::Ordering) where I
    # validate arguments
    eachindex(index) != axes(target, 1) && throw(ArgumentError("eachindex(index) != axes(target, 1)"))
    eachindex(index) != axes(target, 2) && throw(ArgumentError("eachindex(index) != axes(target, 2)"))

    # compute column counts
    total::I = 0
    count = Vector{I}(undef, size(target, 2) + 1)
    count[1] = 1
    count[2:end] .= 0

    for j in axes(target, 2)
        for i in neighbors(j)
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

    # resize target
    resize!(rowvals(target), total)
    resize!(nonzeros(target), total)
    copy!(count, cumsum!(getcolptr(target), count))

    # permute graph
    for j in axes(target, 2)
        for i in neighbors(j)
            if lt(order, i, j)
                u = index[i]
                v = index[j]
                
                if lt(order, v, u)
                    u, v = v, u
                end

                rowvals(target)[count[v]] = u
                count[v] += 1
            end
        end
    end

    target
end


# Algorithms for Sparse Linear Systems
# Scott and Tuma
# Algorithm 8.3: CM and RCM algorithms for band and profile reduction
#
# Compute the reverse Cuthill-Mckee ordering of a connected simple graph.
function rcm(matrix::SparseMatrixCSC)
    rcm!(copy(matrix))
end


# Compute the reverse Cuthill-Mckee ordering of a connected simple graph.
function rcm!(matrix::SparseMatrixCSC{<:Any, I}) where I
    # validate argument
    size(matrix, 1) != size(matrix, 2) && throw(ArgumentError("size(matrix, 1) != size(matrix, 2)"))
    
    # sort neighbors
    degree = diff(getcolptr(matrix))
    
    function neighbors(j)
        @view rowvals(matrix)[nzrange(matrix, j)]
    end
    
    for j in eachindex(degree)
        sort!(neighbors(j); by=i -> degree[i])
    end
    
    # run algorithm
    vertices::OneTo{I} = axes(matrix, 2)
    root::I = argmin(degree)
    reverse!(bfs(neighbors, vertices, root))
end


# Permform a breadth-first search of a connected simple graph.
function bfs(neighbors::Function, vertices::AbstractVector{I}, root::I) where I
    label = fill(false, length(vertices))
    order = sizehint!(I[], length(vertices))

    label[root] = true
    push!(order, root)
    
    for v in order
        for w in neighbors(v)
            if !label[w]
                label[w] = true
                push!(order, w)
            end
        end
    end
    
   order
end


# Compute the union of sorted sets `source1` and `source2`.
# The result is appended to `target`.
function mergesorted!(target, source1, source2, order::Ordering=ForwardOrdering())
    i1 = iterate(source1)
    i2 = iterate(source2)
   
    while !isnothing(i1) && !isnothing(i2)
        x1, s1 = i1
        x2, s2 = i2

        if isequal(x1, x2)
            push!(target, x1)
            i1 = iterate(source1, s1)
            i2 = iterate(source2, s2)
        elseif lt(order, x1, x2)
            push!(target, x1)
            i1 = iterate(source1, s1)
        else
            push!(target, x2)
            i2 = iterate(source2, s2)
        end
       
    end
   
    while !isnothing(i1)
        x1, s1 = i1
        push!(target, x1)
        i1 = iterate(source1, s1)
    end

    while !isnothing(i2)
        x2, s2 = i2
        push!(target, x2)
        i2 = iterate(source2, s2)
    end
   
    target
end
