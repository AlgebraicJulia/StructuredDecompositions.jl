# Permute a symmetric matrix.
function sympermute(matrix::SparseMatrixCSC, index::AbstractVector, order::Ordering=ForwardOrdering())
    sympermute!(similar(matrix), matrix, index, order)
end


# Direct Methods for Sparse Linear Systems §2.11
# Davis
# cs_symperm
function sympermute!(target::SparseMatrixCSC, source::SparseMatrixCSC, index::AbstractVector, order::Ordering=ForwardOrdering())
    count = similar(source.colptr)
    count[1] = 1
    count[2:end] .= 0

    for j in axes(source, 2)
        for p in nzrange(source, j)
            i = rowvals(source)[p]

            if lt(order, i, j) || isequal(i, j)
                if lt(order, index[i], index[j])
                    v = index[j]
                else
                    v = index[i]
                end

                count[v + 1] += 1
            end
        end
    end

    copy!(count, cumsum!(target.colptr, count))

    for j in axes(source, 2)
        for p in nzrange(source, j)
            i = rowvals(source)[p]
            x = nonzeros(source)[p]

            if lt(order, i, j) || isequal(i, j)
                if lt(order, index[i], index[j])
                    u = index[i]
                    v = index[j]
                else
                    u = index[j]
                    v = index[i]
                end

                rowvals(target)[count[v]] = u
                nonzeros(target)[count[v]] = x
                count[v] += 1
            end
        end
    end

    target
end


# Simple Linear-Time Algorithms to Test Chordality of Graphs, Test Acyclicity of Hypergraphs, and Selectively Reduce Acyclic Hypergraphs
# Tarjan and Yannakakis
# Maximum Cardinality Search
function mcs(matrix::SparseMatrixCSC)
    α = Vector{Int}(undef, size(matrix, 2))
    num = Vector{Int}(undef, size(matrix, 2))
    set = DisjointLists(size(matrix, 2))

    for v in axes(matrix, 2)
        num[v] = 1
        push!(set, 1, v)
    end

    i = size(matrix, 2)
    j = 1

    while i >= 1
        v = pop!(set, j)
        α[i] = v
        num[v] = 0

        for w in view(rowvals(matrix), nzrange(matrix, v))
            if num[w] >= 1
                delete!(set, num[w], w)
                num[w] += 1
                push!(set, num[w], w)
            end
        end

        i -= 1
        j += 1

        while j >= 1 && isempty(set, j)
            j -= 1
        end
    end

    α
end
