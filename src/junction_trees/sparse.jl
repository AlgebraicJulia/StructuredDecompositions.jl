# Permute a symmetric matrix.
function symperm(matrix::SparseMatrixCSC, permutation::AbstractVector, order::Ordering=ForwardOrdering())
    invsymperm(matrix, invperm(permutation), order)
end


# Permute a symmetric matrix.
function invsymperm(matrix::SparseMatrixCSC, index::AbstractVector, order::Ordering=ForwardOrdering())
    invsymperm!(similar(matrix), matrix, index, order)
end


# Permute a symmetric matrix.
function symperm!(target::SparseMatrixCSC, source::SparseMatrixCSC, permutation::AbstractVector, order::Ordering=ForwardOrdering())
    invsymperm!(target, source, invperm(permutation), order)
end


# Direct Methods for Sparse Linear Systems §2.11
# Davis
# cs_symperm
function invsymperm!(target::SparseMatrixCSC, source::SparseMatrixCSC, index::AbstractVector, order::Ordering=ForwardOrdering())
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
    siz = Vector{Int}(undef, size(matrix, 2))
    set = Vector{LinkedList{Int}}(undef, size(matrix, 2))
    pointer = Vector{ListNode{Int}}(undef, size(matrix, 2))

    for i in axes(matrix, 2)
        siz[i] = 1
        set[i] = LinkedList{Int}()
        pointer[i] = push!(set[1], i)
    end

    i = size(matrix, 2)
    j = 1

    while i >= 1
        v = first(set[j])
        deleteat!(set[j], pointer[v])
        α[i] = v
        siz[v] = 0

        for w in view(rowvals(matrix), nzrange(matrix, v))
            if siz[w] >= 1
                deleteat!(set[siz[w]], pointer[w])
                siz[w] += 1
                pointer[w] = push!(set[siz[w]], w)
            end
        end

        i -= 1
        j += 1

        while j >= 1 && isempty(set[j])
            j -= 1
        end
    end

    α
end
