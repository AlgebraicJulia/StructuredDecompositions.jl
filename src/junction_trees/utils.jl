# See sympermute!.
function sympermute(matrix::SparseMatrixCSC, index::AbstractVector, order::Ordering=ForwardOrdering())
    sympermute!(similar(matrix), matrix, index, order)
end


# Direct Methods for Sparse Linear Systems §2.11
# Davis
# cs_symperm
#
# Permute the rows and columns of a symmetric matrix, represented by its lower / upper triangular part `source`.
# The permutation is represented by its inverse `index`.
# Use ForwardOrdering() if `source` is upper triangular and ReverseOrdering() if it is lower triangular.
# The result is stored in `target`.
function sympermute!(target::SparseMatrixCSC, source::SparseMatrixCSC, index::AbstractVector, order::Ordering=ForwardOrdering())
    count = similar(getcolptr(source))
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

    resize!(rowvals(target), nnz(source))
    resize!(nonzeros(target), nnz(source))
    copy!(count, cumsum!(getcolptr(target), count))

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


# Compute the union of sorted sets `source1` and `source2`.
# The union is stored in `target`.
# `source1` and `source2` should support iteration, and `target` should support `push!(target, v)`.
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
