# See sympermute!.
function sympermute(matrix::SparseMatrixCSC, index::AbstractVector, order::Ordering=ForwardOrdering())
    sympermute!(similar(matrix), matrix, index, order)
end


# Direct Methods for Sparse Linear Systems §2.11
# Davis
# cs_symperm
#
# Permute the rows and columns of a symmetric matrix.
# - `target`: overwritten by the permuted matrix
# - `source`: upper / lower triangular part of the matrix to be permuted
# - `index`: inverse permutation
# - `order`: `ForwardOrdering()` if `matrix` is upper triangular, `ReverseOrdering()` if it is lower triangular
function sympermute!(target::SparseMatrixCSC, source::SparseMatrixCSC, index::AbstractVector, order::Ordering=ForwardOrdering())
    # validate arguments
    eachindex(index) != axes(target, 1) && throw(ArgumentError("eachindex(index) != axes(target, 1)"))
    eachindex(index) != axes(target, 2) && throw(ArgumentError("eachindex(index) != axes(target, 2)"))
    size(target) != size(source) && throw(ArgumentError("size(target) != size(source)"))

    # resize target array
    resize!(rowvals(target), nnz(source))
    resize!(nonzeros(target), nnz(source))
    
    # run algorithm
    count = similar(getcolptr(source))
    count[1] = 1
    count[2:end] .= 0

    for j in axes(source, 2)
        for p in nzrange(source, j)
            i = rowvals(source)[p]

            if lt(order, i, j) || isequal(i, j)
                u = index[i]
                v = index[j]
                
                if lt(order, v, u)
                    u, v = v, u
                end
                
                count[v + 1] += 1
            end
        end
    end

    copy!(count, cumsum!(getcolptr(target), count))

    for j in axes(source, 2)
        for p in nzrange(source, j)
            i = rowvals(source)[p]
            x = nonzeros(source)[p]

            if lt(order, i, j) || isequal(i, j)
                u = index[i]
                v = index[j]
                
                if lt(order, v, u)
                    u, v = v, u
                end

                q = count[v]
                rowvals(target)[q] = u
                nonzeros(target)[q] = x
                count[v] += 1
            end
        end
    end

    target
end


# A Spectral Algorithm for Envelope Reduction of Sparse Matrices
# Barnard, Pothen, and Simon
# Algorithm 1: Spectral Algorithm
#
# Compute the spectral ordering of a graph.
function spectralorder(matrix::SparseMatrixCSC; tol=0.0)
    matrix = SparseMatrixCSC{Float64, Int}(matrix)
    fill!(nonzeros(matrix), 1)
    fkeep!((i, j, v) -> i != j, matrix)
    value, vector = Laplacians.fiedler(matrix; tol)
    sortperm(reshape(vector, size(matrix, 2)))
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


# Compute the difference of sorted sets `source1` and `source2`.
# The result is appended to `target`.
function diffsorted!(target, source1, source2, order::Ordering=ForwardOrdering())
    i1 = iterate(source1)
    i2 = iterate(source2)

    while !isnothing(i1) && !isnothing(i2)
        x1, s1 = i1
        x2, s2 = i2

        if isequal(x1, x2)
            i1 = iterate(source1, s1)
            i2 = iterate(source2, s2)
        elseif lt(order, x1, x2)
            push!(target, x1)
            i1 = iterate(source1, s1)
        else
            i2 = iterate(source2, s2)
        end

    end

    while !isnothing(i1)
        x1, s1 = i1
        push!(target, x1)
        i1 = iterate(source1, s1)
    end

    target
end
