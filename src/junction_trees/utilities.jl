# Direct Methods for Sparse Linear Systems §2.11
# Davis
# cs_symperm
function symperm(matrix::SparseMatrixCSC, index::AbstractVector)
    count = similar(matrix.colptr)
    rowval = similar(matrix.rowval)
    nzval = similar(matrix.nzval)

    count[1] = 1
    count[2:end] .= 0

    for j in axes(matrix, 2)
        for p in nzrange(matrix, j)
            i = rowvals(matrix)[p]

            if i >= j
                v = min(index[i], index[j])
                count[v + 1] += 1
            end
        end
    end

    colptr = cumsum(count)
    copy!(count, colptr)

    for j in axes(matrix, 2)
        for p in nzrange(matrix, j)
            i = rowvals(matrix)[p]

            if i >= j
                v = min(index[i], index[j])
                rowval[count[v]] = max(index[i], index[j])
                nzval[count[v]] = nonzeros(matrix)[p]
                count[v] += 1
            end
        end
    end

    SparseMatrixCSC(size(matrix)..., colptr, rowval, nzval)
end


# ntz(i) =  { i if i > 0
#           { 0 if i = nothing
function ntz(i::Int)
    i
end


# ntz(i) =  { i if i > 0
#           { 0 if i = nothing
function ntz(i::Nothing)
    0
end


# ntz(i) =  { i       if i > 0
#           { nothing if i = 0
function ztn(i::Int)
    if !iszero(i)
        i
    end
end
