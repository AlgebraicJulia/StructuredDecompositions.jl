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


function flowcutter(matrix::SparseMatrixCSC, time::Integer, seed::Integer)
    mktempdir(dirname(FlowCutterPACE17_jll.flow_cutter_pace17_path)) do directory
        input = directory * "/input.gr"
        output = directory * "/output.td"
        open(io -> writegr(io, matrix), input; write=true)

        command = `$(FlowCutterPACE17_jll.flow_cutter_pace17()) -s $seed`
        process = run(pipeline(input, command, output); wait=false)

        while !process_running(process)
            sleep(1)
        end

        sleep(time)
        kill(process)
        open(readtd, output)
    end
end


# Write a simple graph to a .gr file.
# https://pacechallenge.org/2017/treewidth/
function writegr(io::IO, matrix::SparseMatrixCSC)
    m = size(matrix, 2)
    n = nnz(matrix) ÷ 2
    write(io, "p tw $m $n\n")

    for j in axes(matrix, 2)
        for i in @view rowvals(matrix)[nzrange(matrix, j)]
            if j < i
                write(io, "$j $i\n")
            end
        end
    end
end


# Read a tree decomposition from a .td file.
# https://pacechallenge.org/2017/treewidth/
function readtd(io::IO) 
    # comments
    comment = String[]
    line = readline(io)
    
    while !isempty(line)
        flag = first(line)
        
        if flag == 'c'
            push!(comment, line[3:end])
            line = readline(io)
            continue
        end
        
        break
    end

    # statistics
    if isempty(line) || first(line) != 's'
        throw(EOFError("Missing header."))
    end

    nb, tw, nv = imap(word -> parse(Int, word), drop(split(line), 2))
    line = readline(io)
    
    # bags
    bagptr = sizehint!(Int[], nb + 1)    
    bagval = sizehint!(Int[], tw * nb ÷ 2)
    push!(bagptr, 1)
    
    while !isempty(line)
        flag = first(line)
        
        if flag == 'c'
            push!(comment, line[3:end])
            line = readline(io)
            continue
        end
        
        if flag == 'b' 
            for v in imap(word -> parse(Int, word), drop(split(line), 2))         
                push!(bagval, v)
            end
            
            push!(bagptr, length(bagval) + 1)
            line = readline(io)
            continue
        end
        
        break
    end

    # edges
    tree = map(i -> Int[], 1:nb)

    while !isempty(line)
        flag = first(line)
   
        if flag == 'c'
            push!(comment, line[3:end])
            line = readline(io)
            continue
        end

        i, j = imap(word -> parse(Int, word), split(line))
        push!(tree[i], j)
        push!(tree[j], i)
        line = readline(io)
    end
   
    comment, nb, tw, nv, bagptr, bagval, tree
end
