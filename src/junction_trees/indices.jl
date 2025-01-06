# A stack of unique integers.
struct Index
    top::Array{Int, 0}
    index::Vector{Int}

    function Index(n::Integer)
        new(zeros(Int), Vector{Int}(undef, n))
    end
end


function Base.push!(stack::Index, v::Integer)
    stack.top[] += 1
    stack.index[v] = stack.top[]
    stack
end
