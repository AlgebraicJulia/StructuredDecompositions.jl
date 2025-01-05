# A stack of fixed size.
# This type implements the abstract vector interface.
struct Stack{T} <: AbstractVector{T}
    top::Array{Int, 0}
    items::Vector{T}

    function Stack{T}(n::Integer) where T
        new(zeros(Int), Vector{T}(undef, n))
    end
end


function Base.push!(stack::Stack, v)
    stack.top[] += 1
    stack.items[stack.top[]] = v
    v
end


function Base.pop!(stack::Stack)
    v = stack.items[stack.top[]]
    stack.top[] -= 1
    v
end


#############################
# Abstract Vector Interface #
#############################


function Base.getindex(stack::Stack, i::Integer)
    stack.items[i]
end


function Base.setindex!(stack::Stack, v, i::Integer)
    stack.items[i] = v
end


function Base.IndexStyle(::Type{Stack})
    IndexLinear()
end


function Base.size(stack::Stack)
    (stack.top[],)
end
