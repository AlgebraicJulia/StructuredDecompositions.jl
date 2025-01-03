# A stack of fixed size.
# This type implements the abstract vector interface.
struct FixedStack{T} <: AbstractVector{T}
    top::Array{Int, 0}
    items::Vector{T}

    function FixedStack{T}(n::Integer) where T
        new(zeros(Int), Vector{T}(undef, n))
    end
end


function Base.push!(stack::FixedStack, v)
    stack.top[] += 1
    stack.items[stack.top[]] = v
    v
end


function Base.pop!(stack::FixedStack)
    v = stack.items[stack.top[]]
    stack.top[] -= 1
    v
end


#############################
# Abstract Vector Interface #
#############################


function Base.getindex(stack::FixedStack, i)
    stack.items[i]
end


function Base.IndexStyle(::Type{FixedStack})
    IndexLinear()
end


function Base.size(stack::FixedStack)
    (stack.top[],)
end
