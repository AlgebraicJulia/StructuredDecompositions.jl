"""
    Order <: AbstractVector{Int}

A [permutation](https://en.wikipedia.org/wiki/Permutation) ``\\sigma`` of the set ``\\{1, \\dots, n\\}``.
This type implements the [abstract array interface](https://docs.julialang.org/en/v1/manual/interfaces/#man-interface-array).
"""
struct Order <: AbstractVector{Int}
    order::Vector{Int} # permutation
    index::Vector{Int} # inverse permutation
end


"""
    Order(order::AbstractVector)

Construct a permutation ``\\sigma`` from a sequence ``(\\sigma(1), \\dots, \\sigma(n))``.
"""
function Order(order::AbstractVector)
    n = length(order)
    index = Vector{Int}(undef, n)

    for i in 1:n
        index[order[i]] = i
    end
    
    Order(order, index)
end


# Construct a copy of a permutation.
function Order(order::Order)
    Order(order.order, order.index)
end


# Compose two permutations.
function compose(left::Order, right::Order)
    Order(right.order[left.order], left.index[right.index])
end


"""
    permutation(order::Order, i::Integer)

Get the element ``\\sigma(i)``.
"""
function permutation(order::Order, i::Integer)
    order[i]
end


"""
    permutation(order::Order, i::AbstractVector)

Get the elements ``(\\sigma(i_1), \\dots, \\sigma(i_n))``.
"""
function permutation(order::Order, i::AbstractVector)
    view(order, i)
end


"""
    inverse(order::Order, v::Integer)

Get the index ``\\sigma^{-1}(v)``.
"""
function inverse(order::Order, v::Integer)
    order.index[v]
end


"""
    inverse(order::Order, v::AbstractVector)

Get the indices ``(\\sigma^{-1}(v_1), \\dots, \\sigma^{-1}(v_n))``.
"""
function inverse(order::Order, v::AbstractVector)
    view(order.index, v)
end


"""
    inverse(order::Order)

Construct the inverse permutation ``\\sigma^{-1}``.
"""
function inverse(order::Order)
    Order(order.index, order.order)
end


#############################
# Abstract Vector Interface #
#############################


function Base.getindex(order::Order, i)
    order.order[i]
end


function Base.IndexStyle(::Type{Order})
    IndexLinear()
end


function Base.size(order::Order)
    (length(order.order),)
end
