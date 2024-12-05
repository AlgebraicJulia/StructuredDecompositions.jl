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
    index[order] = 1:n
    Order(order, index)
end


# Compose two permutations.
function compose(left::Order, right::Order)
    Order(right.order[left.order], left.index[right.index])
end


# Construct a copy of a permutation.
function Base.copy(order::Order)
    Order(order.order, order.index)
end


"""
    inverse(order::Order)

Construct the inverse permutation ``\\sigma^{-1}``.
"""
function Base.inv(order::Order)
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
