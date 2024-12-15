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
    index = Vector{Int}(undef, length(order))
    index[order] = eachindex(order)
    Order(order, index)
end


function Order(::UndefInitializer, n::Integer)
    Order(Vector{Int}(undef, n), Vector{Int}(undef, n))
end


function Base.permute!(order::Order, permutation::AbstractVector)
    permute!(order.order, permutation)
    order.index[order.order] = eachindex(order.order)
    order
end


# Construct a copy of a permutation.
function Base.copy(order::Order)
    Order(order.order, order.index)
end


function Base.deepcopy(order::Order)
    Order(copy(order.order), copy(order.index))
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


function Base.setindex!(order::Order, v, i)
    order.index[v] = i
    order.order[i] = v
end


function Base.IndexStyle(::Type{Order})
    IndexLinear()
end


function Base.size(order::Order)
    (length(order.order),)
end
