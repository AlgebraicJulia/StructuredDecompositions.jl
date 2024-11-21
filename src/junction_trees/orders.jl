"""
    Order <: AbstractVector{Int}

A permutation of the set ``\\{1, \\dots, n\\}.``
"""
struct Order <: AbstractVector{Int}
    order::Vector{Int} # permutation
    index::Vector{Int} # inverse permutation
end


"""
    Order(order::AbstractVector)

Construct a permutation ``\\sigma`` from a sequence ``(\\sigma(1), \\dots, \\sigma(n)).``
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


# Determine if i < j, where
#    u = σ(i)
#    v = σ(j)
function (order::Order)(u, v)
    inverse(order, u) < inverse(order, v)
end


# Compose two permutations.
function compose(left::Order, right::Order)
    Order(right.order[left.order], left.index[right.index])
end


# Construct the inverse permutation.
function inverse(order::Order)
    Order(order.index, order.order)
end


# Get the index σ⁻¹(v),
function inverse(order::Order, v)
    order.index[v]
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
