"""
    Bag{I} <: AbstractVector{I}

A bag of a junction tree.
"""
struct Bag{I <: Integer} <: AbstractVector{I}
    residual::UnitRange{I}
    separator::SubArray{I, 1, Vector{I}, Tuple{UnitRange{I}}, true}
end


"""
    residual(bag::Bag)

Get the residual of a bag.
"""
function residual(bag::Bag)
    bag.residual
end


"""
    separator(bag::Bag)

Get the separator of a bag.
"""
function separator(bag::Bag)
    bag.separator
end


#############################
# Abstract Vector Interface #
#############################


function Base.getindex(bag::Bag, i::Integer)
    r = residual(bag)
    s = separator(bag)
    i in eachindex(r) ? r[i] : s[i - length(r)]
end


function Base.IndexStyle(::Type{Bag})
    IndexLinear()
end


function Base.size(bag::Bag)
    (length(residual(bag)) + length(separator(bag)),)
end


function Base.in(v, bag::Bag)
    v in residual(bag) || insorted(v, separator(bag))
end


function Base.hasfastin(::Type{Bag})
    true
end



