"""
    Bag{V, E} <: AbstractVector{V}

A bag of a junction tree.
"""
struct Bag{V <: Integer, E <: Integer} <: AbstractVector{V}
    residual::UnitRange{V}
    separator::SubArray{V, 1, Vector{V}, Tuple{UnitRange{E}}, true}
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
    res = residual(bag)
    sep = separator(bag)

    if i in eachindex(res)
        res[i]
    else
        sep[i - length(res)]
    end
end


function Base.IndexStyle(::Type{<:Bag})
    IndexLinear()
end


function Base.size(bag::Bag)
    (length(residual(bag)) + length(separator(bag)),)
end


function Base.in(v, bag::Bag)
    v in residual(bag) || insorted(v, separator(bag))
end


function Base.hasfastin(::Type{<:Bag})
    true
end



