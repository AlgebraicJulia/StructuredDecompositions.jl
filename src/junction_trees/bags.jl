struct Bag <: AbstractVector{Int}
    residual::UnitRange{Int}
    separator::SubArray{Int, 1, Vector{Int}, Tuple{UnitRange{Int}}, true}
end


function residual(bag::Bag)
    bag.residual
end


function separator(bag::Bag)
    bag.separator
end


#############################
# Abstract Vector Interface #
#############################


function Base.getindex(bag::Bag, i::Integer)
    res = residual(bag)
    sep = separator(bag)
    i in eachindex(res) ? res[i] : sep[i - length(res)]
end


function Base.IndexStyle(::Type{Bag})
    IndexLinear()
end


function Base.size(bag::Bag)
    (length(residual(bag)) + length(separator(bag)),)
end

