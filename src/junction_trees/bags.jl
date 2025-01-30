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


function indexinsorted!(target::AbstractVector{I}, source1::Bag{I}, source2::AbstractVector{I}) where I
    res = residual(source1)
    sep = separator(source1)
    s1 = s2 = 1

    while s2 in eachindex(source2)
        x2 = source2[s2]

        if x2 in res
            target[s2] = x2 - res[begin] + 1
            s2 += 1
        else
            break
        end
    end

    while s2 in eachindex(source2)
        x1 = sep[s1]
        x2 = source2[s2]

        if x1 >= x2
            target[s2] = s1 + length(res)
            s2 += 1
        end

        s1 += 1
    end

    target
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



