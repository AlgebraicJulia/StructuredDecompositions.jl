# The direct sum of one or more vectors.
# This type implements the abstract array interface.
struct SumVector{N, T, Vectors <: NTuple{N, AbstractVector{T}}} <: AbstractVector{T}
    vectors::Vectors
end


function SumVector(vectors::AbstractVector{T}...) where T
    SumVector(vectors)
end


#############################
# Abstract Vector Interface #
#############################


function Base.getindex(sv::SumVector, i::Integer)
    for vector in sv.vectors
        n = length(vector)

        if i <= n
            return vector[i]
        else
            i -= n
        end
    end

    error()
end


function Base.IndexStyle(::Type{SumVector})
    IndexLinear()
end


function Base.size(sv::SumVector)
    (sum(length, sv.vectors),)
end
