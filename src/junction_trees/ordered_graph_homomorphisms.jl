struct OrderedGraphHomomorphism{Source <: AbstractOrderedGraph, Target <: AbstractOrderedGraph}
    source::Source
    target::Target
    mapping::Vector{Int}
    indptr::Vector{Int}
end
