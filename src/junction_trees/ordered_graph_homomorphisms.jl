struct OrderedGraphHomomorphism{Source <: AbstractOrderedGraph, Target <: AbstractOrderedGraph}
    source::Source
    target::Target
    mapping::Vector{Int}
    index::Vector{Int}
end


struct JunctionTree
    order::Order
    homomorphism::OrderedGraphHomomorphism{OrderedGraph, FilledGraph}
end
