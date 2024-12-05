struct OrderedGraphHomomorphism{Domain, Codomain}
    domain::Domain
    codomain::Codomain
    map::Vector{Int}
    index::Vector{Int}
end
