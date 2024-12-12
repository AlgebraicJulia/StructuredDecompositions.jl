struct DisjointSets
    sets::IntDisjointSets{Int}
    index::Vector{Int}
    root::Vector{Int}

    function DisjointSets(n::Integer)
        new(IntDisjointSets(n), collect(1:n), collect(1:n))
    end
end


function find!(sets::DisjointSets, u::Integer)
    sets.index[find_root!(sets.sets, u)]
end


function Base.union!(sets::DisjointSets, u::Integer, v::Integer)
    w = max(u, v)
    sets.root[w] = root_union!(sets.sets, sets.root[u], sets.root[v])
    sets.index[sets.root[w]] = w
end
