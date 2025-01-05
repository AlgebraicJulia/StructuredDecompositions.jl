# A collection of disjoint sets.
# Each set is represented by its largest element.
struct DisjointTrees
    sets::IntDisjointSets{Int}
    root::Vector{Int}
    index::Vector{Int}
end


# Construct the sets {1}, {2}, ..., {n}.
function DisjointTrees(n::Integer)
    DisjointTrees(IntDisjointSets(n), collect(1:n), collect(1:n))
end


# Get the index of the set represented by the element u.
function find!(sets::DisjointTrees, u::Integer)
    v = find_root!(sets.sets, u)
    sets.index[v]
end


# Merge the sets represented by the elements u < v.
function Base.union!(sets::DisjointTrees, u::Integer, v::Integer)
    sets.root[v] = root_union!(sets.sets, sets.root[u], sets.root[v])
    sets.index[sets.root[v]] = v
end
