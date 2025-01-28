"""
    SupernodeType

A type of supernode partition. The options are

| type                  | name                            |
| :-------------------- | :------------------------------ |
| [`Nodal`](@ref)       | nodal supernode partition       |
| [`Maximal`](@ref)     | maximal supernode partition     |
| [`Fundamental`](@ref) | fundamental supernode partition |
"""
abstract type SupernodeType end


"""
    Nodal <: SupernodeType

A nodal  supernode partition.
"""
struct Nodal <: SupernodeType end


"""
    Maximal <: SupernodeType

A maximal supernode partition.
"""
struct Maximal <: SupernodeType end


"""
    Fundamental <: SupernodeType

A fundamental supernode partition.
"""
struct Fundamental <: SupernodeType end


# Compact Clique Tree Data Structures in Sparse Matrix Factorizations
# Pothen and Sun
# Figure 4: The Clique Tree Algorithm 2
#
# Compute the maximal supernode partition of the montone transitive extension of an ordered graph.
# The complexity is O(n), where n = |V|.
function stree(tree::Tree{I}, colcount::AbstractVector{I}, snd::Maximal) where I
    # validate arguments
    tree != eachindex(colcount) && throw(ArgumentError("tree != eachindex(colcount)"))

    # run algorithm
    new = sizehint!(I[], length(tree))
    parent = sizehint!(I[], length(tree))
    ancestor = sizehint!(I[], length(tree))
    new_in_clique = Vector{I}(undef, length(tree))

    for v in tree
        u = nothing

        for s in childindices(tree, v)
            if colcount[s] == colcount[v] + 1
                u = s
                break
            end
        end

        if !isnothing(u)
            new_in_clique[v] = new_in_clique[u]

            for s in childindices(tree, v)
                if s !== u
                    parent[new_in_clique[s]] = new_in_clique[v]
                    ancestor[new_in_clique[s]] = v
                end
            end
        else
            push!(new, v)
            push!(parent, 0)
            push!(ancestor, 0)
            new_in_clique[v] = length(new)

            for s in childindices(tree, v)
                parent[new_in_clique[s]] = new_in_clique[v]
                ancestor[new_in_clique[s]] = v
            end
        end
    end

    new, ancestor, Tree(parent)
end


# Compute the fundamental supernode partition of the montone transitive extension of an ordered graph.
# The complexity is O(n), where n = |V|.
function stree(tree::Tree{I}, colcount::AbstractVector{I}, snd::Fundamental) where I
    # validate arguments
    tree != eachindex(colcount) && throw(ArgumentError("tree != eachindex(colcount)"))

    # run algorithm
    new = sizehint!(I[], length(tree))
    parent = sizehint!(I[], length(tree))
    ancestor = sizehint!(I[], length(tree))
    new_in_clique = Vector{I}(undef, length(tree))

    for v in tree
        u = firstchildindex(tree, v)

        if !isnothing(u) && colcount[u] == colcount[v] + 1 && isnothing(nextsiblingindex(tree, u))
            new_in_clique[v] = new_in_clique[u]
        else
            push!(new, v)
            push!(parent, 0)
            push!(ancestor, 0)
            new_in_clique[v] = length(new)

            for s in childindices(tree, v)
                parent[new_in_clique[s]] = new_in_clique[v]
                ancestor[new_in_clique[s]] = v
            end
        end
    end

    new, ancestor, Tree(parent)
end


# Compute the nodal supernode partition of the montone transitive extension of an ordered graph.
# The complexity is O(1).
function stree(tree::Tree, colcount::AbstractVector, snd::Nodal)
    # validate arguments
    tree != eachindex(colcount) && throw(ArgumentError("tree != eachindex(colcount)"))

    # run algorithm
    tree, tree.parent, tree
end


""" 
    DEFAULT_SUPERNODE_TYPE = Maximal()

The default supernode partition.
"""
const DEFAULT_SUPERNODE_TYPE = Maximal()
