# A rooted tree.
# This type implements the abstract graph and indexed tree interfaces.
struct Tree <: AbstractTree
    parent::Vector{Int}  # vector of parents
    child::Vector{Int}   # vector of left-children
    brother::Vector{Int} # vector of right-siblings
    root::Array{Int, 0}  # root
end


function Tree()
    Tree(Int[])
end


# Construct a rooted tree from a list of parents.
function Tree(parent::AbstractVector)
    child = zeros(Int, length(parent))
    brother = Vector{Int}(undef, length(parent))
    root = zeros(Int)

    for (i, j) in Iterators.reverse(enumerate(parent))
        if iszero(j)
            brother[i] = 0
            root[] = i
        else
            brother[i] = child[j]
            child[j] = i
        end
    end

    Tree(parent, child, brother, root)
end


# Permute the vertices of a rooted tree.
function Base.permute!(tree::Tree, permutation::AbstractVector)
    index = invperm(permutation)
    fill!(tree.child, 0)

    tree.parent .= map(permutation) do i
        j = tree.parent[i]
        iszero(j) ? 0 : index[j]
    end

    for (i, j) in Iterators.reverse(enumerate(tree.parent))
        if iszero(j)
            tree.brother[i] = 0
            tree.root[] = i
        else
            tree.brother[i] = tree.child[j]
            tree.child[j] = i
        end
    end

    tree
end


# Compact Clique AbstractTree Data Structures in Sparse Matrix Factorizations
# Pothen and Sun
# Figure 4: The Clique AbstractTree Algorithm 2
function pothensun(etree::AbstractTree, colcount::AbstractVector, stype::SupernodeType)
    n = nv(etree)
    new_in_clique = Vector{Int}(undef, n)
    new = sizehint!(Int[], n)
    parent = sizehint!(Int[], n)
    ancestor = sizehint!(Int[], n)

    for v in 1:n
        u = child_in_supernode(etree, colcount, stype, v)

        if !isnothing(u)
            new_in_clique[v] = new_in_clique[u]
        else
            new_in_clique[v] = length(new) + 1
            push!(new, v)
            push!(parent, 0)
            push!(ancestor, 0)
        end

        for s in childindices(etree, v)
            if s !== u
                parent[new_in_clique[s]] = new_in_clique[v]
                ancestor[new_in_clique[s]] = v
            end
        end
    end

    new, ancestor, Tree(parent)
end


function Base.:(==)(left::Tree, right::Tree)
    left.parent == right.parent
end


##########################
# Indexed Tree Interface #
##########################


function firstchildindex(tree::Tree, i::Integer)
    ztn(tree.child[i])
end


function AbstractTrees.rootindex(tree::Tree)
    ztn(tree.root[])
end


function AbstractTrees.parentindex(tree::Tree, i::Integer)
    ztn(tree.parent[i])
end


function AbstractTrees.nextsiblingindex(tree::Tree, i::Integer)
    ztn(tree.brother[i])
end


function AbstractTrees.nodevalue(tree::Tree, i::Integer)
    i
end


############################
# Abstract Graph Interface #
############################


function Graphs.nv(tree::Tree)
    length(tree.parent)
end
