# A rooted tree.
# This type implements the abstract graph and indexed tree interfaces.
struct Tree <: AbstractUnitRange{Int}
    parent::Vector{Int}  # vector of parents
    child::Vector{Int}   # vector of left-children
    brother::Vector{Int} # vector of right-siblings
    root::Array{Int, 0}  # root
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


function Children(tree::Tree, i::Integer)
    Children(tree.child[i], tree.brother)
end


# Compact Clique Tree Data Structures in Sparse Matrix Factorizations
# Pothen and Sun
# Figure 4: The Clique Tree Algorithm 2
function stree(tree::Tree, colcount::AbstractVector, type::Maximal)
    new = Stack{Int}(length(tree))
    parent = Stack{Int}(length(tree))
    ancestor = Stack{Int}(length(tree))
    new_in_clique = Vector{Int}(undef, length(tree))

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


# Compute the fundamental supernode partition.
function stree(tree::Tree, colcount::AbstractVector, type::Fundamental)
    new = Stack{Int}(length(tree))
    parent = Stack{Int}(length(tree))
    ancestor = Stack{Int}(length(tree))
    new_in_clique = Vector{Int}(undef, length(tree))

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


# A Compact Row Storage Scheme for Cholesky Factors Using Elimination Trees
# Liu
# Algorithm 4.2: Elimination Tree by Path Compression.
function etree(upper::SparseMatrixCSC)
    parent = Vector{Int}(undef, size(upper, 1))
    ancestor = Vector{Int}(undef, size(upper, 1))

    for i in axes(upper, 1)
        parent[i] = 0
        ancestor[i] = 0

        for k in view(rowvals(upper), nzrange(upper, i))
            r = k

            while !iszero(ancestor[r]) && ancestor[r] != i
                t = ancestor[r]
                ancestor[r] = i
                r = t
            end

            if iszero(ancestor[r])
                ancestor[r] = i
                parent[r] = i
            end
        end
    end

    Tree(parent)
end


# An Efficient Algorithm to Compute Row and Column Counts for Sparse Cholesky Factorization
# Gilbert, Ng, and Peyton
# Figure 3: Implementation of algorithm to compute row and column counts.
function supcnt(lower::SparseMatrixCSC, tree::Tree, level::AbstractVector=levels(tree), fdesc::AbstractVector=firstdescendants(tree))
    sets = DisjointTrees(size(lower, 1))
    prev_p = zeros(Int, size(lower, 1))
    prev_nbr = zeros(Int, size(lower, 1))
    rc = ones(Int, size(lower, 1))
    wt = ones(Int, size(lower, 1))

    for p in tree[1:end - 1]
        wt[parentindex(tree, p)] = 0
    end

    for p in tree[1:end - 1]
        wt[parentindex(tree, p)] -= 1

        for u in view(rowvals(lower), nzrange(lower, p))
            if fdesc[p] > prev_nbr[u]
                wt[p] += 1
                pp = prev_p[u]

                if iszero(pp)
                    rc[u] += level[p] - level[u]
                else
                    q = find!(sets, pp)
                    rc[u] += level[p] - level[q]
                    wt[q] -= 1
                end

                prev_p[u] = p
            end

            prev_nbr[u] = p
        end

        union!(sets, p, parentindex(tree, p))
    end

    cc = wt

    for p in tree[1:end - 1]
        cc[parentindex(tree, p)] += cc[p]
    end

    rc, cc
end


# Compute a postordering of a rooted tree.
function dfs(tree::Tree)
    head = copy(tree.child)
    order = Index(length(tree))
    stack = Stack{Int}(length(tree))
    push!(stack, tree.root[])

    while !isempty(stack)
        j = last(stack)
        i = head[j]

        if iszero(i)
            push!(order, pop!(stack))
        else
            head[j] = tree.brother[i]
            push!(stack, i)
        end
    end

    order.index
end


# Get the level of every vertex in a topologically ordered rooted tree.
function levels(tree::Tree)
    level = Vector{Int}(undef, length(tree))
    level[end] = 0

    for i in reverse(tree[1:end - 1])
        j = parentindex(tree, i)
        level[i] = level[j] + 1
    end

    level
end


# Get the first descendant of every vertex in a postordered rooted tree.
function firstdescendants(tree::Tree)
    fdesc = collect(tree)

    for i in tree[1:end - 1]
        j = parentindex(tree, i)
        fdesc[j] = min(fdesc[i], fdesc[j])
    end

    fdesc
end


# Permute the vertices of a rooted tree.
function Base.permute!(tree::Tree, permutation::AbstractVector)
    invpermute!(tree, invperm(permutation))
end


# Permute the vertices of a rooted tree.
function Base.invpermute!(tree::Tree, index::AbstractVector)
    fill!(tree.child, 0)

    tree.parent[index] = map(tree.parent) do i
        iszero(i) ? i : index[i]
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


function Base.show(io::IO, ::MIME"text/plain", tree::Tree)
    print_tree(io, IndexNode(tree))
end


##########################
# Indexed Tree Interface #
##########################


function firstchildindex(tree::Tree, i::Integer)
    j = tree.child[i]
    iszero(j) ? nothing : j
end


function AbstractTrees.rootindex(tree::Tree)
    j = tree.root[]
    iszero(j) ? nothing : j
end


function AbstractTrees.parentindex(tree::Tree, i::Integer)
    j = tree.parent[i]
    iszero(j) ? nothing : j
end


function AbstractTrees.nextsiblingindex(tree::Tree, i::Integer)
    j = tree.brother[i]
    iszero(j) ? nothing : j
end


function AbstractTrees.childindices(tree::Tree, i::Integer)
    Children(tree, i)
end


function AbstractTrees.ParentLinks(::Type{IndexNode{Tree, Int}})
    StoredParents()
end


function AbstractTrees.SiblingLinks(::Type{IndexNode{Tree, Int}})
    StoredSiblings()
end


function AbstractTrees.NodeType(::Type{IndexNode{Tree, Int}})
    HasNodeType()
end


function AbstractTrees.nodetype(::Type{IndexNode{Tree, Int}})
    IndexNode{Tree, Int}
end


#################################
# Abstract Unit Range Interface #
#################################


function Base.first(tree::Tree)
    1
end


function Base.last(tree::Tree)
    length(tree.parent)
end
