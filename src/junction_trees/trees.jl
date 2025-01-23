"""
    Tree <: AbstractUnitRange{Int}

A rooted forest. This type implements the
[indexed tree interface](https://juliacollections.github.io/AbstractTrees.jl/stable/#The-Indexed-Tree-Interface).
"""
struct Tree <: AbstractUnitRange{Int}
    parent::Vector{Int}  # vector of parents
    root::Scalar{Int}    # root
    child::Vector{Int}   # vector of left-children
    brother::Vector{Int} # vector of right-siblings

    function Tree(parent::AbstractVector)
        root = Scalar{Int}(undef)
        child = Vector{Int}(undef, length(parent))
        brother = Vector{Int}(undef, length(parent))
        lcrs!(new(parent, root, child, brother))
    end

    function Tree(parent::AbstractVector, root::AbstractScalar, child::AbstractVector, brother::AbstractVector)
        # validate arguments
        eachindex(parent) != eachindex(child) && throw(ArgumentError("eachindex(parent) != eachindex(child)"))
        eachindex(parent) != eachindex(brother) && throw(ArgumentError("eachindex(parent) != eachindex(brother)"))
        isempty(parent) && !iszero(root[]) && throw(ArgumentError("isempty(parent) && !iszero(root[])"))
        !isempty(parent) && root[] ∉ eachindex(parent) && throw(ArgumentError("!isempty(parent) && root[] ∉ eachindex(parent)"))

        # construct tree
        new(parent, root, child, brother)
    end
end


function Tree(tree::Tree)
    Tree(tree.parent, tree.root, tree.child, tree.brother)
end


"""
    eliminationtree(matrix::AbstractMatrix;
        alg::PermutationOrAlgorithm=DEFAULT_ELIMINATION_ALGORITHM)

A non-mutating version of [`eliminationtree!`](@ref).
"""
function eliminationtree(matrix::AbstractMatrix; alg::PermutationOrAlgorithm=DEFAULT_ELIMINATION_ALGORITHM)
    eliminationtree!(sparse(matrix); alg)
end


"""
    eliminationtree!(matrix::AbstractMatrix;
        alg::PermutationOrAlgorithm=DEFAULT_ELIMINATION_ALGORITHM)

Construct a [tree-depth decomposition](https://en.wikipedia.org/wiki/Trémaux_tree) of a simple graph.
See [`junctiontree!`](@ref) for the meaning of `alg`.
```julia
julia> using StructuredDecompositions

julia> graph = [
           0 1 1 0 0 0 0 0
           1 0 1 0 0 1 0 0
           1 1 0 1 1 0 0 0
           0 0 1 0 1 0 0 0
           0 0 1 1 0 0 1 1
           0 1 0 0 0 0 1 0
           0 0 0 0 1 1 0 1
           0 0 0 0 1 0 1 0
       ];

julia> label, tree = eliminationtree(graph);

julia> tree
8-element Tree:
8
└─ 7
   ├─ 5
   │  ├─ 1
   │  └─ 4
   │     └─ 3
   │        └─ 2
   └─ 6
```
"""
function eliminationtree!(matrix::SparseMatrixCSC; alg::PermutationOrAlgorithm=DEFAULT_ELIMINATION_ALGORITHM)
    label, tree, upper, cache = eliminationtree!(matrix, alg)
    label, tree
end


function eliminationtree!(matrix::SparseMatrixCSC, alg::PermutationOrAlgorithm)
    label, index = permutation(matrix, alg)
    cache = triu(matrix, 1)
    upper = sympermute!(matrix, cache, index)
    label, etree(upper), upper, cache
end


"""
    treedepth(tree::Tree)

Compute the depth of a topologically ordered forest.
"""
function treedepth(tree::Tree)
    maximum(levels(tree))
end


"""
    treedepth(matrix::AbstractMatrix;
        alg::PermutationOrAlgorithm=DEFAULT_ELIMINATION_ALGORITHM)

A non-mutating version of [`treedepth!`](@ref).
"""
function treedepth(matrix::AbstractMatrix; alg::PermutationOrAlgorithm=DEFAULT_ELIMINATION_ALGORITHM)
    treedepth!(sparse(matrix); alg)
end


"""
    treedepth!(matrix::SparseMatrixCSC;
        alg::PermutationOrAlgorithm=DEFAULT_ELIMINATION_ALGORITHM)

Compute an upper bound to the [tree-depth](https://en.wikipedia.org/wiki/Tree-depth) of a simple graph.
See [`junctiontree!`](@ref) for the meaning of `alg`.
"""
function treedepth!(matrix::AbstractMatrix; alg::PermutationOrAlgorithm=DEFAULT_ELIMINATION_ALGORITHM)
    treedepth(eliminationtree!(matrix, alg))
end


function etree(upper::SparseMatrixCSC)
    # validate argument
    size(upper, 1) != size(upper, 2) && throw(ArgumentError("size(upper, 1) != size(upper, 2)"))

    # run algorithm
    etree(size(upper, 2)) do j
        @view rowvals(upper)[nzrange(upper, j)]
    end
end


# A Compact Row Storage Scheme for Cholesky Factors Using Elimination Trees
# Liu
# Algorithm 4.2: Elimination Tree by Path Compression.
#
# Construct the elimination tree of an ordered graph.
# The complexity is O(mlogn), where m = |E| and n = |V|.
function etree(neighbors::Function, n::Integer)
    # validate arguments
    n < 0 && throw(ArgumentError("n < 0"))

    # run algorithm
    parent = Vector{Int}(undef, n)
    ancestor = Vector{Int}(undef, n)

    for i in 1:n
        parent[i] = 0
        ancestor[i] = 0

        for k in neighbors(i)
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


function supcnt(lower::SparseMatrixCSC, tree::Tree, index::AbstractVector=dfs(tree))
    # validate arguments
    tree != axes(lower, 1) && throw(ArgumentError("tree != axes(lower, 1)"))
    tree != axes(lower, 2) && throw(ArgumentError("tree != axes(lower, 2)"))

    # run algorithm
    supcnt(tree, index) do j
        @view rowvals(lower)[nzrange(lower, j)]
    end
end


# An Efficient Algorithm to Compute Row and Column Counts for Sparse Cholesky Factorization
# Gilbert, Ng, and Peyton
# Figure 3: Implementation of algorithm to compute row and column counts.
#
# Compute the lower and higher degrees of the monotone transitive extesion of an ordered graph.
# The complexity is O(mα(m, n)), where m = |E|, n = |V|, and α is the inverse Ackermann function.
function supcnt(neighbors::Function, tree::Tree, index::AbstractVector=dfs(tree))
    # validate arguments
    tree != eachindex(index) && throw(ArgumentError("tree != eachindex(index)"))

    # find first descendants and levels
    order = Perm(ForwardOrdering(), index)
    fdesc = firstdescendants(tree, order)
    level = levels(tree)

    # construct disjoint set forest
    sets = IntDisjointSets(length(tree))
    root = collect(tree)
    repr = collect(tree)
    
    function find(u)
        v = find_root!(sets, u)
        repr[v]
    end
    
    function union(u, v)
        root[v] = root_union!(sets, root[u], root[v])
        repr[root[v]] = v
    end
    
    # run algorithm
    prev_p = zeros(Int, length(tree))
    prev_nbr = zeros(Int, length(tree))
    rc = ones(Int, length(tree))
    cc = ones(Int, length(tree))

    for p in tree
        r = parentindex(tree, p)
        
        if !isnothing(r)
            cc[r] = 0
        end
    end

    for p in invperm(index)
        r = parentindex(tree, p)

        for u in neighbors(p)
            if iszero(prev_nbr[u]) || lt(order, prev_nbr[u], fdesc[p])
                cc[p] += 1
                pp = prev_p[u]

                if iszero(pp)
                    rc[u] += level[p] - level[u]
                else
                    q = find(pp)
                    rc[u] += level[p] - level[q]
                    cc[q] -= 1
                end

                prev_p[u] = p
            end

            prev_nbr[u] = p
        end

        if !isnothing(r)
            cc[r] -= 1
            union(p, r)
        end
    end

    for p in tree
        r = parentindex(tree, p)

        if !isnothing(r)
            cc[r] += cc[p]
        end
    end

    rc, cc
end


# Compute a postordering of a forest.
function dfs(tree::Tree)
    # construct disjoint sets data structure
    child = copy(tree.child)

    function set(i)
        head = @view child[i]
        SinglyLinkedList(head, tree.brother)
    end
 
    # run algorithm
    index = Vector{Int}(undef, length(tree))
    stack = sizehint!(Int[], length(tree))
    append!(stack, rootindices(tree))

    for i in tree
        j = pop!(stack)

        while !isempty(set(j))
            push!(stack, j)
            j = popfirst!(set(j))
        end

        index[j] = i
    end

    index
end


# Get the level of every vertex in a topologically ordered tree.
function levels(tree::Tree)
    level = Vector{Int}(undef, length(tree))

    for i in reverse(tree)
        j = parentindex(tree, i)
        level[i] = isnothing(j) ? 0 : level[j] + 1
    end

    level
end


# Get the first descendant of every vertex in a topologically ordered forest.
function firstdescendants(tree::Tree, order::Ordering=ForwardOrdering())
    fdesc = collect(tree)

    for j in tree
        for i in childindices(tree, j)
            u = fdesc[i]
            v = fdesc[j]

            if lt(order, v, u)
                u, v = v, u
            end

            fdesc[j] = u
        end
    end

    fdesc
end


# Compute the `root`, `child`, and `brother` fields of a forest.
function lcrs!(tree::Tree)
    fill!(tree.root, 0)
    fill!(tree.child, 0)

    for i in reverse(tree)
        j = parentindex(tree, i)

        if isnothing(j)
            tree.brother[i] = tree.root[]
            tree.root[] = i
        else
            tree.brother[i] = tree.child[j]
            tree.child[j] = i
        end
    end

    tree
end


# Make the node `j` a root.
function setrootindex!(tree::Tree, j::Integer)
    # validate arguments
    j ∉ tree && throw(ArgumentError("j ∉ tree"))

    # run algorithm
    i = 0

    while !iszero(j)
        k = tree.parent[j]
        tree.parent[j] = i
        i, j = j, k
    end

    lcrs!(tree)
end


# Permute the vertices of a forest.
function Base.invpermute!(tree::Tree, index::AbstractVector)
    # validate arguments
    tree != eachindex(index) && throw(ArgumentError("tree != eachindex(index)"))

    # run algorithm
    tree.parent[index] = map(tree.parent) do i
        iszero(i) ? 0 : index[i]
    end

    lcrs!(tree)
end


##########################
# Indexed Tree Interface #
##########################


function AbstractTrees.rootindex(tree::Tree)
    j = tree.root[]
    iszero(j) ? nothing : j
end


function AbstractTrees.parentindex(tree::Tree, i::Integer)
    j = tree.parent[i]
    iszero(j) ? nothing : j
end


function firstchildindex(tree::Tree, i::Integer)
    j = tree.child[i]
    iszero(j) ? nothing : j
end


function AbstractTrees.nextsiblingindex(tree::Tree, i::Integer)
    j = tree.brother[i]
    iszero(j) ? nothing : j
end


function rootindices(tree::Tree)
     SinglyLinkedList(tree.root, tree.brother)
end


function AbstractTrees.childindices(tree::Tree, i::Integer)
    head = @view tree.child[i]
    SinglyLinkedList(head, tree.brother)
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
