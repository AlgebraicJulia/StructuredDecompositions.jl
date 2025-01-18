"""
    Tree <: AbstractUnitRange{Int}

A rooted tree. This type implements the
[indexed tree interface](https://juliacollections.github.io/AbstractTrees.jl/stable/#The-Indexed-Tree-Interface).
"""
struct Tree <: AbstractUnitRange{Int}
    parent::Vector{Int}  # vector of parents
    root::Scalar{Int}    # root

    # left-child right-sibling tree
    child::Vector{Int}   # vector of left-children
    brother::Vector{Int} # vector of right-siblings

    function Tree(parent::AbstractVector, root::AbstractScalar, child::AbstractVector, brother::AbstractVector)
        # validate arguments
        eachindex(parent) != eachindex(brother) && throwthrow(ArgumentError("eachindex(parent) != eachindex(brother)"))
        eachindex(parent) != eachindex(child) && throwthrow(ArgumentError("eachindex(parent) != eachindex(child)"))
        root[] ∉ eachindex(parent) && throw(ArgumentError("root[] ∉ eachindex(parent)"))

        # construct tree
        new(parent, root, child, brother)
    end
end


# Construct a rooted tree.
function Tree(parent::AbstractVector)
    # validate argument
    isempty(parent) && throw(ArgumentError("isempty(parent)"))

    # construct tree
    root = zeros(Int)
    child = zeros(Int, length(parent))
    brother = Vector{Int}(undef, length(parent))

    for (i, j) in ireverse(enumerate(parent))
        if iszero(j)
            root[] = i
            brother[i] = 0
        else
            brother[i] = child[j]
            child[j] = i
        end
    end

    Tree(parent, root, child, brother)
end


function Tree(tree::Tree)
    Tree(tree.parent, tree.root, tree.child, tree.brother)
end


# Construct a rooted tree.
function Tree(tree::Tree, root::Integer)
    # validate arguments
    root ∉ tree && throw(ArgumentError("root ∉ tree"))

    # construct tree
    parent = copy(tree.parent)
    i = root
    j = 0

    while !isnothing(i)
        parent[i] = j
        j = i
        i = parentindex(tree, j)
    end

    Tree(parent)
end


# Construct a rooted tree.
function Tree(matrix::SparseMatrixCSC, root::Integer)
    # validate arguments
    size(matrix, 1) != size(matrix, 2) && throw(ArgumentError("size(matrix, 1) != size(matrix, 2)"))
    root ∉ axes(matrix, 2) && throw(ArgumentError("root ∉ axes(matrix, 2)"))

    # construct tree
    parent = zeros(Int, size(matrix, 2))
    parent[root] = 0

    seen = zeros(Bool, size(matrix, 2))
    seen[root] = true

    stack = sizehint!(Int[], size(matrix, 2))
    push!(stack, root)

    while !isempty(stack)
        v = last(stack)
        u = 0

        for n in @view rowvals(matrix)[nzrange(matrix, v)]
            if !seen[n]
                u = n
                break
            end
        end

        if iszero(u)
            pop!(stack)
        else
            seen[u] = true
            push!(stack, u)
            parent[u] = v
        end
    end

    Tree(parent)
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

Construct a [tree-depth decomposition](https://en.wikipedia.org/wiki/Trémaux_tree) of a connected simple graph.
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

Compute the depth of a topologically ordered tree.
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

Compute an upper bound to the [tree-depth](https://en.wikipedia.org/wiki/Tree-depth) of a connected simple graph.
See [`junctiontree!`](@ref) for the meaning of `alg`.
"""
function treedepth!(matrix::AbstractMatrix; alg::PermutationOrAlgorithm=DEFAULT_ELIMINATION_ALGORITHM)
    treedepth(eliminationtree!(matrix, alg))
end


# A Compact Row Storage Scheme for Cholesky Factors Using Elimination Trees
# Liu
# Algorithm 4.2: Elimination Tree by Path Compression.
#
# Construct the elimination tree of an ordered graph.
# The complexity is O(mlogn), where m = |E| and n = |V|.
function etree(upper::SparseMatrixCSC)
    # validate argument
    size(upper, 1) != size(upper, 2) && throw(ArgumentError("size(upper, 1) != size(upper, 2)"))

    # run algorithm
    parent = Vector{Int}(undef, size(upper, 2))
    ancestor = Vector{Int}(undef, size(upper, 2))

    for i in axes(upper, 2)
        parent[i] = 0
        ancestor[i] = 0

        for k in @view rowvals(upper)[nzrange(upper, i)]
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
#
# Compute the lower and higher degrees of the monotone transitive extesion of an ordered graph.
# The complexity is O(mα(m, n)), where m = |E|, n = |V|, and α is the inverse Ackermann function.
function supcnt(lower::SparseMatrixCSC, tree::Tree, index::AbstractVector=dfs(tree))
    # validate arguments
    tree != axes(lower, 1) && throw(ArgumentError("tree != axes(lower, 1)"))
    tree != axes(lower, 2) && throw(ArgumentError("tree != axes(lower, 2)"))
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
    wt = ones(Int, length(tree))

    for p in tree[1:end - 1]
        wt[parentindex(tree, p)] = 0
    end

    for p in @view invperm(index)[1:end - 1]
        wt[parentindex(tree, p)] -= 1

        for u in @view rowvals(lower)[nzrange(lower, p)]
            if iszero(prev_nbr[u]) || lt(order, prev_nbr[u], fdesc[p])
                wt[p] += 1
                pp = prev_p[u]

                if iszero(pp)
                    rc[u] += level[p] - level[u]
                else
                    q = find(pp)
                    rc[u] += level[p] - level[q]
                    wt[q] -= 1
                end

                prev_p[u] = p
            end

            prev_nbr[u] = p
        end

        union(p, parentindex(tree, p))
    end

    cc = wt

    for p in tree[1:end - 1]
        cc[parentindex(tree, p)] += cc[p]
    end

    rc, cc
end


# Compute a postordering of a rooted tree.
function dfs(tree::Tree)
    stack = Vector{Int}(undef, length(tree))
    index = Vector{Int}(undef, length(tree))
    head = copy(tree.child)
    p = 1
    q = 1
    
    stack[p] = rootindex(tree)

    while !iszero(p)
        j = stack[p]
        i = head[j]

        if iszero(i)
            index[stack[p]] = q
            p -= 1
            q += 1
        else
            head[j] = tree.brother[i]
            p += 1
            stack[p] = i
        end
    end

    index
end


# Get the level of every vertex in a topologically ordered tree.
function levels(tree::Tree)
    level = Vector{Int}(undef, length(tree))
    level[end] = 0

    for i in reverse(tree[1:end - 1])
        j = parentindex(tree, i)
        level[i] = level[j] + 1
    end

    level
end


# Get the first descendant of every vertex in a topologically ordered tree.
function firstdescendants(tree::Tree, order::Ordering=ForwardOrdering())
    fdesc = collect(tree)

    for i in tree[1:end - 1]
        j = parentindex(tree, i)
        u = fdesc[i]
        v = fdesc[j]

        if lt(order, v, u)
            u, v = v, u
        end

        fdesc[j] = u
    end

    fdesc
end


# Permute the vertices of a rooted tree.
function Base.invpermute!(tree::Tree, index::AbstractVector)
    # validate arguments
    tree != eachindex(index) && throw(ArgumentError("tree != eachindex(index)"))

    # run algorithm
    fill!(tree.child, 0)

    tree.parent[index] = map(tree.parent) do i
        iszero(i) ? i : index[i]
    end

    for (i, j) in ireverse(enumerate(tree.parent))
        if iszero(j)
            tree.root[] = i
            tree.brother[i] = 0
        else
            tree.brother[i] = tree.child[j]
            tree.child[j] = i
        end
    end

    tree
end


function Base.show(io::IO, ::MIME"text/plain", tree::Tree)
    println(io, "$(length(tree))-element Tree:")
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
    tree.root[]
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
    SinglyLinkedList(view(tree.child, i), tree.brother)
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
