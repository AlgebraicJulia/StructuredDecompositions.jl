"""
    Tree{V <: Integer} <: AbstractUnitRange{V}

    Tree(tree::AbstractTree)

    Tree{V}(tree::AbstractTree) where V

A rooted forest. This type implements the
[indexed tree interface](https://juliacollections.github.io/AbstractTrees.jl/stable/#The-Indexed-Tree-Interface).
"""
struct Tree{V <: Integer} <: AbstractUnitRange{V}
    parent::Vector{V}  # vector of parents
    root::Scalar{V}    # root
    child::Vector{V}   # vector of left-children
    brother::Vector{V} # vector of right-siblings

    function Tree{V}(parent::AbstractVector) where V
        root = Scalar{V}(undef)
        child = Vector{V}(undef, length(parent))
        brother = Vector{V}(undef, length(parent))
        tree = new{V}(parent, root, child, brother)
        lcrs!(tree)
    end

    function Tree{V}(parent::AbstractVector, root::AbstractScalar, child::AbstractVector, brother::AbstractVector) where V
        # validate arguments
        eachindex(parent) != eachindex(child) && throw(ArgumentError("eachindex(parent) != eachindex(child)"))
        eachindex(parent) != eachindex(brother) && throw(ArgumentError("eachindex(parent) != eachindex(brother)"))
        isempty(parent) && !iszero(root[]) && throw(ArgumentError("isempty(parent) && !iszero(root[])"))
        !isempty(parent) && root[] ∉ eachindex(parent) && throw(ArgumentError("!isempty(parent) && root[] ∉ eachindex(parent)"))

        # construct tree
        new{V}(parent, root, child, brother)
    end
end


function Tree(parent::AbstractVector{V}) where V
    Tree{V}(parent)
end


function Tree(parent::AbstractVector{V}, root::AbstractScalar{V}, child::AbstractVector{V}, brother::AbstractVector{V}) where V
    Tree{I}(parent, root, child, brother)
end


function Tree(tree::Tree)
    Tree(tree.parent, tree.root, tree.child, tree.brother)
end


function Tree{V}(tree::Tree) where V
    Tree{V}(tree.parent, tree.root, tree.child, tree.brother)
end


"""
    eliminationtree(graph;
        alg::PermutationOrAlgorithm=DEFAULT_ELIMINATION_ALGORITHM)

Construct a [tree-depth decomposition](https://en.wikipedia.org/wiki/Trémaux_tree) of a simple graph.
```julia
julia> using SparseArrays, StructuredDecompositions

julia> graph = sparse([
           0 1 1 0 0 0 0 0
           1 0 1 0 0 1 0 0
           1 1 0 1 1 0 0 0
           0 0 1 0 1 0 0 0
           0 0 1 1 0 0 1 1
           0 1 0 0 0 0 1 0
           0 0 0 0 1 1 0 1
           0 0 0 0 1 0 1 0
       ]);

julia> label, tree = eliminationtree(graph);

julia> tree
8-element Tree{Int64}:
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
function eliminationtree(graph; alg::PermutationOrAlgorithm=DEFAULT_ELIMINATION_ALGORITHM)
    label, tree, upper = eliminationtree(graph, alg)
    label, tree
end


function eliminationtree(graph, alg::PermutationOrAlgorithm)
    label, index = permutation(graph, alg)
    upper = sympermute(graph, index, ForwardOrdering())
    label, etree(upper), upper
end


# A Compact Row Storage Scheme for Cholesky Factors Using Elimination Trees
# Liu
# Algorithm 4.2: Elimination Tree by Path Compression.
#
# Construct the elimination tree of an ordered graph.
# The complexity is O(mlogn), where m = |E| and n = |V|.
function etree(upper::Graph{V}) where V
    parent = Vector{V}(undef, nv(upper))
    ancestor = Vector{V}(undef, nv(upper))

    for i in vertices(upper)
        parent[i] = 0
        ancestor[i] = 0

        for k in neighbors(upper, i)
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
function supcnt(lower::Graph{V}, tree::Tree{V}) where V
    # validate arguments
    vertices(lower) != tree && throw(ArgumentError("vertices(lower) != tree"))

    # find postordering, first descendants, and levels
    index = postorder(tree)
    order = Perm(ForwardOrdering(), index)
    fdesc = firstdescendants(tree, order)
    level = levels(tree)

    # construct disjoint set forest
    sets = IntDisjointSets{V}(length(tree))
    root::Vector{V} = tree
    repr::Vector{V} = tree
    
    function find(u)
        v = find_root!(sets, u)
        repr[v]
    end
    
    function union(u, v)
        root[v] = root_union!(sets, root[u], root[v])
        repr[root[v]] = v
    end
    
    # run algorithm
    prev_p = zeros(V, length(tree))
    prev_nbr = zeros(V, length(tree))
    rc = ones(V, length(tree))
    wt = ones(Int, length(tree))

    for p in tree
        r = parentindex(tree, p)
        
        if !isnothing(r)
            wt[r] = 0
        end
    end

    for p in invperm(index)
        r = parentindex(tree, p)

        for u in neighbors(lower, p)
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

        if !isnothing(r)
            wt[r] -= 1
            union(p, r)
        end
    end

    for p in tree
        r = parentindex(tree, p)

        if !isnothing(r)
            wt[r] += wt[p]
        end
    end

    cc::Vector{V} = wt
    rc, cc 
end


# Compute a postordering of a forest.
function postorder(tree::Tree{V}) where V
    index = Vector{V}(undef, length(tree))

    dfs(tree) do i, j
        index[j] = i
    end

    index
end


# Perform a depth-first search of a forest.
function dfs(f::Function, tree::Tree{V}) where V
    # construct disjoint sets data structure
    child = copy(tree.child)

    function set(i)
        head = @view child[i]
        SinglyLinkedList(head, tree.brother)
    end

    # construct stack 
    stack = sizehint!(V[], length(tree))

    # run algorithm
    append!(stack, rootindices(tree))

    for i in tree
        j = pop!(stack)

        while !isempty(set(j))
            push!(stack, j)
            j = popfirst!(set(j))
        end

        f(i, j)
    end
end


# Get the level of every vertex in a topologically ordered tree.
function levels(tree::Tree{V}) where V
    level = Vector{V}(undef, length(tree))
    
    for i in reverse(tree)
        j = parentindex(tree, i)
        level[i] = isnothing(j) ? 0 : level[j] + 1
    end

    level
end


# Get the first descendant of every vertex in a topologically ordered forest.
function firstdescendants(tree::Tree{V}, order::Ordering=ForwardOrdering()) where V
    fdesc = Vector{V}(undef, length(tree))

    for j in tree
        v = j

        for i in childindices(tree, j)
            u = fdesc[i]

            if lt(order, u, v)
                v = u
            end
        end

        fdesc[j] = v
    end

    fdesc
end


# Compute the `root`, `child`, and `brother` fields of a forest.
function lcrs!(tree::Tree{V}) where V
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


# Permute the vertices of a forest.
function Base.invpermute!(tree::Tree{V}, index::AbstractVector{V}) where V
    # validate arguments
    tree != eachindex(index) && throw(ArgumentError("tree != eachindex(index)"))

    # run algorithm
    tree.parent[index] = map(tree.parent) do i
        iszero(i) ? i : index[i]
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


function ancestorindices(tree::Tree, i::Integer)
    head = @view tree.parent[i]
    SinglyLinkedList(head, tree.parent)
end


function setrootindex!(tree::Tree{V}, i::Integer) where V
    # validate arguments
    i ∉ tree && throw(ArgumentError("i ∉ tree"))

    # run algorithm
    j::V = 0
    k::V = i

    for l in ancestorindices(tree, k)
        tree.parent[k] = j
        j, k = k, l
    end

    tree.parent[k] = j
    lcrs!(tree)
end


#################################
# Abstract Unit Range Interface #
#################################


function Base.first(tree::Tree{V}) where V
    one(V)
end


function Base.last(tree::Tree{V}) where V
    convert(V, length(tree.parent))
end
