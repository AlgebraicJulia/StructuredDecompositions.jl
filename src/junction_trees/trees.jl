"""
    Tree{I <: Integer} <: AbstractUnitRange{I}

    Tree(tree::AbstractTree)

    Tree{I}(tree::AbstractTree) where I

A rooted forest. This type implements the
[indexed tree interface](https://juliacollections.github.io/AbstractTrees.jl/stable/#The-Indexed-Tree-Interface).
"""
struct Tree{I <: Integer} <: AbstractUnitRange{I}
    parent::Vector{I}  # vector of parents
    root::Scalar{I}    # root
    child::Vector{I}   # vector of left-children
    brother::Vector{I} # vector of right-siblings

    function Tree{I}(parent::AbstractVector) where I
        root = Scalar{I}(undef)
        child = Vector{I}(undef, length(parent))
        brother = Vector{I}(undef, length(parent))
        tree = new{I}(parent, root, child, brother)
        lcrs!(tree)
    end

    function Tree{I}(parent::AbstractVector, root::AbstractScalar, child::AbstractVector, brother::AbstractVector) where I
        # validate arguments
        eachindex(parent) != eachindex(child) && throw(ArgumentError("eachindex(parent) != eachindex(child)"))
        eachindex(parent) != eachindex(brother) && throw(ArgumentError("eachindex(parent) != eachindex(brother)"))
        isempty(parent) && !iszero(root[]) && throw(ArgumentError("isempty(parent) && !iszero(root[])"))
        !isempty(parent) && root[] ∉ eachindex(parent) && throw(ArgumentError("!isempty(parent) && root[] ∉ eachindex(parent)"))

        # construct tree
        new{I}(parent, root, child, brother)
    end
end


function Tree(parent::AbstractVector{I}) where I
    Tree{I}(parent)
end


function Tree(parent::AbstractVector{I}, root::AbstractScalar{I}, child::AbstractVector{I}, brother::AbstractVector{I}) where I
    Tree{I}(parent, root, child, brother)
end


function Tree(tree::Tree)
    Tree(tree.parent, tree.root, tree.child, tree.brother)
end


function Tree{I}(tree::Tree) where I
    Tree{I}(tree.parent, tree.root, tree.child, tree.brother)
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


function etree(upper::SparseMatrixCSC{<:Any, I}) where I
    # validate argument
    size(upper, 1) != size(upper, 2) && throw(ArgumentError("size(upper, 1) != size(upper, 2)"))

    # run algorithm
    vertices::OneTo{I} = axes(upper, 2)

    etree(vertices) do j
        @view rowvals(upper)[nzrange(upper, j)]
    end
end


# A Compact Row Storage Scheme for Cholesky Factors Using Elimination Trees
# Liu
# Algorithm 4.2: Elimination Tree by Path Compression.
#
# Construct the elimination tree of an ordered graph.
# The complexity is O(mlogn), where m = |E| and n = |V|.
function etree(neighbors::Function, vertices::AbstractVector{I}) where I
    parent = Vector{I}(undef, length(vertices))
    ancestor = Vector{I}(undef, length(vertices))

    for i in vertices
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


function supcnt(lower::SparseMatrixCSC{<:Any, I}, tree::Tree{I}) where I
    # validate arguments
    tree != axes(lower, 1) && throw(ArgumentError("tree != axes(lower, 1)"))
    tree != axes(lower, 2) && throw(ArgumentError("tree != axes(lower, 2)"))

    # run algorithm
    supcnt(tree) do j
        @view rowvals(lower)[nzrange(lower, j)]
    end
end


# An Efficient Algorithm to Compute Row and Column Counts for Sparse Cholesky Factorization
# Gilbert, Ng, and Peyton
# Figure 3: Implementation of algorithm to compute row and column counts.
#
# Compute the lower and higher degrees of the monotone transitive extesion of an ordered graph.
# The complexity is O(mα(m, n)), where m = |E|, n = |V|, and α is the inverse Ackermann function.
function supcnt(neighbors::Function, tree::Tree{I}) where I
    # find postordering, first descendants, and levels
    index = postorder(tree)
    order = Perm(ForwardOrdering(), index)
    fdesc = firstdescendants(tree, order)
    level = levels(tree)

    # construct disjoint set forest
    sets = IntDisjointSets{I}(length(tree))
    root::Vector{I} = tree
    repr::Vector{I} = tree
    
    function find(u)
        v = find_root!(sets, u)
        repr[v]
    end
    
    function union(u, v)
        root[v] = root_union!(sets, root[u], root[v])
        repr[root[v]] = v
    end
    
    # run algorithm
    prev_p = zeros(I, length(tree))
    prev_nbr = zeros(I, length(tree))
    rc = ones(I, length(tree))
    cc = ones(I, length(tree))

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
function postorder(tree::Tree{I}) where I
    index = Vector{I}(undef, length(tree))

    dfs(tree) do i, j
        index[j] = i
    end

    index
end


# Perform a depth-first search of a forest.
function dfs(f::Function, tree::Tree{I}) where I
    # construct disjoint sets data structure
    child = copy(tree.child)

    function set(i)
        head = @view child[i]
        SinglyLinkedList(head, tree.brother)
    end

    # construct stack 
    stack = sizehint!(I[], length(tree))

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
function levels(tree::Tree{I}) where I
    level = Vector{I}(undef, length(tree))
    
    for i in reverse(tree)
        j = parentindex(tree, i)
        level[i] = isnothing(j) ? 0 : level[j] + 1
    end

    level
end


# Get the first descendant of every vertex in a topologically ordered forest.
function firstdescendants(tree::Tree{I}, order::Ordering=ForwardOrdering()) where I
    fdesc = Vector{I}(undef, length(tree))

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
function lcrs!(tree::Tree{I}) where I
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
function Base.invpermute!(tree::Tree{I}, index::AbstractVector{I}) where I
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


function setrootindex!(tree::Tree{I}, i::Integer) where I
    # validate arguments
    i ∉ tree && throw(ArgumentError("i ∉ tree"))

    # run algorithm
    j::I = 0
    k::I = i

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


function Base.first(tree::Tree{I}) where I
    one(I)
end


function Base.last(tree::Tree{I}) where I
    convert(I, length(tree.parent))
end
