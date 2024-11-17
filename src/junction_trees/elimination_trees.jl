# An ordered graph (G, σ) equipped with the elimination tree T of its elimination graph.
# Nodes i in T correspond to vertices σ(i) in G.
struct EliminationTree{T <: AbstractTree} <: AbstractTree
    tree::T
    ograph::OrderedGraph
end


# Construct an elimination tree using an elimination algorithm.
function EliminationTree(
        graph::AbstractSymmetricGraph,
        ealg::Union{Order, EliminationAlgorithm}=DEFAULT_ELIMINATION_ALGORITHM)
    EliminationTree(OrderedGraph(graph, ealg))
end


# Construct the elimination tree of an ordered graph.
function EliminationTree(ograph::OrderedGraph)
    EliminationTree(Tree(etree(ograph)), ograph)
end


# Postorder an elimination tree.
function EliminationTree{PostorderTree}(etree::EliminationTree, order::Order)
    EliminationTree(PostorderTree(etree.tree, order), OrderedGraph(etree.ograph, order))
end


# A Compact Row Storage Scheme for Cholesky Factors Using Elimination Trees
# Liu
# Algorithm 4.2: Elimination Tree by Path Compression.
function etree(ograph::OrderedGraph)
    n = nv(graph)
    parent = collect(1:n)
    ancestor = collect(1:n)

    for i in 1:n
        for k in inneighbors(ograph, i)
            r = k

            while ancestor[r] != r && ancestor[r] != i
                t = ancestor[r]
                ancestor[r] = i
                r = t
            end

            if ancestor[r] == r
                ancestor[r] = i
                parent[r] = i
            end
        end
    end

    parent
end


# An Efficient Algorithm to Compute Row and Column Counts for Sparse Cholesky Factorization
# Gilbert, Ng, and Peyton
# Figure 3: Implementation of algorithm to compute row and column counts.
function supcnt(etree::EliminationTree)
    order = postorder(etree)
    index = inverse(order)
    rc, cc = supcnt(EliminationTree{PostorderTree}(etree, order))
    rc[index], cc[index]
end


# An Efficient Algorithm to Compute Row and Column Counts for Sparse Cholesky Factorization
# Gilbert, Ng, and Peyton
# Figure 3: Implementation of algorithm to compute row and column counts.
function supcnt(etree::EliminationTree{PostorderTree})
    n = length(etree)
    
    #### Disjoint Set Union ####
    
    rvert = collect(1:n)
    index = collect(1:n)
    forest = IntDisjointSets(n)
    
    function find(u)
        index[find_root!(forest, u)]
    end
    
    function union(u, v)
        w = max(u, v)
        rvert[w] = root_union!(forest, rvert[u], rvert[v])
        index[rvert[w]] = w
    end
    
    ############################
    
    prev_p = zeros(Int, n)
    prev_nbr = zeros(Int, n)
    rc = ones(Int, n)
    wt = ones(Int, n)

    for u in 1:n - 1
        wt[parentindex(etree, u)] = 0
    end
    
    for p in 1:n - 1
        wt[parentindex(etree, p)] -= 1

        for u in outneighbors(etree, p)
            if firstdescendant(etree, p) > prev_nbr[u]
                wt[p] += 1
                pp = prev_p[u]
                
                if iszero(pp)
                    rc[u] += level(etree, p) - level(etree, u)
                else
                    q = find(pp)
                    rc[u] += level(etree, p) - level(etree, q)
                    wt[q] -= 1
                end
    
                prev_p[u] = p
            end

            prev_nbr[u] = p
        end

        union(p, parentindex(etree, p))
    end

    cc = wt

    for v in 1:n - 1
        cc[parentindex(etree, v)] += cc[v]
    end

    rc, cc
end


# Compute higher degree of every vertex in the elimination graph of
#    (G, σ).
function outdegrees(etree::EliminationTree)
    rc, cc = supcnt(etree)
    cc .- 1
end


# Get the number of nodes in T.
function Base.length(etree::EliminationTree)
    length(etree.tree)
end


# Get the level of a node i.
function level(etree::EliminationTree, i::Integer)
    level(etree.tree, i)
end


# Get the first descendant of a node i.
function firstdescendant(etree::EliminationTree, i::Integer)
    firstdescendant(etree.tree, i)
end


# Determine whether a vertex i is a descendant of a node j.
function isdescendant(etree::EliminationTree, i::Integer, j::Integer)
    isdescendant(etree.tree, i, j)
end


# Get the vertex σ(i).
function order(etree::EliminationTree, i)
    order(etree.ograph, i)
end


# Get the index σ⁻¹(v),
function inverse(etree::EliminationTree, i)
    inverse(etree.ograph, i)
end


##########################
# Indexed Tree Interface #
##########################


function AbstractTrees.rootindex(etree::EliminationTree)
    rootindex(etree.tree)
end


function AbstractTrees.parentindex(etree::EliminationTree, i::Integer)
    parentindex(etree.tree, i)
end


function AbstractTrees.childindices(etree::EliminationTree, i::Integer)
    childindices(etree.tree, i)
end


function AbstractTrees.NodeType(::Type{IndexNode{EliminationTree, Int}})
    HasNodeType()
end


function AbstractTrees.nodetype(::Type{IndexNode{EliminationTree{T}, Int}}) where T
    IndexNode{EliminationTree{T}, Int}
end


############################
# Abstract Graph Interface #
############################


function BasicGraphs.inneighbors(etree::EliminationTree, i::Integer)
    inneighbors(etree.ograph, i)
end


function BasicGraphs.outneighbors(etree::EliminationTree, i::Integer)
    outneighbors(etree.ograph, i)
end
