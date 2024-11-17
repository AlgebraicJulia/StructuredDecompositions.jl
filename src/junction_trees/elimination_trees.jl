# An ordered graph (G, σ) equipped with the elimination tree T of its elimination graph.
# Nodes i in T correspond to vertices σ(i) in G.
struct EliminationTree{T <: AbstractTree}
    tree::T              # elimination tree
    ograph::OrderedGraph # ordered graph
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
    n = nv(ograph)
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
    order = postorder(etree.tree)
    index = inverse(order)
    rc, cc = supcnt(EliminationTree{PostorderTree}(etree, order))
    rc[index], cc[index]
end


# An Efficient Algorithm to Compute Row and Column Counts for Sparse Cholesky Factorization
# Gilbert, Ng, and Peyton
# Figure 3: Implementation of algorithm to compute row and column counts.
function supcnt(etree::EliminationTree{PostorderTree})
    n = length(etree.tree)
    
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
        wt[parentindex(etree.tree, u)] = 0
    end
    
    for p in 1:n - 1
        wt[parentindex(etree.tree, p)] -= 1

        for u in outneighbors(etree.ograph, p)
            if firstdescendant(etree.tree, p) > prev_nbr[u]
                wt[p] += 1
                pp = prev_p[u]
                
                if iszero(pp)
                    rc[u] += level(etree.tree, p) - level(etree.tree, u)
                else
                    q = find(pp)
                    rc[u] += level(etree.tree, p) - level(etree.tree, q)
                    wt[q] -= 1
                end
    
                prev_p[u] = p
            end

            prev_nbr[u] = p
        end

        union(p, parentindex(etree.tree, p))
    end

    cc = wt

    for v in 1:n - 1
        cc[parentindex(etree.tree, v)] += cc[v]
    end

    rc, cc
end


# Compute higher degree of every vertex in the elimination graph of
#    (G, σ).
function outdegrees(etree::EliminationTree)
    rc, cc = supcnt(etree)
    cc .- 1
end
