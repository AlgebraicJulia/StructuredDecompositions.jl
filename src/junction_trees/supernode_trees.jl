# An ordered graph (G, σ) equipped with a supernodal elimination tree T.
struct SupernodeTree
    tree::PostorderTree         # supernodal elimination tree
    graph::OrderedGraph         # ordered graph
    representative::Vector{Int} # representative vertex
    cardinality::Vector{Int}    # supernode cardinality
    ancestor::Vector{Int}       # first ancestor
    degree::Vector{Int}         # higher degrees
end


# Construct a supernodal elimination tree using an elimination algorithm.
function SupernodeTree(
        graph::AbstractSymmetricGraph,
        ealg::Union{Order, EliminationAlgorithm}=DEFAULT_ELIMINATION_ALGORITHM,
        stype::SupernodeType=DEFAULT_SUPERNODE_TYPE)

        SupernodeTree(EliminationTree(graph, ealg), stype)
end
    

# Construct a supernodal elimination tree.
function SupernodeTree(etree::EliminationTree, stype::SupernodeType=DEFAULT_SUPERNODE_TYPE)
    degree = outdegrees(etree)
    snode, parent, ancestor = stree(etree, degree, stype)
    tree = Tree(parent)

    treeorder = postorder(tree)
    graphorder = Order(vcat(snode[treeorder]...))

    tree = PostorderTree(tree, treeorder)
    graph = OrderedGraph(etree.graph, graphorder)

    permute!(snode, treeorder)
    permute!(ancestor, treeorder)
    permute!(degree, graphorder)

    n = length(tree)
    representative = zeros(Int, n)
    cardinality = zeros(Int, n)

    for i in 1:n - 1
        representative[i] = inverse(graphorder, snode[i][1])
        cardinality[i] = length(snode[i])
        ancestor[i] = inverse(graphorder, ancestor[i])
    end

    representative[n] = inverse(graphorder, snode[n][1])
    cardinality[n] = length(snode[n])

    SupernodeTree(tree, graph, representative, cardinality, ancestor, degree)
end


# Construct an elimination graph.
function eliminationgraph(stree::SupernodeTree)
    graph = deepcopy(stree.graph)
    n = length(stree.tree)

    for i in 1:n - 1
        for u in supernode(stree, i)[1:end - 1]
            v = u + 1

            for w in outneighbors(graph, u)
                if v < w
                    add_edge!(graph, v, w)
                end
            end
        end

        u = last(supernode(stree, i))
        v = first(supernode(stree, parentindex(stree.tree, i)))
        
        for w in outneighbors(graph, u)
            if v < w
                add_edge!(graph, v, w)
            end
        end
    end

    graph
end


# Compute the width of a supernodal elimination tree.
function width(stree::SupernodeTree)
    maximum(stree.degree[stree.representative])
end


# Get the (sorted) supernode at node i.
function supernode(stree::SupernodeTree, i::Integer)
    v = stree.representative[i]
    n = stree.cardinality[i]
    v:v + n - 1
end


# Compute the (unsorted) seperators of every node in T.
function seperators(stree::SupernodeTree)
    n = length(stree.tree)
    seperator = Vector{Vector{Int}}(undef, n)
    graph = eliminationgraph(stree)
    
    for i in 1:n
        clique = collect(outneighbors(graph, stree.representative[i]))
        filter!(j -> stree.ancestor[i] <= j, clique)
        seperator[i] = clique
    end
    
    seperator
end
