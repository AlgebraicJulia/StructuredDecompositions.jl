# An ordered graph (G, σ) equipped with a supernodal elimination tree T.
struct SupernodeTree
    tree::PostorderTree         # supernodal elimination tree
    graph::OrderedGraph         # ordered graph
    partition::Vector{Int}      # supernode partition
    representative::Vector{Int} # vector of representative vertices
    cardinality::Vector{Int}    # vector of supernode cardinalities
    ancestor::Vector{Int}       # vector of first ancestors
    degree::Vector{Int}         # vector of higher degrees
end


# Construct a supernodal elimination tree using an elimination algorithm.
# ----------------------------------------
#    graph    simple connected graph
#    ealg     elimination algorithm
#    stype    supernode type
# ----------------------------------------
function SupernodeTree(
    graph,
    ealg::Union{Order, EliminationAlgorithm}=DEFAULT_ELIMINATION_ALGORITHM,
    stype::SupernodeType=DEFAULT_SUPERNODE_TYPE)

    SupernodeTree(EliminationTree(graph, ealg), stype)
end


# Construct a supernodal elimination tree.
# ----------------------------------------
#    etree    elimination tree
#    stype    supernode type
# ----------------------------------------
function SupernodeTree(etree::EliminationTree, stype::SupernodeType=DEFAULT_SUPERNODE_TYPE)
    degree = outdegrees(etree)
    partition, supernode, parent, ancestor = stree(etree, degree, stype)
    tree = Tree(parent)

    order = postorder(tree)
    tree = PostorderTree(tree, order)
    partition = inverse(order, partition)
    supernode = view(supernode, order)
    ancestor = view(ancestor, order)

    order = Order(vcat(supernode...))
    graph = OrderedGraph(etree.graph, order)
    degree = view(degree, order)
    partition = view(partition, order)
    ancestor = inverse(order, ancestor)
    representative = inverse(order, map(first, supernode))
    cardinality = map(length, supernode)

    SupernodeTree(tree, graph, partition, representative, cardinality, ancestor, degree)
end


# Compute the width of a supernodal elimination tree.
function treewidth(stree::SupernodeTree)
    maximum(stree.degree[stree.representative])
end


# Get the (sorted) supernode at node i.
function supernode(stree::SupernodeTree, i::Integer)
    v = stree.representative[i]
    n = stree.cardinality[i]
    v:v + n - 1
end


# Get the unique node j satisfying i ∈ supernode(j).
function find_node(stree::SupernodeTree, i::Integer)
    stree.partition[i]
end


# Compute the (unsorted) seperators of every node in T.
function seperators(stree::SupernodeTree)
    n = treesize(stree.tree)
    seperator = Vector{Set{Int}}(undef, n)

    for i in 1:n - 1
        seperator[i] = Set(stree.ancestor[i])

        for v in outneighbors(stree.graph, stree.representative[i])
            if stree.ancestor[i] < v
                push!(seperator[i], v)
            end
        end
    end

    for j in 1:n - 1
        for i in childindices(stree.tree, j)
            for v in seperator[i]
                if stree.ancestor[j] < v
                    push!(seperator[j], v)
                end
            end
        end
    end

    seperator[n] = Set()
    seperator
end
