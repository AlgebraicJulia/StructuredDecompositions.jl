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
# ----------------------------------------
#    graph    simple connected graph
#    ealg     elimination algorithm
#    stype    supernode type
# ----------------------------------------
function SupernodeTree(
        graph::AbstractSymmetricGraph,
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
    seperator = Vector{Set{Int}}(undef, n)

    for i in 1:n - 1
        seperator[i] = Set(stree.ancestor[i])

        for v in outneighbors(stree.graph, stree.representative[i])
            if stree.ancestor[i] < v
                push!(seperator[i], v)
            end
        end
    end

    for i in 1:n - 2
        j = parentindex(stree.tree, i)

        for v in seperator[i]
            if stree.ancestor[j] < v
                push!(seperator[j], v)
            end
        end
        
    end

    seperator[n] = Set()
    seperator
end


function Base.show(io::IO, stree::SupernodeTree)
    n = width(stree)
    print(io, "width: $n\nsupernodal elimination tree:\n")
    
    print_tree(io, IndexNode(stree.tree)) do io, node
        show(IOContext(io, :compact => true, :limit => true), permutation(stree.graph, supernode(stree, node.index)))
    end
end
