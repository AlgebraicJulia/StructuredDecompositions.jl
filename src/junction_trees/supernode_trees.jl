# An ordered graph (G, σ) equipped with a supernodal elimination tree T.
struct SupernodeTree
    tree::PostorderTree         # supernodal elimination tree
    graph::OrderedGraph        # ordered graph
    representative::Vector{Int} # representative vertex
    cardinality::Vector{Int}    # supernode cardinality
    ancestor::Vector{Int}       # first ancestor
    degrees::Vector{Int}        # higher degrees
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
    sorder = postorder(tree)
    tree = PostorderTree(tree, sorder)
    
    order = Order(vcat(snode[sorder]...))
    graph = OrderedGraph(etree.graph, order)

    n = length(tree)
    representative = zeros(Int, n)
    cardinality = zeros(Int, n)
    _ancestor = zeros(Int, n)
    _degree = zeros(Int, n)
    
    for i in 1:n - 1
        j = sorder[i]
        representative[i] = inverse(order, snode[j][1])
        cardinality[i] = length(snode[j])
        _degree[i] = degree[j]
        _ancestor[i] = inverse(order, ancestor[j])
    end

    representative[n] = inverse(order, snode[n][1])
    cardinality[n] = length(snode[n])
    _degree[n] = degree[n]
    _ancestor[n] = ancestor[n]
    
    SupernodeTree(tree, graph, representative, cardinality, _ancestor, _degree)
end

    
# Chordal Graphs and Semidefinite Optimization
# Vanderberghe and Andersen
# Algorithm 4.1: Maximal supernodes and supernodal elimination tree.
function stree(etree::EliminationTree, degree::AbstractVector, stype::SupernodeType)
    n = length(etree.tree)
    index = zeros(Int, n)
    snd = Vector{Int}[]
    q = Int[]
    a = Int[]

    for v in 1:n
        ww = findchild(etree, degree, stype, v)
        
        if isnothing(ww)
            i = length(snd) + 1
            index[v] = i
            push!(snd, [v])
            push!(q, length(snd))
            push!(a, n + 1)
        else
            i = index[ww]
            index[v] = i
            push!(snd[i], v)
        end

        for w in childindices(etree.tree, v)
            if w !== ww
                j = index[w]
                q[j] = i
                a[j] = v
            end
        end
    end

    snd, q, a
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
