# An ordered graph (G, σ) equipped with a supernodal elimination tree T.
struct SupernodeTree <: AbstractTree
    tree::PostorderTree         # tree
    ograph::OrderedGraph        # ordered graph
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
    ograph = OrderedGraph(etree.ograph, order)

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
    
    SupernodeTree(tree, ograph, representative, cardinality, _ancestor, _degree)
end

    
# Chordal Graphs and Semidefinite Optimization
# Vanderberghe and Andersen
# Algorithm 4.1: Maximal supernodes and supernodal elimination tree.
function stree(etree::EliminationTree, degree::AbstractVector, stype::SupernodeType)
    n = length(etree)
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

        for w in childindices(etree, v)
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
    ograph = deepcopy(stree.ograph)
    n = length(stree)

    for i in 1:n - 1
        for u in supernode(stree, i)[1:end - 1]
            v = u + 1

            for w in outneighbors(ograph, u)
                if v < w
                    add_edge!(ograph, v, w)
                end
            end
        end

        u = last(supernode(stree, i))
        v = first(supernode(stree, parentindex(stree, i)))
        
        for w in outneighbors(ograph, u)
            if v < w
                add_edge!(ograph, v, w)
            end
        end
    end

    ograph
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
    n = length(stree)
    seperator = Vector{Vector{Int}}(undef, n)
    ograph = eliminationgraph(stree)
    
    for i in 1:n
        clique = collect(outneighbors(ograph, stree.representative[i]))
        filter!(j -> stree.ancestor[i] <= j, clique)
        seperator[i] = clique
    end
    
    seperator
end


# Get the number of nodes in T.
function Base.length(stree::SupernodeTree)
    length(stree.tree)
end


# Get the level of node i.
function level(stree::SupernodeTree, i)
    level(stree.tree, i)
end


# Get the vertex σ(i).
function order(stree::SupernodeTree, i)
    order(stree.ograph, i)
end


# Get the index σ⁻¹(v),
function inverse(stree::SupernodeTree, i)
    inverse(stree.ograph, i)
end


##########################
# Indexed Tree Interface #
##########################


function AbstractTrees.rootindex(stree::SupernodeTree)
    rootindex(stree.tree)
end


function AbstractTrees.parentindex(stree::SupernodeTree, i)
    parentindex(stree.tree, i)
end


function AbstractTrees.childindices(stree::SupernodeTree, i)
    childindices(stree.tree, i)
end


function AbstractTrees.NodeType(::Type{IndexNode{SupernodeTree, Int}})
    HasNodeType()
end


function AbstractTrees.nodetype(::Type{IndexNode{SupernodeTree, Int}})
    IndexNode{SupernodeTree, Int}
end
