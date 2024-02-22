# A junction tree.
struct JunctionTree
    stree::EliminationTree
    seperatorlist::Vector{Vector{Int}}
end


# Construct a tree decomposition.
function JunctionTree(graph::AbstractSymmetricGraph, stree::EliminationTree)
    graph = makeeliminationgraph(graph, stree) 

    n = length(stree)
    seperatorlist = Vector{Vector{Int}}(undef, n) 
    seperatorlist[n] = []

    for i in 1:n - 1
        v₁ = stree.firstsupernodelist[i]
        v₂ = stree.lastsupernodelist[i]
        bag = collect(neighbors(graph, v₁))
        sort!(bag)
        seperatorlist[i] = bag[v₂ - v₁ + 1:end]
    end

    JunctionTree(stree, seperatorlist)
end


# Reorient a juncton tree towards the given root.
function JunctionTree(root::Integer, jtree::JunctionTree)
    m = length(jtree.stree.order)
    n = length(jtree)
    seperatorlist = Vector{Vector{Int}}(undef, n)
    supernodelist = Vector{Vector{Int}}(undef, n)
    subtreelist = Vector{Int}(undef, m)

    v₁ = jtree.stree.firstsupernodelist[root]
    v₂ = jtree.stree.lastsupernodelist[root]
    seperatorlist[n] = []
    supernodelist[n] = [v₁:v₂; jtree.seperatorlist[root]]
    subtreelist[supernodelist[n]] .= n 

    tree = Tree(root, jtree.stree.tree)
    postorder, tree = makepostorder(tree)

    for i in 1:n - 1
        j = postorder[i]
        v₁ = jtree.stree.firstsupernodelist[j]
        v₂ = jtree.stree.lastsupernodelist[j]

        if isdescendant(jtree, root, j)
            seperatorlist[i] = jtree.seperatorlist[postorder[parentindex(tree, i)]]
            supernodelist[i] = [v₁:v₂; jtree.seperatorlist[j]] 
            deletesorted!(supernodelist[i], seperatorlist[i])
        else
            seperatorlist[i] = jtree.seperatorlist[j]
            supernodelist[i] = v₁:v₂
        end

        subtreelist[supernodelist[i]] .= i
    end 

    order = jtree.stree.order
    width = jtree.stree.width
    stree = EliminationTree(order, tree, supernodelist, subtreelist, width)

    for i in 1:n
        seperatorlist[i] = stree.order.index[order[seperatorlist[i]]]
        sort!(seperatorlist[i])
    end

    JunctionTree(stree, seperatorlist)
end


# Construct a tree decomposition, first computing an elimination order and a supernodal
# elimination tree.
function JunctionTree(
    graph::AbstractSymmetricGraph,
    algorithm::Union{Order, EliminationAlgorithm}=DEFAULT_ELIMINATION_ALGORITHM,
    supernode::Supernode=DEFAULT_SUPERNODE)

    stree = EliminationTree(graph, algorithm, supernode)
    JunctionTree(graph, stree)
end


# Get the number of nodes in a junction tree.
function Base.length(jtree::JunctionTree)
    length(jtree.stree)
end


# Get the width of a junction tree.
function getwidth(jtree::JunctionTree)
    getwidth(jtree.stree)
end


# Get the seperator at node i.
function getseperator(jtree::JunctionTree, i::Integer)
    jtree.stree.order[jtree.seperatorlist[i]]
end


# Get the residual at node i.
function getresidual(jtree::JunctionTree, i::Integer)
    getsupernode(jtree.stree, i)
end


# Get the highest node containing the vertex v.
function getsubtree(jtree::JunctionTree, v::Union{Integer, AbstractVector})
    getsubtree(jtree.stree, v)
end


# Get the level of node i.
function getlevel(jtree::JunctionTree, i::Integer)
    getlevel(jtree.stree, i)
end


# Evaluate whether node i₁ is a descendant of node i₂.
function AbstractTrees.isdescendant(jtree::JunctionTree, i₁::Integer, i₂::Integer)
    isdescendant(jtree.stree, i₁, i₂)
end


# Construct an elimination graph.
function makeeliminationgraph(graph::AbstractSymmetricGraph, stree::EliminationTree)
    n = length(stree)
    graph = Graph(graph, stree.order)

    for i in 1:n - 1
        u₁ = stree.firstsupernodelist[i]
        u₂ = stree.lastsupernodelist[i]

        for u in u₁:u₂ - 1
            v = u + 1

            for w in neighbors(graph, u)
                if v != w && !has_edge(graph, v, w)
                    add_edge!(graph, v, w)
                end
            end
        end

        u = u₂
        v = stree.firstsupernodelist[parentindex(stree, i)]

        for w in neighbors(graph, u)
            if v != w && !has_edge(graph, v, w)
                add_edge!(graph, v, w)
            end
        end
    end

    graph
end


##########################
# Indexed Tree Interface #
##########################


function AbstractTrees.rootindex(jtree::JunctionTree)
    rootindex(jtree.stree)
end


function AbstractTrees.parentindex(jtree::JunctionTree, i::Integer)
    parentindex(jtree.stree, i)
end


function AbstractTrees.childindices(jtree::JunctionTree, i::Integer)
    childindices(jtree.stree, i)
end


function AbstractTrees.NodeType(::Type{IndexNode{JunctionTree, Int}})
    HasNodeType()
end


function AbstractTrees.nodetype(::Type{IndexNode{JunctionTree, Int}})
    IndexNode{JunctionTree, Int}
end
