# A rooted tree.
struct Tree
    root::Int
    parentlist::Vector{Int}
    childrenlist::Vector{Vector{Int}}
    levellist::Vector{Int}
    firstdescendantlist::Vector{Int}
end


# Orient a tree towards the given root.
function Tree(root::Integer, tree::Tree)
    i = root
    parent = parentindex(tree, i)
    parentlist = copy(tree.parentlist)
    childrenlist = deepcopy(tree.childrenlist)

    while !isnothing(parent)
        parentlist[parent] = i 
        push!(childrenlist[i], parent)
        deletesorted!(childrenlist[parent], i)
        i = parent
        parent = parentindex(tree, i)
    end

    Tree(root, parentlist, childrenlist) 
end


# Construct a tree from a list of parent and a list of children.
function Tree(root::Integer, parentlist::AbstractVector, childrenlist::AbstractVector)
    n = length(parentlist)
    levellist = Vector{Int}(undef, n)
    firstdescendantlist = Vector{Int}(undef, n)
    Tree(root, parentlist, childrenlist, levellist, firstdescendantlist)    
end


# Construct a tree from a list of parents.
function Tree(root::Integer, parentlist::AbstractVector)
    n = length(parentlist)
    childrenlist = Vector{Vector{Int}}(undef, n)
    childrenlist .= [[]]

    for i in 1:n
        if i != root
            push!(childrenlist[parentlist[i]], i)
        end
    end

    Tree(root, parentlist, childrenlist)
end


# Construct an elimination tree.
function Tree(graph::AbstractSymmetricGraph, order::Order)
    n = nv(graph)
    parentlist = makeetree(graph, order)
    @assert count(parentlist .== 0) == 1
    Tree(n, parentlist)
end


function Base.length(tree::Tree)
    length(tree.parentlist)
end


# Compute the parent vector of the elimination tree of the elimination graph of a ordered
# graph.
#
# The complexity is
# 𝒪(m log(n))
# where m = |E| and n = |V|.
#
# doi:10.1145/6497.6499
# Algorithm 4.2: Elimination Tree by Path Compression
function makeetree(graph::AbstractSymmetricGraph, order::Order)
    graph = Graph(graph, order)

    n = nv(graph)
    parent = Vector{Int}(undef, n)
    ancestor = Vector{Int}(undef, n)

    for i in 1:n
        parent[i] = 0
        ancestor[i] = 0

        for k in inneighbors(graph, i)
            r = k

            while ancestor[r] != 0 && ancestor[r] != i
                t = ancestor[r]
                ancestor[r] = i
                r = t
            end

            if ancestor[r] == 0
                ancestor[r] = i
                parent[r] = i
            end
        end
    end

    parent
end


# Given an ordered graph
# (G, σ),
# construct a directed graph by ordering the edges in G from lower to higher index.
#
# The complexity is
# 𝒪(m)
# where m = |E|.
function BasicGraphs.Graph(graph::AbstractSymmetricGraph, order::Order)
    n = nv(graph)
    digraph = Graph(n)

    for v in vertices(graph)
        i = order.index[v]

        for w in neighbors(graph, v)
            j = order.index[w]

            if i < j
                add_edge!(digraph, i, j)
            end
        end
    end

    digraph
end


##############
# Postorders #
##############


# Get the level of node i.
# This function only works on postordered trees.
function getlevel(tree::Tree, i::Integer)
    tree.levellist[i]
end


# Get the first descendant of node i.
# This function only works on postordered trees.
function getfirstdescendant(tree::Tree, i::Integer)
    tree.firstdescendantlist[i]
end


# Evaluate whether node i₁ is a descendant of node i₂.
# This function only works on postordered trees.
function AbstractTrees.isdescendant(tree::Tree, i₁::Integer, i₂::Integer)
    getfirstdescendant(tree, i₂) <= i₁ < i₂
end


# Compute a postordering of a tree.
#
# The complexity is
# 𝒪(n)
# where n = |V|.
function makepostorder(tree::Tree)
    n = length(tree)
    order = Order(n)
    parentlist = Vector{Int}(undef, n)
    childrenlist = Vector{Vector{Int}}(undef, n)
    levellist = Vector{Int}(undef, n)
    firstdescendantlist = Vector{Int}(undef, n)

    root, nodes... = PreOrderDFS(IndexNode(tree))

    order[n] = root.index
    parentlist[n] = 0
    childrenlist[n] = []
    levellist[n] = 0

    for (i, node) in enumerate(nodes)
        j = n - i
        order[j] = node.index

        k = order.index[parentindex(tree, node.index)]
        parentlist[j] = k
        childrenlist[j] = []
        pushfirst!(childrenlist[k], j)
        levellist[j] = 1 + levellist[k]
    end 

    for i in 1:n
        init = i
        firstdescendantlist[i] = minimum(firstdescendantlist[childrenlist[i]]; init)
    end

    tree = Tree(n, parentlist, childrenlist, levellist, firstdescendantlist)
    order, tree
end


##########################
# Indexed Tree Interface #
##########################


function AbstractTrees.rootindex(tree::Tree)
    tree.root
end


function AbstractTrees.parentindex(tree::Tree, i::Integer)
    if i != rootindex(tree)
        tree.parentlist[i]
    end
end


function AbstractTrees.childindices(tree::Tree, i::Integer)
    tree.childrenlist[i]
end


function AbstractTrees.NodeType(::Type{IndexNode{Tree, Int}})
    HasNodeType()
end


function AbstractTrees.nodetype(::Type{IndexNode{Tree, Int}})
    IndexNode{Tree, Int}
end
