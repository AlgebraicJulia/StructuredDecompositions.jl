"""
    JunctionTree

A [tree decomposition](https://en.wikipedia.org/wiki/Tree_decomposition) of a graph ``G``.
This type implements the [indexed tree interface](https://juliacollections.github.io/AbstractTrees.jl/stable/#The-Indexed-Tree-Interface).
"""
struct JunctionTree
    order::Order                # elimination order
    tree::PostorderTree         # supernodal elimination tree
    representative::Vector{Int} # vector of representative vertices
    seperator::Vector{Vector{Int}} # vector of seperators

    # cache
    partition::Vector{Int}      # supernode partition
    degree::Vector{Int}         # vector of higher degrees
end


"""
    JunctionTree(graph[, ealg::Union{Order, EliminationAlgorithm}[, stype::SupernodeType]])

Construct a tree decomposition of a connected simple graph, optionally specifying an elimination algorithm and
a supernode type.
"""
function JunctionTree(
    graph,
    ealg::EliminationAlgorithm=DEFAULT_ELIMINATION_ALGORITHM,
    stype::SupernodeType=DEFAULT_SUPERNODE_TYPE)

    graph = adjacencymatrix(graph)
    order = Order(graph, ealg)
    JunctionTree(graph, order, stype)
end


function JunctionTree(
    graph,
    order::Order,
    stype::SupernodeType=DEFAULT_SUPERNODE_TYPE)

    graph = OrderedGraph(graph, order)
    tree = Tree(etree(graph))
    JunctionTree(order, graph, tree, stype)
end


# Construct a junction tree.
# ----------------------------------------
#    stree    supernodal elimination tree
# ----------------------------------------
function JunctionTree(order::Order, graph::OrderedGraph, tree::Tree, stype::SupernodeType=DEFAULT_SUPERNODE_TYPE)
    degree = outdegrees(graph, tree)
    partition, supernode, parent, ancestor = stree(tree, degree, stype)
    tree = Tree(parent)

    order = postorder(tree)
    tree = PostorderTree(tree, order)
    partition = view(inv(order), partition)
    supernode = view(supernode, order)
    ancestor = view(ancestor, order)

    order = Order(vcat(supernode...))
    graph = OrderedGraph(graph, order)
    degree = view(degree, order)
    partition = view(partition, order)
    ancestor = view(inv(order), ancestor)

    representative = Vector{Int}(undef, treesize(tree) + 1)
    representative[1:end - 1] .= view(inv(order), map(first, supernode))
    representative[end] = representative[end - 1] + length(supernode[end])

    seperator = map(sort ∘ collect, seperators(graph, tree, representative, ancestor))
    JunctionTree(order, tree, representative, seperator, partition, degree)
end


# Compute the (unsorted) seperators of every node in T.
function seperators(graph::OrderedGraph, tree::PostorderTree, representative::AbstractVector, ancestor::AbstractVector)
    n = treesize(tree)
    seperator = Vector{Set{Int}}(undef, n)

    for i in 1:n - 1
        seperator[i] = Set(ancestor[i])

        for v in outneighbors(graph, representative[i])
            if ancestor[i] < v
                push!(seperator[i], v)
            end
        end
    end

    for j in 1:n - 1
        for i in childindices(tree, j)
            for v in seperator[i]
                if ancestor[j] < v
                    push!(seperator[j], v)
                end
            end
        end
    end

    seperator[n] = Set()
    seperator
end


"""
    Order(jtree::JunctionTree)

Construct a perfect elimination ordering.
"""
function Order(jtree::JunctionTree)
    Order(jtree.order) 
end


"""
    clique(jtree::JunctionTree, i::Integer)

Get the clique at node ``i``.
"""
function clique(jtree::JunctionTree, i::Integer)
    view(jtree.order, [residualindices(jtree, i); seperatorindices(jtree, i)])
end


"""
    treewidth(jtree::JunctionTree)

Compute the width of a junction tree.
"""
function treewidth(jtree::JunctionTree)
    maximum(jtree.degree)
end


# Get the (sorted) supernode at node i.
function residualindices(jtree::JunctionTree, i::Integer)
    jtree.representative[i]:jtree.representative[i + 1] - 1
end


function seperatorindices(jtree::JunctionTree, i::Integer)
    jtree.seperator[i]
end


"""
    seperator(jtree::JunctionTree, i::Integer)

Get the seperator at node ``i``.
"""
function seperator(jtree::JunctionTree, i::Integer)
    view(jtree.order, seperatorindices(jtree, i))
end


"""
    residual(jtree::JunctionTree, i::Integer)

Get the residual at node ``i``.
"""
function residual(jtree::JunctionTree, i::Integer)
    view(jtree.order, residualindices(jtree, i))
end


"""
    find_clique(jtree::JunctionTree, v::Integer)

Find a node `i` safisfying `v ∈ clique(jtree, i)`.
"""
function find_clique(jtree::JunctionTree, v::Integer)
    jtree.partition[inv(jtree.order)[v]]
end


"""
    find_clique(jtree::JunctionTree, set::AbstractVector)

Find a node `i` satisfying `vertices ⊆ clique(jtree, i)`.
"""
function find_clique(jtree::JunctionTree, vertices::AbstractVector)
    jtree.partition[minimum(view(jtree.order, vertices))]
end


# Multiline printing.
function Base.show(io::IO, ::MIME"text/plain", jtree::JunctionTree)
    n = treewidth(jtree)
    print(io, "width: $n\njunction tree:\n")
    
    print_tree(io, IndexNode(jtree)) do io, node
        show(IOContext(io, :compact => true, :limit => true), clique(jtree, node.index))
    end
end


#########
# Lifts #
#########
# In principal, cliques are subsets C ⊆ V. In practice, we represent them by vectors
#    C: n → V
# Given another vector S: m → V, we may wish to find a vector L: m → n satisfying
#      C
#    n → V
#      ↖ ↑ S
#      L m
# The vector L is called a lift, and we write L: S → C.


# Compute the lift L: seperator(i) → clique(i). This satisfies
#    seperator(jtree, i) == clique(jtree, i)[lift_sep(jtree, i)]
function lift_sep(jtree::JunctionTree, i::Integer)
    residual = residualindices(jtree, i)
    seperator = seperatorindices(jtree, i)
    length(residual) + 1:length(residual) + length(seperator)
end

# Compute the lift L: seperator(i) → clique(parent(i)). This satisfies
#    seperator(jtree, i) == clique(jtree, parentindex(jtree, i)[lift_sep_par(jtree, i)]
function lift_par(jtree::JunctionTree, i::Integer)
    lift_ind(jtree, seperatorindices(jtree, i), parentindex(jtree, i))
end


# Compute the lift L: vertices → clique(i). This satisfies
#    vertices == clique(jtree, i)[lift(jtree, vertices, i)]
function lift(jtree::JunctionTree, vertices::AbstractVector, i::Integer)
    lift_ind(jtree, view(inv(jtree.order), vertices), i)
end


function lift_ind(jtree::JunctionTree, indices::AbstractVector, i::Integer)
    residual = residualindices(jtree, i)
    seperator = seperatorindices(jtree,i)

    map(indices) do v
        if v in residual
            v - first(residual) + 1
        else
            length(residual) + searchsortedfirst(seperator, v)
        end
    end
end


##########################
# Indexed Tree Interface #
##########################


function AbstractTrees.treesize(jtree::JunctionTree)
    treesize(jtree.tree)
end


function AbstractTrees.treeheight(jtree::JunctionTree)
    treeheight(jtree.tree)
end


function AbstractTrees.rootindex(jtree::JunctionTree)
    rootindex(jtree.tree)
end


function AbstractTrees.parentindex(jtree::JunctionTree, i::Integer)
    parentindex(jtree.tree, i)
end


function AbstractTrees.childindices(jtree::JunctionTree, i::Integer)
    childindices(jtree.tree, i)
end


function AbstractTrees.NodeType(::Type{IndexNode{JunctionTree, Int}})
    HasNodeType()
end


function AbstractTrees.nodetype(::Type{IndexNode{JunctionTree, Int}})
    IndexNode{JunctionTree, Int}
end
