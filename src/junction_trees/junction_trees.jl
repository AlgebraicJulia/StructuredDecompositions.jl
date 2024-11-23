"""
    JunctionTree

A [tree decomposition](https://en.wikipedia.org/wiki/Tree_decomposition) of a graph ``G``.
This type implements the [indexed tree interface](https://juliacollections.github.io/AbstractTrees.jl/stable/#The-Indexed-Tree-Interface).
"""
struct JunctionTree
    stree::SupernodeTree           # supernodal elimination tree
    seperator::Vector{Vector{Int}} # vector of seperators
end


"""
    JunctionTree(graph[, ealg::Union{Order, EliminationAlgorithm}[, stype::SupernodeType]])

Construct a tree decomposition of a connected simple graph, optionally specifying an elimination algorithm and
a supernode type.
"""
function JunctionTree(
    graph,
    ealg::Union{Order, EliminationAlgorithm}=DEFAULT_ELIMINATION_ALGORITHM,
    stype::SupernodeType=DEFAULT_SUPERNODE_TYPE)

    JunctionTree(SupernodeTree(graph, ealg, stype))
end


# Construct a junction tree.
# ----------------------------------------
#    stree    supernodal elimination tree
# ----------------------------------------
function JunctionTree(stree::SupernodeTree)
    JunctionTree(stree, map(sort ∘ collect, seperators(stree)))
end


"""
    Order(jtree::JunctionTree)

Construct a perfect elimination ordering.
"""
function Order(jtree::JunctionTree)
    Order(jtree.stree.graph) 
end


"""
    OrderedGraph(jtree::JunctionTree)

Construct the ordered graph ``(G, \\sigma)``, where ``\\sigma`` is a perfect elimination ordering.
"""
function OrderedGraph(jtree::JunctionTree)
    copy(jtree.stree.graph)
end


"""
    clique(jtree::JunctionTree, i::Integer)

Get the clique at node ``i``.
"""
function clique(jtree::JunctionTree, i::Integer)
    view(Order(jtree), [supernode(jtree.stree, i); jtree.seperator[i]])
end


"""
    seperator(jtree::JunctionTree, i::Integer)

Get the seperator at node ``i``.
"""
function seperator(jtree::JunctionTree, i::Integer)
    view(Order(jtree), jtree.seperator[i])
end


"""
    residual(jtree::JunctionTree, i::Integer)

Get the residual at node ``i``.
"""
function residual(jtree::JunctionTree, i::Integer)
    view(Order(jtree), supernode(jtree.stree, i))
end


"""
    find_clique(jtree::JunctionTree, v::Integer)

Find a node `i` safisfying `v ∈ clique(jtree, i)`.
"""
function find_clique(jtree::JunctionTree, v::Integer)
    find_supernode(jtree.stree, inv(Order(jtree))[v])
end


"""
    find_clique(jtree::JunctionTree, set::AbstractVector)

Find a node `i` satisfying `vertices ⊆ clique(jtree, i)`.
"""
function find_clique(jtree::JunctionTree, vertices::AbstractVector)
    find_supernode(jtree.stree, minimum(view(inv(Order(jtree)), vertices)))
end


"""
    treewidth(jtree::JunctionTree)

Compute the width of a junction tree.
"""
function treewidth(jtree::JunctionTree)
    treewidth(jtree.stree)
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
    residual = supernode(jtree.stree, i)
    seperator = jtree.seperator[i]
    length(residual) + 1:length(residual) + length(seperator)
end

# Compute the lift L: seperator(i) → clique(parent(i)). This satisfies
#    seperator(jtree, i) == clique(jtree, parentindex(jtree, i)[lift_sep_par(jtree, i)]
function lift_par(jtree::JunctionTree, i::Integer)
    lift_ind(jtree, jtree.seperator[i], parentindex(jtree, i))
end


# Compute the lift L: vertices → clique(i). This satisfies
#    vertices == clique(jtree, i)[lift(jtree, vertices, i)]
function lift(jtree::JunctionTree, vertices::AbstractVector, i::Integer)
    lift_ind(jtree, view(inv(Order(jtree)), vertices), i)
end


function lift_ind(jtree::JunctionTree, indices::AbstractVector, i::Integer)
    residual = supernode(jtree.stree, i)
    seperator = jtree.seperator[i]

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
    treesize(jtree.stree.tree)
end


function AbstractTrees.treeheight(jtree::JunctionTree)
    treeheight(jtree.stree.tree)
end


function AbstractTrees.rootindex(jtree::JunctionTree)
    rootindex(jtree.stree.tree)
end


function AbstractTrees.parentindex(jtree::JunctionTree, i::Integer)
    parentindex(jtree.stree.tree, i)
end


function AbstractTrees.childindices(jtree::JunctionTree, i::Integer)
    childindices(jtree.stree.tree, i)
end


function AbstractTrees.NodeType(::Type{IndexNode{JunctionTree, Int}})
    HasNodeType()
end


function AbstractTrees.nodetype(::Type{IndexNode{JunctionTree, Int}})
    IndexNode{JunctionTree, Int}
end
