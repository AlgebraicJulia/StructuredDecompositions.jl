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
    OrderedGraph(jtree.stree.graph)
end


"""
    clique(jtree::JunctionTree, i::Integer)

Get the clique at node ``i``.
"""
function clique(jtree::JunctionTree, i::Integer)
    [residual(jtree, i); seperator(jtree, i)]
end


"""
    seperator(jtree::JunctionTree, i::Integer)

Get the seperator at node ``i``.
"""
function seperator(jtree::JunctionTree, i::Integer)
    permutation(jtree.stree.graph, jtree.seperator[i])
end


"""
    residual(jtree::JunctionTree, i::Integer)

Get the residual at node ``i``.
"""
function residual(jtree::JunctionTree, i::Integer)
    permutation(jtree.stree.graph, supernode(jtree.stree, i))
end


"""
    in_clique(jtree::JunctionTree, v::Integer)

Find a node `i` safisfying `v ∈ clique(jtree, i)`.
"""
function in_clique(jtree::JunctionTree, v::Integer)
    in_supernode(jtree.stree, inverse(jtree.stree.graph, v))
end


"""
    in_clique(jtree::JunctionTree, set::AbstractVector)

Find a node `i` satisfying `set ⊆ clique(jtree, i)`. In none exists, return `nothing`.
"""
function in_clique(jtree::JunctionTree, set::AbstractVector)
    minimum = length(jtree)
    maximum = 1

    for i in inverse(jtree.stree.graph, set)
        minimum = min(minimum, i)
        maximum = max(maximum, i)
    end

    if maximum <= last(jtree.seperator[minimum])
        i
    end
end


"""
    seperator_to_clique(jtree::JunctionTree, i::Integer)

Construct a vector `index` satisfying `seperator(jtree, i) == clique(jtree, i)[index]`
"""
function seperator_to_clique(jtree::JunctionTree, i::Integer)
    residual = supernode(jtree.stree, i)
    seperator = jtree.seperator[i]
    length(residual) + 1:length(residual) + length(seperator)
end


"""
    seperator_to_parent(jtree::JunctionTree, i::Integer)

Construct a vector `index` satisfying `seperator(jtree, i) == clique(jtree, parent(jtree, i))[index]`.
"""
function seperator_to_parent(jtree::JunctionTree, i::Integer)
    j = parentindex(jtree, i)
    residual = supernode(jtree.stree, j)
    seperator = jtree.seperator[j]

    map(jtree.seperator[i]) do v
        if v in residual
            v - first(residual) + 1
        else
            length(residual) + searchsortedfirst(seperator, v)
        end
    end
end


"""
    set_to_clique(jtree::JunctionTree, i::Integer, set::AbstractVector)

Construct a vector `index` satisfying `set == clique(jtree, i)[index]`.
"""
function set_to_clique(jtree::JunctionTree, i::Integer, set::AbstractVector)
    res = supernode(jtree.stree, i)
    sep = jtree.seperator[i]

    map(inverse(jtree.stree.graph, set)) do v
        if v in res
            v - first(res) + 1
        else
            length(res) + searchsortedfirst(sep, v)
        end
    end
end


"""
    width(jtree::JunctionTree)

Compute the width of a junction tree.
"""
function width(jtree::JunctionTree)
    width(jtree.stree)
end


# Multiline printing.
function Base.show(io::IO, ::MIME"text/plain", jtree::JunctionTree)
    n = width(jtree)
    print(io, "width: $n\njunction tree:\n")
    
    print_tree(io, IndexNode(jtree)) do io, node
        show(IOContext(io, :compact => true, :limit => true), clique(jtree, node.index))
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
