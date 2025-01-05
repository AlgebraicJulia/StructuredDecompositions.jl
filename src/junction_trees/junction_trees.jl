"""
    JunctionTree

A junction tree.
"""
struct JunctionTree
    stree::SupernodeTree           # supernodal elimination tree
    seperator::Vector{Vector{Int}} # seperator
end
export JunctionTree


"""
    JunctionTree(graph::AbstractSymmetricGraph[, ealg::Union{Order, EliminationAlgorithm}[, stype::SupernodeType]])

Construct a tree decomposition of a connected simple graph, optionally specifying an elimination algorithm and
a supernode type.
"""
function JunctionTree(
    graph::AbstractSymmetricGraph,
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
    clique(jtree::JunctionTree, i::Integer)

Get the clique at node i.
"""
function clique(jtree::JunctionTree, i::Integer)
    [residual(jtree, i); seperator(jtree, i)]
end
export clique

"""
    seperator(jtree::JunctionTree, i::Integer)

Get the seperator at node i.
"""
function seperator(jtree::JunctionTree, i::Integer)
    permutation(jtree.stree.graph, jtree.seperator[i])
end
export seperator

"""
    residual(jtree::JunctionTree, i::Integer)

Get the residual at node i.
"""
function residual(jtree::JunctionTree, i::Integer)
    permutation(jtree.stree.graph, supernode(jtree.stree, i))
end
export residual


# Construct the inclusion seperator(i) → clique(parent(i)).
function seperator_to_parent(jtree::JunctionTree, i::Integer)
    j = parentindex(jtree, i)
    sep = jtree.seperator[i]
    sep_parent = jtree.seperator[j]
    res_parent = supernode(jtree.stree, j)

    i = 0
    index = Vector{Int}(undef, length(sep))
    
    for (j, v) in enumerate(sep)
        if v in res_parent
            index[j] = v - first(res_parent) + 1
        else
            i += searchsortedfirst(view(sep_parent, i + 1:length(sep_parent)), v)
            index[j] = length(res_parent) + i
        end
    end

    index
end
export seperator_to_parent


# Construct the inclusion seperator(i) → clique(i).
function seperator_to_self(jtree::JunctionTree, i::Integer)
    sep = jtree.seperator[i]
    res = supernode(jtree.stree, i)
    length(res) + 1:length(res) + length(sep)
end
export seperator_to_self

"""
    length(jtree::JunctionTree)

Get the number of nodes in a junction tree.
"""
function Base.length(jtree::JunctionTree)
    length(jtree.stree.tree)
end


"""
    height(jtree::JunctionTree)

Compute the height of a junction tree.
"""
function height(jtree::JunctionTree)
    height(jtree.stree.tree)
end
export height

"""
    width(jtree::JunctionTree)

Compute the width of a junction tree.
"""
function width(jtree::JunctionTree)
    width(jtree.stree)
end
export width

function Base.show(io::IO, jtree::JunctionTree)
    n = width(jtree)
    print(io, "width: $n\njunction tree:\n")
    
    print_tree(io, IndexNode(jtree)) do io, node
        show(IOContext(io, :compact => true, :limit => true), clique(jtree, node.index))
    end
end


##########################
# Indexed Tree Interface #
##########################


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
