# A rooted tree.
# This type implements the indexed tree interface.
struct Tree
    root::Int                     # root
    parent::Vector{Int}           # vector of parents
    children::Vector{Vector{Int}} # vector of children
end


# Construct a tree from a list of parents.
# ----------------------------------------
#    parent    list of parents
# ----------------------------------------
function Tree(parent::AbstractVector)
    n = root = length(parent)
    children = Vector{Vector{Int}}(undef, n)
    
    for i in 1:n
        children[i] = []
    end

    for i in 1:n
        j = parent[i]
        
        if i == j
            root = i
        else
            push!(children[j], i)
        end
    end

    Tree(root, parent, children)
end


# Compute a postordering of tree's vertices.
function postorder(tree::Tree)
    n = treesize(tree)
    order = Vector{Int}(undef, n)
    index = Vector{Int}(undef, n)
    
    for node in PreOrderDFS(IndexNode(tree))
        order[n] = node.index
        index[node.index] = n
        n -= 1
    end
    
    Order(order, index)
end


##########################
# Indexed Tree Interface #
##########################


function AbstractTrees.treesize(tree::Tree)
    length(tree.parent)
end


function AbstractTrees.rootindex(tree::Tree)
    tree.root
end


function AbstractTrees.parentindex(tree::Tree, i::Integer)
    if i != rootindex(tree)
        tree.parent[i]
    end
end


function AbstractTrees.childindices(tree::Tree, i::Integer)
    tree.children[i]
end


function AbstractTrees.NodeType(::Type{IndexNode{Tree, Int}})
    HasNodeType()
end


function AbstractTrees.nodetype(::Type{IndexNode{Tree, Int}})
    IndexNode{Tree, Int}
end
