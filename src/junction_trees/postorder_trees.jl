# A postordered rooted tree.
# This type implements the indexed tree interface.
struct PostorderTree
    parent::Vector{Int}           # vector of parents
    children::Vector{Vector{Int}} # vector of children
    level::Vector{Int}            # vector of levels
    descendant::Vector{Int}       # vector of first descendants
end


# Construct a tree from a postordered list of parents.
# ----------------------------------------
#    parent    list of parents
# ----------------------------------------
function PostorderTree(parent::AbstractVector)
    n = length(parent)
    children = Vector{Vector{Int}}(undef, n)
    level = Vector{Int}(undef, n)
    descendant = Vector{Int}(undef, n)
    
    for i in 1:n
        children[i] = []
        level[i] = 0
        descendant[i] = i
    end
    
    for i in 1:n - 1
        j = parent[i]
        push!(children[j], i)
        descendant[j] = min(descendant[i], descendant[j])
    end
    
    for i in n - 1:-1:1
        j = parent[i]
        level[i] = level[j] + 1
    end
    
    PostorderTree(parent, children, level, descendant)
end


# Postorder a tree.
# ----------------------------------------
#    tree    tree
#    order   postorder
# ----------------------------------------
function PostorderTree(tree::Tree, order::Order)
    n = treesize(tree)
    parent = Vector{Int}(undef, n)
    
    for i in 1:n - 1
        parent[i] = inverse(order, parentindex(tree, order[i]))
    end
    
    parent[n] = n
    PostorderTree(parent)
end


# Get the level of a node i.
function level(tree::PostorderTree, i::Integer)
    tree.level[i]
end


# Get the first descendant of a node i.
function firstdescendant(tree::PostorderTree, i::Integer)
    tree.descendant[i]
end


##########################
# Indexed Tree Interface #
##########################


function AbstractTrees.treesize(tree::PostorderTree)
    length(tree.parent)
end


function AbstractTrees.treeheight(tree::PostorderTree)
    maximum(tree.level)
end


function AbstractTrees.parentindex(tree::PostorderTree, i::Integer)
    if i != rootindex(tree)
        tree.parent[i]
    end
end


function AbstractTrees.childindices(tree::PostorderTree, i::Integer)
    tree.children[i]
end


function AbstractTrees.rootindex(tree::PostorderTree)
    treesize(tree)
end


function AbstractTrees.NodeType(::Type{IndexNode{PostorderTree, Int}})
    HasNodeType()
end


function AbstractTrees.nodetype(::Type{IndexNode{PostorderTree, Int}})
    IndexNode{PostorderTree, Int}
end
