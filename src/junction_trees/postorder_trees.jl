# A postordered rooted tree.
# This type implements the indexed tree interface.
struct PostorderTree
    tree::Tree         # rooted tree
    level::Vector{Int} # vector of levels
    fdesc::Vector{Int} # vector of first descendants
end


# Construct a tree from a postordered list of parents.
# ----------------------------------------
#    parent    list of parents
# ----------------------------------------
function PostorderTree(parent::AbstractVector)
    tree = Tree(parent)
    n = treesize(tree)
    level = Vector{Int}(undef, n)
    fdesc = Vector{Int}(undef, n)
    level[n] = 0
    fdesc[n] = n

    for i in n - 1:-1:1
        level[i] = level[parentindex(tree, i)] + 1
        fdesc[i] = i
    end

    for i in 1:n - 1
        fdesc[parentindex(tree, i)] = min(fdesc[i], fdesc[parentindex(tree, i)])
    end
    
    PostorderTree(tree, level, fdesc)
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
        parent[i] = inv(order)[parentindex(tree, order[i])]
    end
    
    parent[n] = 0
    PostorderTree(parent)
end


# Get the level of a node i.
function level(tree::PostorderTree, i::Integer)
    tree.level[i]
end


# Get the first descendant of a node i.
function fdesc(tree::PostorderTree, i::Integer)
    tree.fdesc[i]
end


##########################
# Indexed Tree Interface #
##########################


function AbstractTrees.treesize(tree::PostorderTree)
    treesize(tree.tree)
end


function AbstractTrees.treeheight(tree::PostorderTree)
    maximum(tree.level)
end


function AbstractTrees.parentindex(tree::PostorderTree, i::Integer)
    parentindex(tree.tree, i)
end


function AbstractTrees.childindices(tree::PostorderTree, i::Integer)
    childindices(tree.tree, i)
end


function AbstractTrees.rootindex(tree::PostorderTree)
    rootindex(tree.tree)
end


function AbstractTrees.NodeType(::Type{IndexNode{PostorderTree, Int}})
    HasNodeType()
end


function AbstractTrees.nodetype(::Type{IndexNode{PostorderTree, Int}})
    IndexNode{PostorderTree, Int}
end
