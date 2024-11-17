# A rooted tree.
struct Tree
    root::Int                     # root
    parent::Vector{Int}           # parent
    children::Vector{Vector{Int}} # children
end


# Construct a tree from a list of parents.
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


# Get the number of nodes in a tree.
function Base.length(tree::Tree)
    length(tree.parent)
end


# Compute a postordering of tree's vertices.
function postorder(tree::Tree)
    n = length(tree)
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
