# A left-child right-sibling binary tree.
# This type implements the indexed tree interface.
struct Tree
    parent::Vector{Int} # vector of parents
    head::Vector{Int}   # vector of first children
    next::Vector{Int}   # vector of next siblings
    root::Int           # root
end


# Construct a tree from a list of parents.
# ----------------------------------------
#    parent    list of parents
# ----------------------------------------
function Tree(parent::AbstractVector)
    n = length(parent)
    head = zeros(Int, n)
    next = zeros(Int, n)
    root = n

    for i in n:-1:1
        if parent[i] == i
            root = i
        else
            next[i] = head[parent[i]]
            head[parent[i]] = i
        end
    end

    Tree(parent, head, next, root)
end


#=
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
=#


# Compute a postordering of tree's vertices.
function postorder(tree::Tree)
    #### Stack ####

    top = Ref{Int}(0)
    items = Vector{Int}(undef, treesize(tree))

    function empty()
        iszero(top[])
    end

    function push(i)
        top[] += 1
        items[top[]] = i
    end

    function pop()
        top[] -= 1
        items[top[] + 1]
    end

    ###############

    push(rootindex(tree))
    order = Vector{Int}(undef, treesize(tree))
    index = Vector{Int}(undef, treesize(tree))
    head = copy(tree.head)
    count = 1

    while !empty()
        j = pop()
        i = head[j]

        if iszero(i)
            order[count] = j
            index[j] = count
            count += 1
        else
            head[j] = tree.next[i]
            push(j)
            push(i)
        end
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
    j = tree.parent[i]
    i != j ? j : nothing
end


function AbstractTrees.nextsiblingindex(tree::Tree, i::Integer)
    j = tree.next[i]
    !iszero(j) ? j : nothing
end


function AbstractTrees.childindices(tree::Tree, i::Integer)
    index = Int[]
    j = tree.head[i]

    while !iszero(j)
        push!(index, j)
        j = tree.next[j]
    end

    index
end


function AbstractTrees.SiblingLinks(::Type{IndexNode{Tree, Int}})
    StoredSiblings()
end


function AbstractTrees.NodeType(::Type{IndexNode{Tree, Int}})
    HasNodeType()
end


function AbstractTrees.nodetype(::Type{IndexNode{Tree, Int}})
    IndexNode{Tree, Int}
end
