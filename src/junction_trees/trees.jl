# A rooted tree.
# This type implements the indexed tree interface.
struct Tree
    parent::Vector{Int}  # vector of parents
    child::Vector{Int}   # vector of left-children
    brother::Vector{Int} # vector of right-siblings
    root::Array{Int, 0}  # root
end


# Construct a tree from a list of parents.
# ----------------------------------------
#    parent    list of parents
# ----------------------------------------
function Tree(parent::AbstractVector)
    n = length(parent)
    child = zeros(Int, n)
    brother = zeros(Int, n)
    root = fill(n)

    for i in n:-1:1
        if iszero(parent[i])
            root .= i
        else
            brother[i] = child[parent[i]]
            child[parent[i]] = i
        end
    end

    Tree(parent, child, brother, root)
end


# Compute a postordering of tree's vertices.
function postorder(tree::Tree)
    n = treesize(tree)
    stack = FixedStack{Int}(n)
    push!(stack, rootindex(tree))
    order = Order(undef, n)
    child = copy(tree.child)
    count = 1

    while !isempty(stack)
        j = pop!(stack)
        i = child[j]

        if iszero(i)
            order[count] = j
            count += 1
        else
            child[j] = tree.brother[i]
            push!(stack, j)
            push!(stack, i)
        end
    end
   
    order
end


function postorder!(tree::Tree)
    order = postorder(tree)
    parent = copy(tree.parent)

    map!(tree.parent, order) do i
        j = parent[i]
        iszero(j) ? 0 : inv(order)[j]
    end

    n = treesize(tree)
    tree.child .= 0
    tree.root .= n

    for i in n - 1:-1:1
        tree.brother[i] = tree.child[tree.parent[i]]
        tree.child[tree.parent[i]] = i
    end

    order
end


function level(tree::Tree)
    n = treesize(tree)
    level = Vector{Int}(undef, n)
    level[n] = 0

    for i in n - 1:-1:1
        j = parentindex(tree, i)
        level[i] = level[j] + 1
    end

    level
end


function fdesc(tree::Tree)
    n = treesize(tree)
    fdesc = collect(1:n)

    for i in 1:n - 1
        j = parentindex(tree, i)
        fdesc[j] = min(fdesc[i], fdesc[j])
    end

    fdesc
end


##########################
# Indexed Tree Interface #
##########################


function AbstractTrees.treesize(tree::Tree)
    length(tree.parent)
end


function AbstractTrees.rootindex(tree::Tree)
    tree.root[]
end


function AbstractTrees.parentindex(tree::Tree, i::Integer)
    j = tree.parent[i]

    if !iszero(j)
        j
    end
end


function AbstractTrees.nextsiblingindex(tree::Tree, i::Integer)
    j = tree.brother[i]

    if !iszero(j)
        j
    end
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
