# A rooted tree.
abstract type AbstractTree end


# Compute a postordering of tree's vertices.
function postorder(tree::AbstractTree)
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


# Get the height of a tree.
function height(tree::AbstractTree)
    n = length(tree)
    maximum(map(i -> level(tree, i), 1:n))
end
