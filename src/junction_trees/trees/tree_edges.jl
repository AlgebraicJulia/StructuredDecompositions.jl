# The edges of a rooted tree.
struct TreeEdges{T <: AbstractTree} <: AbstractVector{Graphs.SimpleEdge{Int}}
    tree::T
end


function Base.show(io::IO, edges::TreeEdges)
    n = length(edges)
    print(io, "TreeEdges $n")
end 


######################
# Iterator Interface #
######################


function Base.hasfastin(::Type{TreeEdges{T}}) where T <: AbstractTree
    true
end


function Base.in(edge, edges::TreeEdges)
    has_edge(edges.tree, edge)
end


#############################
# Abstract Vector Interface #
#############################


function Base.getindex(edges::TreeEdges, i::Integer)
    if rootindex(edges.tree) <= i
        i += 1
    end

    Graphs.SimpleEdge{Int}(parentindex(edges.tree, i), i)
end


function Base.IndexStyle(::Type{TreeEdges{T}}) where T <: AbstractTree
    IndexLinear()
end


function Base.size(edges::TreeEdges)
    (ne(edges.tree),)
end


############################
# Abstract Graph Interface #
############################


function Graphs.edges(tree::AbstractTree)
    TreeEdges(tree)
end
