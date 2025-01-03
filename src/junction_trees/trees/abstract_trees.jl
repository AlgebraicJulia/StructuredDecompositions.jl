# A rooted tree.
# This type implements the abstract graph and abstract tree interfaces.
abstract type AbstractTree <: Graphs.AbstractSimpleGraph{Int} end


# Compute a postordering of a rooted tree.
function dfs(tree::AbstractTree)
    head = Vector{Int}(undef, nv(tree))
    order = FixedStack{Int}(nv(tree))
    stack = FixedStack{Int}(nv(tree))
    push!(stack, tree.root[])

    for i in vertices(tree)
        head[i] = ntz(firstchildindex(tree, i))
    end

    while !isempty(stack)
        j = last(stack)
        i = head[j]

        if iszero(i)
            push!(order, pop!(stack))
        else
            head[j] = ntz(nextsiblingindex(tree, i))
            push!(stack, i)
        end
    end

    order.items
end


# Get the level of every vertex in a topologically ordered rooted tree.
function levels(tree::AbstractTree)
    level = Vector{Int}(undef, nv(tree))
    level[end] = 0

    for edge in Iterators.reverse(Graphs.edges(tree))
        i = Graphs.dst(edge)
        j = Graphs.src(edge)
        level[i] = level[j] + 1
    end

    level
end


# Get the first descendant of every vertex in a postordered rooted tree.
function firstdescendants(tree::AbstractTree)
    fdesc = collect(vertices(tree))

    for edge in Graphs.edges(tree)
        i = Graphs.dst(edge)
        j = Graphs.src(edge)
        fdesc[j] = min(fdesc[i], fdesc[j])
    end

    fdesc
end


# Find a child w of v such that v ∈ supernode(w).
# If no such child exists, return nothing.
function child_in_supernode(tree::AbstractTree, colcount::AbstractVector, stype::Node, v::Integer) end


# Find a child w of v such that v ∈ supernode(w).
# If no such child exists, return nothing.
function child_in_supernode(tree::AbstractTree, colcount::AbstractVector, stype::Maximal, v::Integer)
    for u in childindices(tree, v)
        if colcount[u] == colcount[v] + 1
            return u
        end
    end
end


# Find a child w of v such that v ∈ supernode(w).
# If no such child exists, return nothing.
function child_in_supernode(tree::AbstractTree, colcount::AbstractVector, stype::Fundamental, v::Integer)
    u = firstchildindex(tree, v)

    if !isnothing(u) && colcount[u] == colcount[v] + 1 && isnothing(nextsiblingindex(tree, u))
        return u
    end
end


function Base.show(io::IO, ::MIME"text/plain", tree::AbstractTree)
    print_tree(io, IndexNode(tree))
end


###########################
# Abstract Tree Interface #
###########################


function AbstractTrees.ParentLinks(::Type{IndexNode{T, Int}}) where T <: AbstractTree
    StoredParents()
end


function AbstractTrees.SiblingLinks(::Type{IndexNode{T, Int}}) where T <: AbstractTree
    StoredSiblings()
end


function AbstractTrees.NodeType(::Type{IndexNode{T, Int}}) where T <: AbstractTree
    HasNodeType()
end


function AbstractTrees.nodetype(::Type{IndexNode{T, Int}}) where T <: AbstractTree
    IndexNode{T, Int}
end


############################
# Abstract Graph Interface #
############################


function Graphs.ne(tree::AbstractTree)
    nv(tree) - 1
end


function Graphs.is_directed(::Type{T}) where T <: AbstractTree
    true
end


function Graphs.outneighbors(tree::AbstractTree, i::Integer)
    childindices(tree, i)
end


function Graphs.edgetype(tree::AbstractTree)
    Graphs.SimpleEdge{Int}
end


function Graphs.has_edge(tree::AbstractTree, i::Integer, j::Integer)
    i == parentindex(tree, j)
end
