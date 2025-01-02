# A directed graph whose edges i → j satisfy i > j.
# This type implements the abstract graph interface.
struct OrderedGraph <: Graphs.AbstractSimpleGraph{Int}
    lower::SparseMatrixCSC{Bool, Int}
    upper::SparseMatrixCSC{Bool, Int}
end


function OrderedGraph()
    OrderedGraph(spzeros(Bool, 0, 0), spzeros(Bool, 0, 0))
end


function OrderedGraph(graph::SparseMatrixCSC)
    OrderedGraph(tril(graph), triu(graph))
end


# Permute the vertices of an ordered graph.
function Base.permute!(graph::OrderedGraph, permutation::AbstractVector)
    transpose!(graph.lower, transpose!(graph.upper, symperm(graph.lower, invperm(permutation))))
    graph
end


# Construct a postordered elimination tree.
function etree!(labels::AbstractVector, graph::OrderedGraph)
    tree = etree(graph)
    order = Permutation(tree, DFS())
    permute!(labels, order)
    permute!(graph, order)
    permute!(tree, order)
end


# Construct a postordered supernodal elimination tree.
function stree!(labels::AbstractVector, graph::OrderedGraph, stype::SupernodeType, etree::Tree=etree!(labels, graph))
    rowcount, colcount = supcnt(graph, etree)
    new, ancestor, stree = pothensun(etree, colcount, stype)
    order = Permutation(stree, DFS())
 
    sndptr = Vector{Int}(undef, nv(stree) + 1)
    sepptr = Vector{Int}(undef, nv(stree) + 1)
    sndval = Vector{Int}(undef, nv(graph))
    sndptr[1] = sepptr[1] = 1

    for (i, j) in enumerate(order)
        v = new[j]
        p = sndptr[i]

        while !isnothing(v) && v != ancestor[j]
            sndval[p] = v
            v = parentindex(etree, v)
            p += 1
        end

        sndptr[i + 1] = p
        sepptr[i + 1] = sndptr[i] + sepptr[i] + colcount[new[j]] - p
    end

    permute!(labels, sndval)
    permute!(graph, sndval)
    permute!(stree, order)
    stree, sndptr, sepptr
end


function Base.show(io::IO, ::MIME"text/plain", graph::OrderedGraph)
    SparseArrays._show_with_braille_patterns(io, graph.lower)
end 


############################
# Abstract Graph Interface #
############################


function Graphs.ne(graph::OrderedGraph)
    nnz(graph.lower)
end


function Graphs.nv(graph::OrderedGraph)
    size(graph.lower, 1)
end


function Graphs.inneighbors(graph::OrderedGraph, i::Integer)
    view(rowvals(graph.lower), nzrange(graph.lower, i))
end


function Graphs.outneighbors(graph::OrderedGraph, i::Integer)
    view(rowvals(graph.upper), nzrange(graph.upper, i))
end


function Graphs.is_directed(::Type{OrderedGraph})
    true
end


function Graphs.edgetype(graph::OrderedGraph)
    Graphs.SimpleEdge{Int}
end


function Graphs.has_edge(graph::OrderedGraph, i::Integer, j::Integer)
    i > j ? insorted(i, outneighbors(graph, j)) : false
end
