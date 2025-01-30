module GraphsExt


using Base.Order
using Graphs
using SparseArrays
using StructuredDecompositions.JunctionTrees


function JunctionTrees.permutation(graph::AbstractGraph{I}, alg::Union{RCM, AMD, SymAMD, NodeND, Spectral, BT}) where I
    matrix = spzeros(Bool, I, nv(graph), nv(graph))
    sizehint!(rowvals(matrix), 2ne(graph))

    for v in vertices(graph)
        append!(rowvals(matrix), neighbors(graph, v))
        getcolptr(matrix)[v + 1] = length(rowvals) + 1
    end

    resize!(nonzeros(matrix), length(rowvals(matrix)))
    permutation(matrix, alg)
end


function JunctionTrees.ischordal(graph::AbstractGraph)
    ischordal(vertices(graph)) do v
        neighbors(graph, v)
    end
end


function JunctionTrees.isfilled(graph::AbstractGraph)
    isfilled(vertices(graph)) do v
        neighbors(graph, v)
    end
end


function JunctionTrees.isperfect(graph::AbstractGraph{I}, order::AbstractVector{I}, index::AbstractVector{I}) where I
    # validate arguments
    vertices(graph) != eachindex(order) && throw(ArgumentError("vertices(graph) != eachindex(order)"))

    # run algorithm
    isperfect(order, index) do v
        neighbors(graph, v)
    end
end


function JunctionTrees.mcs(graph::AbstractGraph)
    mcs(vertices(graph)) do v
        neighbors(graph, v)
    end
end


function JunctionTrees.sympermute!(target::SparseMatrixCSC{Nothing, I}, graph::AbstractGraph{I}, index::AbstractVector{I}, order::Ordering) where I
    # validate arguments
    size(target, 1) != nv(graph) && throw(ArgumentError("size(target, 1) != nv(graph)"))
    size(target, 2) != nv(graph) && throw(ArgumentError("size(target, 2) != nv(graph)"))

    # run algorithm
    sympermute!(target, index, order) do v
        neighbors(graph, v)
    end
end


end
