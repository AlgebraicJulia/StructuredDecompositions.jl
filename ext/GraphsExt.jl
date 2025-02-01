module GraphsExt


using Graphs
using StructuredDecompositions.JunctionTrees


function JunctionTrees.Graph(graph::AbstractGraph{V}) where V
    ptr = sizehint!(V[], nv(graph))
    tgt = sizehint!(V[], 2ne(graph))
    push!(ptr, 1)
    
    for v in vertices(graph)
        append!(tgt, neighbors(graph, v))
        push!(ptr, length(tgt) + 1)
    end

    JunctionTrees.Graph(ptr, tgt)
end


end
