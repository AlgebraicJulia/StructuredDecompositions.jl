"""
    StrDecomp(graph::HasGraph;
        alg::PermutationOrAlgorithm=DEFAULT_ELIMINATION_ALGORITHM,
        snd::SupernodeType=DEFAULT_SUPERNODE_TYPE)

Construct a structured decomposition of a simple graph. See [`junctiontree`](@ref) for the meaning of `alg` and `snd`.
"""
function StrDecomp(graph::HasGraph; alg::PermutationOrAlgorithm=DEFAULT_ELIMINATION_ALGORITHM, snd::SupernodeType=DEFAULT_SUPERNODE_TYPE)
    label, tree = junctiontree(adjacency_matrix(graph); alg, snd)
    relative!(tree)
    n = length(tree)
    shape = Graph(n)
    
    for i in 1:n - 1
        add_edge!(shape, parentindex(tree, i), i)
    end

    diagram = FinDomFunctor(homomorphisms(graph, label, tree)..., ∫(shape))
    StrDecomp(shape, diagram, Decomposition, dom(diagram))
end


function homomorphisms(graph::HasGraph, label::AbstractVector, tree::JunctionTree)
    m = length(collect(rootindices(tree)))
    n = length(tree)
    subgraph = Vector{Any}(undef, 2n - m)
    homomorphism = Vector{Any}(undef, 2n - 2m)
    
    for i in 1:n
        # bag(i)
        subgraph[i] = induced_subgraph(graph, view(label, getindex(tree, i)))
    end
 
    k = 0

    for i in 1:n
        j = parentindex(tree, i)

        if !isnothing(j)
            k += 1

            # separator(i)
            subgraph[n + k] = induced_subgraph(graph, view(label, separator(tree, i)))

            # separator(i) → bag(j)
            homomorphism[k] = induced_homomorphism(subgraph[n + k], subgraph[j], relative(tree, i))

            # separator(i) → bag(i)
            homomorphism[n - m + k] = induced_homomorphism(subgraph[n + k], subgraph[i], length(residual(tree, i)) .+ eachindex(separator(tree, i)))
        end
    end
    
    subgraph, homomorphism
end


function induced_homomorphism(domain::HasGraph, codomain::HasGraph, V::AbstractVector)
    index = Dict{Tuple{Int, Int}, Int}()
    sizehint!(index, ne(codomain))
    
    for e in edges(codomain)
        index[src(codomain, e), tgt(codomain, e)] = e
    end
    
    E = Vector{Int}(undef, ne(domain))
    
    for e in edges(domain)
        E[e] = index[V[src(domain, e)], V[tgt(domain, e)]]
    end
    
    ACSetTransformation(domain, codomain; V, E)
end


# Construct the adjacency matrix of an undirected graph.
function adjacency_matrix(graph::HasGraph)
    symmetric = SymmetricGraph(nv(graph))
    add_edges!(symmetric, src(graph), tgt(graph))
    adjacency_matrix(symmetric)
end


function adjacency_matrix(graph::AbstractSymmetricGraph)
    matrix = spzeros(Bool, Int, nv(graph), nv(graph))
    sizehint!(rowvals(matrix), ne(graph))

    for v in vertices(graph)
        append!(rowvals(matrix), neighbors(graph, v))
        getcolptr(matrix)[v + 1] = length(rowvals(matrix)) + 1
    end

    resize!(nonzeros(matrix), ne(graph))
    fill!(nonzeros(matrix), true)
    matrix
end
