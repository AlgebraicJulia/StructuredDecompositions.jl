"""
    StrDecomp(graph::HasGraph;
        alg::PermutationOrAlgorithm=DEFAULT_ELIMINATION_ALGORITHM,
        snd::SupernodeType=DEFAULT_SUPERNODE_TYPE)

Construct a structured decomposition of a simple graph. See [`cliquetree`](@ref) for the meaning of `alg` and `snd`.
"""
function StrDecomp(graph::HasGraph; alg::PermutationOrAlgorithm=DEFAULT_ELIMINATION_ALGORITHM, snd::SupernodeType=DEFAULT_SUPERNODE_TYPE)
    label, tree = cliquetree(symmetrize(graph); alg, snd)
    n = length(tree)
    shape = Graph(n)
    
    for i in 1:n
        j = parentindex(tree, i)
        
        if !isnothing(j)
            add_edge!(shape, j, i)
        end
    end

    diagram = FinDomFunctor(homomorphisms(graph, label, tree)..., ∫(shape))
    StrDecomp(shape, diagram, Decomposition, dom(diagram))
end


function homomorphisms(graph::HasGraph, label::AbstractVector, tree::CliqueTree)
    m = length(collect(rootindices(tree)))
    n = length(tree)
    relative = relatives(tree)
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
            homomorphism[k] = induced_homomorphism(subgraph[n + k], subgraph[j], CliqueTrees.neighbors(relative, i))

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


function symmetrize(graph::AbstractSymmetricGraph)
    graph
end


function symmetrize(graph::HasGraph)
    symmetric = SymmetricGraph(nv(graph))
    add_edges!(symmetric, src(graph), tgt(graph))
    symmetric
end
