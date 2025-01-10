"""
    StrDecomp(graph::HasGraph;
        alg::PermutationOrAlgorithm=AMD(),
        snd::SupernodeType=Maximal())

Construct a structured decomposition of a simple graph. See [`junctiontree!`](@ref) for the meaning of `alg` and `snd`.
"""
function StrDecomp(graph::HasGraph;
    alg::PermutationOrAlgorithm=DEFAULT_ELIMINATION_ALGORITHM,
    snd::SupernodeType=DEFAULT_SUPERNODE_TYPE)

    merge_decompositions(decompositions(graph, alg, snd))
end


# Construct a tree decomposition.
function StrDecomp(graph::HasGraph, label::AbstractVector, tree::JunctionTree)
    n = length(tree)
    shape = Graph(n)
    
    for i in 1:n - 1
        add_edge!(shape, parentindex(tree, i), i)
    end

    diagram = FinDomFunctor(homomorphisms(graph, label, tree)..., ∫(shape))
    StrDecomp(shape, diagram, Decomposition, dom(diagram))
end


function merge_decompositions(decomposition::AbstractVector)      
    tree = apex(coproduct(map(d -> d.decomp_shape, decomposition)))
    l = length(decomposition)
    m = nv(tree)
    subgraph = Vector(undef, 2m - l)
    homomorphism = Vector(undef, 2m - 2l)

    i = 0

    for j in 1:l
        n = nv(decomposition[j].decomp_shape)

        for k in 1:n
            subgraph[i + k] = ob_map(decomposition[j].diagram, k)
        end

        for k in 1:n - 1
            subgraph[i - j + k + m + 1] = ob_map(decomposition[j].diagram, k + n)
        end

        for k in 1:n - 1
            homomorphism[i - j + k + 1] = hom_map(decomposition[j].diagram, k)
        end

        for k in 1:n - 1
            homomorphism[i - j + k - l + m + 1] = hom_map(decomposition[j].diagram, k + n - 1)
        end

        i += n
    end  

    diagram = FinDomFunctor(subgraph, homomorphism, ∫(tree))
    StrDecomp(tree, diagram, Decomposition, dom(diagram))
end


function decompositions(graph::HasGraph, alg::EliminationAlgorithm, snd::SupernodeType)
    component = connected_components(graph)

    n = length(component)
    decomposition = Vector(undef, n)
    
    @threads for i in 1:n
        subgraph = induced_subgraph(graph, component[i])
        decomposition[i] = StrDecomp(subgraph, junctiontree!(adjacency_matrix(subgraph); alg, snd)...)
    end
    
    decomposition
end


function decompositions(graph::HasGraph, alg::AbstractVector, snd::SupernodeType)
    component = connected_components(graph)

    n = length(component)
    decomposition = Vector(undef, n)
    
    @threads for i in 1:n
        subgraph = induced_subgraph(graph, component[i])
        decomposition[i] = StrDecomp(subgraph, junctiontree!(adjacency_matrix(subgraph); alg=induced_order(invperm(alg), component[i]), snd)...)
    end
    
    decomposition
end


function homomorphisms(graph::HasGraph, label::AbstractVector, tree::JunctionTree)
    n = length(tree)
    subgraph = Vector{Any}(undef, 2n - 1)
    homomorphism = Vector{Any}(undef, 2n - 2)
    
    for i in 1:n
        # bag(i)
        subgraph[i] = induced_subgraph(graph, view(label, getindex(tree, i)))
    end
 
    for i in 1:n - 1 
        # separator(i)
        subgraph[n + i] = induced_subgraph(graph, view(label, separator(tree, i)))
    end
 
    for i in 1:n - 1 
        # separator(i) → bag(j)
        j = parentindex(tree, i)
        homomorphism[i] = induced_homomorphism(subgraph[n + i], subgraph[j], relative(tree, i))
    end

    for i in 1:n - 1    
        # separator(i) → bag(i)
        j = parentindex(tree, i)
        homomorphism[n + i - 1] = induced_homomorphism(subgraph[n + i], subgraph[i], length(residual(tree, i)) .+ eachindex(separator(tree, i)))
    end
    
    subgraph, homomorphism
end


function induced_order(index::AbstractVector, elements::AbstractVector)
    sortperm(elements; by=i -> index[i])
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
function adjacency_matrix(graph::HasGraph, neighbors::Function, ne::Function)
    rowval = Vector{Int}(undef, 2ne(graph))
    colptr = Vector{Int}(undef, nv(graph) + 1)
    colptr[1] = 1

    for i in vertices(graph)
        p = colptr[i]

        for j in neighbors(graph, i)
            rowval[p] = j
            p += 1
        end

        colptr[i + 1] = p
        sort!(rowval, colptr[i], colptr[i + 1] - 1, DEFAULT_STABLE, ForwardOrdering())
    end

    nzval = ones(Bool, 2ne(graph))
    SparseMatrixCSC(nv(graph), nv(graph), colptr, rowval, nzval)
end


function adjacency_matrix(graph::AbstractGraph)
    adjacency_matrix(graph, all_neighbors, ne)
end


function adjacency_matrix(graph::AbstractSymmetricGraph)
    adjacency_matrix(graph, neighbors, g -> ne(g) ÷ 2)
end
