module JunctionTrees

import AMD
import CuthillMcKee
import LinkedLists
import Metis
import TreeWidthSolver

using Base.Threads
using AbstractTrees
using DataStructures
using SparseArrays

using Catlab
using Catlab.CategoricalAlgebra
using Catlab.BasicGraphs
using Catlab.Graphs
using Catlab.ACSetInterface
using Catlab.CategoricalAlgebra.Diagrams

import Catlab.CategoricalAlgebra.Diagrams: ob_map, hom_map, colimit, limit

using ..Decompositions
import ..Decompositions: StrDecomp

include("junction_trees/orders.jl")
include("junction_trees/elimination_algorithms.jl")
include("junction_trees/ordered_graphs.jl")
include("junction_trees/trees.jl")
include("junction_trees/postorder_trees.jl")
include("junction_trees/elimination_trees.jl")
include("junction_trees/supernode_types.jl")
include("junction_trees/supernode_trees.jl")
include("junction_trees/junction_trees.jl")

##################################
# Integration with JunctionTrees #
##################################

"""
    StrDecomp(graph::AbstractSymmetricGraph[, ealg::Union{Order, EliminationAlgorithm}[, stype::SupernodeType]])

Construct a structured decomposition of a simple graph, optionally specifying an elimination algorithm and
supernode type.
"""
function StrDecomp(
    graph::AbstractSymmetricGraph,
    ealg::Union{Order, EliminationAlgorithm}=DEFAULT_ELIMINATION_ALGORITHM,
    stype::SupernodeType=DEFAULT_SUPERNODE_TYPE)

    merge_decompositions(decompositions(graph, ealg, stype))
end

# Construct a tree decomposition.
# ----------------------------------------
#    graph    simple connected graph
#    jtree    junction tree
# ----------------------------------------
function StrDecomp(graph::AbstractSymmetricGraph, jtree::JunctionTree)
    n = length(jtree)
    tree = Graph(n)
    
    for i in 1:n - 1
        add_edge!(tree, i, parentindex(jtree, i))
    end

    diagram = FinDomFunctor(homomorphisms(graph, jtree)..., ∫(tree))
    StrDecomp(tree, diagram, Decomposition, dom(diagram))
end
# TODO this function needs a test

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


function decompositions(graph::AbstractSymmetricGraph, ealg::EliminationAlgorithm, stype::SupernodeType)
    component = connected_components(graph)

    n = length(component)
    decomposition = Vector(undef, n)
   
    # TODO in this case, a Junction Tree
    @threads for i in 1:n
        subgraph = induced_subgraph(graph, component[i])
        decomposition[i] = StrDecomp(subgraph, JunctionTree(subgraph, ealg, stype))
    end
    
    decomposition
end


function decompositions(graph::AbstractSymmetricGraph, order::Order, stype::SupernodeType)
    component = connected_components(graph)

    n = length(component)
    decomposition = Vector(undef, n)
    
    @threads for i in 1:n
        subgraph = induced_subgraph(graph, component[i])
        decomposition[i] = StrDecomp(subgraph, JunctionTree(subgraph, induced_order(order, component[i]), stype))
    end
    
    decomposition
end


function homomorphisms(graph::AbstractSymmetricGraph, jtree::JunctionTree)
    n = length(jtree)
    subgraph = Vector{Any}(undef, 2n - 1)
    homomorphism = Vector{Any}(undef, 2n - 2)
    
    for i in 1:n
        # clique(i)
        subgraph[i] = induced_subgraph(graph, clique(jtree, i))
    end
  
    for i in 1:n - 1
        # seperator(i)
        subgraph[n + i] = induced_subgraph(graph, seperator(jtree, i))
    end
  
    for i in 1:n - 1
        # seperator(i) → clique(parent(i))
        j = parentindex(jtree, i)
        homomorphism[i] = induced_homomorphism(subgraph[n + i], subgraph[j], seperator_to_parent(jtree, i))
    end
    
    for i in 1:n - 1
        # seperator(i) → clique(i)
        homomorphism[n + i - 1] = induced_homomorphism(subgraph[n + i], subgraph[i], seperator_to_self(jtree, i))
    end
    
    subgraph, homomorphism
end


function induced_order(order::Order, elements::AbstractVector)
    Order(sortperm(inverse(order, elements)))
end


function induced_homomorphism(domain::AbstractSymmetricGraph, codomain::AbstractSymmetricGraph, V::AbstractVector)
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


end
