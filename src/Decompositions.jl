module Decompositions

export StructuredDecomposition, StrDecomp, 
      DecompType, Decomposition, CoDecomposition, 
      𝐃, bags, adhesions, adhesionSpans, 
      ∫

using ..JunctionTrees
using ..JunctionTrees: EliminationAlgorithm, SupernodeType, DEFAULT_ELIMINATION_ALGORITHM, DEFAULT_SUPERNODE_TYPE

using PartialFunctions
using MLStyle

using AbstractTrees
using Base: DEFAULT_STABLE, ForwardOrdering
using Base.Threads
using Catlab
using Catlab.CategoricalAlgebra
using Catlab.Graphs
using Catlab.ACSetInterface
using Catlab.CategoricalAlgebra.Diagrams
using SparseArrays

import Catlab.CategoricalAlgebra.Diagrams: ob_map, hom_map, colimit, limit


#####################
#   DATA
#####################
"""Structured decompositions are diagrams.
"""
abstract type StructuredDecomposition{G, C, D} <: Diagram{id, C, D} end

@data DecompType begin
  Decomposition 
  CoDecomposition
end

"""Structured decomposition struct
    -- Consider these as graphs whose vertices are labeled by the objects of some category 
    and whose edges are labeled by spans in this category
"""
struct StrDecomp{G, C, D} <: StructuredDecomposition{G, C, D}  
  decomp_shape ::G 
  diagram      ::D             
  decomp_type  ::DecompType
  domain       ::C  
end

"""One can construct a structured decomposition by simply providing the graph representing the shape of the decompostion and the relevant diagram.
This constructor will default to setting the Decompsotion Type to Decomposition (i.e. we default to viewing C-valued structured decompositions as diagrams into Span C)
"""
StrDecomp(the_decomp_shape, the_diagram) = StrDecomp(the_decomp_shape, the_diagram, Decomposition)

"""If we want to explicitly specify the decomposition type of a structured decomposition, then this can be done by explicitly passing a some t::DecompType as an argument to the constructor.
DecompType has two values: Decomposition and CoDecomposition. These are used to distinguish whether we think of a C-valued structured decomposition d: FG → Span C as a diagram into Span C or whether 
we think of it as a diagram into Cospan C of the form d: FG → Cospan C^op. 
"""
#construct a structured decomposition and check whether the decomposition shape actually makes sense. 
# TODO: check that the domain is also correct...
function StrDecomp(the_decomp_shape, the_diagram, the_decomp_type)
  d  = StrDecomp(the_decomp_shape, the_diagram, the_decomp_type, dom(the_diagram))
  dc = s -> @match the_decomp_type begin
    Decomposition   => dom(s[1])   == dom(s[2])
    CoDecomposition => codom(s[1]) == codom(s[2])
  end
  if all(dc, adhesionSpans(d))
   return d
  else 
    error(str(d) * " is not a " * string(the_decomp_type))
  end
end

ob_map(d::StructuredDecomposition, x)  = ob_map(d.diagram, x)
hom_map(d::StructuredDecomposition, f) = hom_map(d.diagram, f)

function colimit(d::StructuredDecomposition) colimit(FreeDiagram(d.diagram)) end
function limit(d::StructuredDecomposition) limit(FreeDiagram(d.diagram)) end



# BEGIN UTILS
#=Structured decomposition Utils=#
# ShapeVertex is the "vertex-objects" of ∫G; i.e in the form (V, v)
#  ShapeEdge is the "edge-objects" of ∫G; i.e in the form (E, e)
#  ShapeSpan is the "span-objects" of ∫G; i.e. in the form (V,x) <-- (E,e=xy) --> (V,y)
@data ShapeCpt begin
  ShapeVertex 
  ShapeEdge   
  ShapeSpan   
end

function getFromDom(c::ShapeCpt, d::StructuredDecomposition, el::Elements = elements(d.decomp_shape))
  #Get the points in el:Elements corresp to either Vertices of Edges
  get_Ele_cpt(j) = filter(part -> el[part, :πₑ] == j, parts(el, :El)) 
  #Get the points in el:Elements into actual objects of the Category ∫G (for G the shape of decomp)
  get_Cat_cpt(j) = map(grCpt -> ob_generators(d.domain)[grCpt], get_Ele_cpt(j))
  @match c begin
    ShapeVertex => get_Cat_cpt(1)
    ShapeEdge   => get_Cat_cpt(2)
    ShapeSpan   => map( epart -> filter(part -> el[part, :src] == epart,  parts(el, :Arr)), get_Ele_cpt(2))
  end
end

@data MapType begin
  ObMap
  HomMap
end

@data StrDcmpCpt begin
  Bag
  AdhesionApex
  AdhesionSpan
end

function get(c::StrDcmpCpt, d::StructuredDecomposition, indexing::Bool)
  # Either apply the object- or the morphism component of the diagram of d
  evalDiagr(t::MapType, x) = @match t begin 
    ObMap  => ob_map( d.diagram, x) 
    HomMap => hom_map(d.diagram, x)
  end
  el = elements(d.decomp_shape) 
  get_cat_cpt_of_flavor(sc::ShapeCpt) = getFromDom(sc, d, el)

  map_ind(f, x) = indexing == true ? collect(zip(x, map(f, x))) : map(f, x)
  # now just do the actual computation
  @match c begin 
    Bag          => map_ind(evalDiagr $ ObMap     , get_cat_cpt_of_flavor(ShapeVertex) ) 
    AdhesionApex => map_ind(evalDiagr $ ObMap     ,get_cat_cpt_of_flavor(ShapeEdge)   ) 
    AdhesionSpan => map_ind(map $ (evalDiagr $ HomMap), get_cat_cpt_of_flavor(ShapeSpan)   ) 
  end
end

"""Get a vector of indexed bags; i.e. a vector consisting of pairs (x, dx) where x ∈ FG and dx is the value mapped to x under the decompositon d"""
bags(d, ind)          = get(Bag, d, ind)
"""Get a vector of the bags of a decomposition"""
bags(d)               = bags(d, false)

"""Get a vector of indexed adhesions; i.e. a vector consisting of pairs (e, de) where e is an edge in ∫G and de is the value mapped to e under the decompositon d"""
adhesions(d, ind)     = get(AdhesionApex, d, ind)
"""Get a vector of the adhesions of a decomposition"""
adhesions(d)          = adhesions(d, false)     

"""Get a vector of indexed adhesion spans; i.e. a vector consisting of pairs (x₁ <- e -> x₂, dx₁ <- de -> dx₂) where x₁ <- e -> x₂ is span in ∫G and dx₁ <- de -> dx₂ is what its image under the decompositon d"""
adhesionSpans(d, ind) = get(AdhesionSpan, d, ind)
"""Get a vector of the adhesion spans of a decomposition"""
adhesionSpans(d)      = adhesionSpans(d, false)

function elements_graph(el::Elements)
  F = FinFunctor(Dict(:V => :El, :E => :Arr), Dict(:src => :src, :tgt => :tgt), SchGraph, SchElements)
  ΔF = DeltaMigration(F)
  return migrate(Graph, el, ΔF)
end

"""Syntactic sugar for costructing the category of elements of a graph. 
Note that ∫(G) has type Category whereas elements(G) has type Elements
"""
function ∫(G::T) where {T <: ACSet} ∫(elements(G))            end 
function ∫(G::Elements)             FinCat(elements_graph(G)) end 

#Reverse direction of the edges
function op_graph(g::Graph)::Graph
  F = FinFunctor(Dict(:V => :V, :E => :E), Dict(:src => :tgt, :tgt => :src), SchGraph, SchGraph)
  ΔF = DeltaMigration(F)
  return migrate(Graph, g, ΔF)
end

"""
The construction of categories of structured decompostions is functorial; 
it consists of a functor 𝐃: Cat_{pullback} → Cat taking any category C with pullbacks to the category 
𝐃C of C-values structured decompositions. The functoriality of this construction allows us to lift any functor 
F : C → E to a functor 𝐃f : 𝐃 C → 𝐃 E which maps C-valued structured decompositions to E-valued structured decompositions. 
When we think of the functor F as a computational problem (taking inputs in C to solution spaces in E), then 𝐃f should be thought
of as lifting the global comuputation of F to local computation on the constituent parts of C-valued decompositions. 
In particular, given a structured decomposition d: FG → C and a sheaf F: C → FinSet^{op} w.r.t to the decompositon topology, 
we can make a structured decomposition valued in FinSet^{op} by lifting the sheaf to a functor 𝐃_f: 𝐃C → 𝐃(S^{op}) between 
categories of structured decompositions. 
"""
function 𝐃(f, d ::StructuredDecomposition, t::DecompType = d.decomp_type)::StructuredDecomposition 
  flip = t == d.decomp_type ?  x -> x : FinCat ∘ op_graph ∘ graph #work with ( ∫G )^{op}
  δ   = d.diagram 
  X   = dom(δ)
  #Q is the composite Q = F ∘ d : FG → C → S^{op}
  Q   = FinDomFunctor(
          Dict(x => f(ob_map(δ,x))   for x ∈ ob_generators(X) ), #the ob  map Q₀
          Dict(g => f(hom_map(δ, g)) for g ∈ hom_generators(X)), #the hom map Q₁
          flip(X)
        )
  StrDecomp(d.decomp_shape, Q, t) 
end


##################################
# Integration with JunctionTrees #
##################################


"""
    StrDecomp(graph::AbstractSymmetricGraph[, alg::Union{Order, EliminationAlgorithm}[, type::SupernodeType]])

Construct a structured decomposition of a simple graph, optionally specifying an elimination algorithm and
supernode type.
"""
function StrDecomp(
    graph::HasGraph,
    alg::Union{Permutation, EliminationAlgorithm}=DEFAULT_ELIMINATION_ALGORITHM,
    type::SupernodeType=DEFAULT_SUPERNODE_TYPE)

    merge_decompositions(decompositions(graph, alg, type))
end


# Construct a tree decomposition.
function StrDecomp(graph::HasGraph, order::Permutation, tree::JunctionTree)
    n = length(tree)
    shape = Graph(n)
    
    for i in 1:n - 1
        add_edge!(shape, parentindex(tree, i), i)
    end

    diagram = FinDomFunctor(homomorphisms(graph, order, tree)..., ∫(shape))
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


function decompositions(graph::HasGraph, alg::EliminationAlgorithm, type::SupernodeType)
    component = connected_components(graph)

    n = length(component)
    decomposition = Vector(undef, n)
    
    @threads for i in 1:n
        subgraph = induced_subgraph(graph, component[i])
        decomposition[i] = StrDecomp(subgraph, junctiontree(adjacency_matrix(subgraph), alg, type)...)
    end
    
    decomposition
end


function decompositions(graph::HasGraph, order::Permutation, type::SupernodeType)
    component = connected_components(graph)

    n = length(component)
    decomposition = Vector(undef, n)
    
    @threads for i in 1:n
        subgraph = induced_subgraph(graph, component[i])
        decomposition[i] = StrDecomp(subgraph, junctiontree(adjacency_matrix(subgraph), induced_order(order, component[i]), type)...)
    end
    
    decomposition
end


function homomorphisms(graph::HasGraph, order::Permutation, tree::JunctionTree)
    n = length(tree)
    subgraph = Vector{Any}(undef, 2n - 1)
    homomorphism = Vector{Any}(undef, 2n - 2)
    
    for i in 1:n
        # bag(i)
        subgraph[i] = induced_subgraph(graph, view(order, getindex(tree, i)))
    end
 
    for i in 1:n - 1 
        # separator(i)
        subgraph[n + i] = induced_subgraph(graph, view(order, separator(tree, i)))
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


function induced_order(order::Permutation, elements::AbstractVector)
    Permutation(sortperm(view(invperm(order), elements)))
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
    function neighbors(graph, v)
        i = view(view(graph, :tgt), incident(graph, v, :src))
        o = view(view(graph, :src), incident(graph, v, :tgt))
        Iterators.flatten((i, o))
    end

    adjacency_matrix(graph, neighbors, ne)
end


function adjacency_matrix(graph::AbstractSymmetricGraph)
    function neighbors(graph, v)
        view(view(graph, :tgt), incident(graph, v, :src))
    end

    function ne(graph)
        nparts(graph, :E) ÷ 2
    end

    adjacency_matrix(graph, neighbors, ne)
end


end
