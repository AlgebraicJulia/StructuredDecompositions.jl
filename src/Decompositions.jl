module Decompositions

export StructuredDecomposition, StrDecomp, 
      DecompType, Decomposition, CoDecomposition, 
      𝐃, bags, adhesions, adhesionSpans, ∫

using PartialFunctions
using MLStyle

using AbstractTrees
using Base.Threads
using Catlab
using Catlab.CategoricalAlgebra
using Catlab.Graphs
using Catlab.ACSetInterface
using Catlab.CategoricalAlgebra.Diagrams

import Catlab.CategoricalAlgebra.Diagrams: ob_map, hom_map, colimit, limit

#forgetful functor vs: Gr → Set taking G to VG
function vs(X::Graph) FinSet(length(vertices(X))) end
function vs(f::ACSetTransformation) components(f)[1] end

#define the functor skeleton : Set → Skel(Set)
function skeleton(s::FinSet) FinSet(length(s)) end
function skeleton(f::FinFunction)
  (dd, cc) = (dom(f), codom(f))
  #(skel_dom, skel_cod) = (skeleton(dd), skeleton(cc))
  ℓ = isempty(dd) ? Int[] : [findfirst(item -> item == f(x), collect(cc)) for x ∈ collect(dd)]
  FinFunction(ℓ, skeleton(dd), skeleton(cc))
end

export vs, skeleton

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
StrDecomp(the_decomp_shape, the_diagram::FinDomFunctor) = StrDecomp(the_decomp_shape, the_diagram, Decomposition)
# XXX I've added a type here to minimize conflict with the JunctionTrees interop.

"""If we want to explicitly specify the decomposition type of a structured decomposition, then this can be done by explicitly passing a some t::DecompType as an argument to the constructor.
DecompType has two values: Decomposition and CoDecomposition. These are used to distinguish whether we think of a C-valued structured decomposition d: FG → Span C as a diagram into Span C or whether 
we think of it as a diagram into Cospan C of the form d: FG → Cospan C^op. 
"""
function StrDecomp(the_decomp_shape, the_diagram, the_decomp_type)
  d  = StrDecomp(the_decomp_shape, the_diagram, the_decomp_type, dom(the_diagram))
  dc = s -> @match the_decomp_type begin
    Decomposition   => dom(s[1])   == dom(s[2])
    CoDecomposition => codom(s[1]) == codom(s[2])
  end
  all(dc, adhesionSpans(d)) ? d : throw(StrDecompError(d, the_decomp_type))
end
#construct a structured decomposition and check whether the decomposition shape actually makes sense. 
# TODO: check that the domain is also correct...

ob_map(d::StructuredDecomposition, x)  = ob_map(d.diagram, x)
hom_map(d::StructuredDecomposition, f) = hom_map(d.diagram, f)

colimit(d::StructuredDecomposition) = colimit(FreeDiagram(d.diagram))
limit(d::StructuredDecomposition) = limit(FreeDiagram(d.diagram))

struct StrDecompError <: Exception
    d::StrDecomp
    decomptype::DecompType
end

Base.showerror(io::IO, e::StrDecompError) = print(io, "$(str(e.d)) is not a $(string(e.decomptype))")

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

# Get the points in el:Elements corresp to either Vertices of Edges
function get_components(el::Elements, j)
    filter(part -> el[part, :πₑ] == j, parts(el, :El))
end

# Get the points in el:Elements into actual objects of the Category ∫G (for G the shape of decomp)
function get_components(d::StructuredDecomposition, el::Elements, j)
    getindex.(Ref(ob_generators(d.domain)), get_components(el, j))
end

function getFromDom(c::ShapeCpt, d::StructuredDecomposition, el::Elements = elements(d.decomp_shape))
  @match c begin
    ShapeVertex => get_components(d, el, 1)
    ShapeEdge   => get_components(d, el, 2)
    ShapeSpan   => map(get_components(el, 2)) do elpart
        filter(parts(el, :Arr)) do part
            el[part, :src] == elpart
        end
    end
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

function map_index(f::Function, x; indexing::Bool=true)
    indexing ? collect(zip(x, map(f, x))) : map(f, x)
end

function get(c::StrDcmpCpt, d::StructuredDecomposition, indexing::Bool)
  # Either apply the object- or the morphism component of the diagram of d
  el = elements(d.decomp_shape) 
  get_cat_cpt_of_flavor(sc::ShapeCpt) = getFromDom(sc, d, el)

  # XXX map_ind(f, x) = indexing == true ? collect(zip(x, map(f, x))) : map(f, x)
  # now just do the actual computation
  @match c begin 
    Bag          => map_index(x -> ob_map(d.diagram, x), get_cat_cpt_of_flavor(ShapeVertex); indexing=indexing) 
    AdhesionApex => map_index(x -> ob_map(d.diagram, x), get_cat_cpt_of_flavor(ShapeEdge); indexing=indexing) 
    AdhesionSpan => map_index(y -> map(x -> hom_map(d.diagram, x), y), get_cat_cpt_of_flavor(ShapeSpan); indexing=indexing) 
  end
end

"""Get a vector of indexed bags; i.e. a vector consisting of pairs (x, dx) where x ∈ FG and dx is the value mapped to x under the decompositon d"""
bags(d, ind) = get(Bag, d, ind)

"""Get a vector of the bags of a decomposition"""
bags(d) = bags(d, false)

"""
Get a vector of indexed adhesions; i.e. a vector consisting of pairs (e, de) where e is an edge in ∫G and de is the value mapped to e under the decompositon d
"""
adhesions(d, ind) = get(AdhesionApex, d, ind)

"""Get a vector of the adhesions of a decomposition"""
adhesions(d) = adhesions(d, false)     

"""
Get a vector of indexed adhesion spans; i.e. a vector consisting of pairs (x₁ <- e -> x₂, dx₁ <- de -> dx₂) where x₁ <- e -> x₂ is span in ∫G and dx₁ <- de -> dx₂ is what its image under the decompositon d
"""
adhesionSpans(d, ind) = get(AdhesionSpan, d, ind)

"""Get a vector of the adhesion spans of a decomposition"""
adhesionSpans(d) = adhesionSpans(d, false)

function elements_graph(el::Elements)
  F = FinFunctor(Dict(:V => :El, :E => :Arr), Dict(:src => :src, :tgt => :tgt), SchGraph, SchElements)
  ΔF = DeltaMigration(F)
  return migrate(Graph, el, ΔF)
end

"""
Syntactic sugar for costructing the category of elements of a graph. 
Note that ∫(G) has type Category whereas elements(G) has type Elements
"""
function ∫(G::T) where {T <: ACSet} ∫(elements(G))            end 
function ∫(G::Elements)             FinCat(elements_graph(G)) end 

# Reverse direction of the edges
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

When we think of the functor F as a computational problem (taking inputs in C to solution spaces in E), then 𝐃f should be thought of as lifting the global comuputation of F to local computation on the constituent parts of C-valued decompositions. 

In particular, given a structured decomposition d: FG → C and a sheaf F: C → FinSet^{op} w.r.t to the decompositon topology, we can make a structured decomposition valued in FinSet^{op} by lifting the sheaf to a functor 𝐃_f: 𝐃C → 𝐃(S^{op}) between categories of structured decompositions. 
"""
function 𝐃(f, d::StructuredDecomposition, t::DecompType = d.decomp_type)::StructuredDecomposition 
  flip = t == d.decomp_type ?  x -> x : FinCat ∘ op_graph ∘ graph #work with ( ∫G )^{op}
  δ   = d.diagram 
  X   = dom(δ)
  # Q is the composite Q = F ∘ d : FG → C → S^{op}
  Q   = FinDomFunctor(
          Dict(x => f(ob_map(δ,x))   for x ∈ ob_generators(X) ), #the ob  map Q₀
          Dict(g => f(hom_map(δ, g)) for g ∈ hom_generators(X)), #the hom map Q₁
          flip(X))
  StrDecomp(d.decomp_shape, Q, t) 
end

end
