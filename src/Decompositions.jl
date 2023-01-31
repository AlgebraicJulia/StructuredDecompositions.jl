module Decompositions

export StructuredDecomposition, StrDecomp, 
      DecompType, Decomposition, CoDecomposition, 
      𝐃, bags, adhesions, adhesionSpans, 
      ∫ #, op_graph, codecomp

using PartialFunctions
using MLStyle

using Catlab
using Catlab.CategoricalAlgebra
using Catlab.Graphs
using Catlab.ACSetInterface
using Catlab.CategoricalAlgebra.Diagrams
import Catlab.CategoricalAlgebra.Diagrams: ob_map, hom_map, colimit, limit

#####################
#   DATA
#####################
"""Structured decompositions
"""
abstract type StructuredDecomposition{G, C, D} <: Diagram{id, C, D} end

@data DecompType begin
  Decomposition 
  CoDecomposition
end

"""Structrured decomposition struct
    -- think of these are graphs whose vertices are labeled by the objects of some category 
    and whose edges are labeled by SPANS in this category
"""
struct StrDecomp{G, C, D} <: StructuredDecomposition{G, C, D}  
  decomp_shape ::G 
  diagram      ::D             
  decomp_type  ::DecompType
  domain       ::C  
end

StrDecomp(the_decomp_shape, the_diagram) = StrDecomp(the_decomp_shape, the_diagram, Decomposition, dom(the_diagram))

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
function limit(  d::StructuredDecomposition) limit(FreeDiagram(d.diagram))   end



# BEGIN UTILS
"""Structured decomposition Utils"""
# ShapeVertex is the "vertex-objects" of ∫G; i.e in the form (V, v)
#  ShapeEdge is the "edge-objects" of ∫G; i.e in the form (E, e)
#  ShapeSpan is the "span-objects" of ∫G; i.e. in the form (V,x) <-- (E,e=xy) --> (V,y)
@data ShapeCpt begin
  ShapeVertex 
  ShapeEdge   
  ShapeSpan   
end

function getFromDom(c::ShapeCpt, d::StructuredDecomposition, el::Elements = elements(d.decomp_shape))
  #get the points in el:Elements corresp to either Vertices of Edges
  get_Ele_cpt(j) = filter(part -> el[part, :πₑ] == j, parts(el, :El)) 
  #get the points in el:Elements into actual objects of the Category ∫G (for G the shape of decomp)
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
  # either apply the object- or the morphism component of the diagram of d
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

bags(d, ind)          = get(Bag, d, ind)
bags(d)               = bags(d, false)
adhesions(d, ind)     = get(AdhesionApex, d, ind)
adhesions(d)          = adhesions(d, false)     
adhesionSpans(d, ind) = get(AdhesionSpan, d, ind)
adhesionSpans(d)      = adhesionSpans(d, false)

function elements_graph(el::Elements)
  F = FinFunctor(Dict(:V => :El, :E => :Arr), Dict(:src => :src, :tgt => :tgt), SchGraph, SchElements)
  ΔF = DeltaMigration(F, Elements{Symbol}, Graph)
  return ΔF(el)
end

"""∫(G) has type Category whereas elements(G) has type Elements
"""
function ∫(G::T) where {T <: ACSet} ∫(elements(G))            end 
function ∫(G::Elements)             FinCat(elements_graph(G)) end 

#reverse direction of the edges
function op_graph(g::Graph)::Graph
  F = FinFunctor(Dict(:V => :V, :E => :E), Dict(:src => :tgt, :tgt => :src), SchGraph, SchGraph)
  ΔF = DeltaMigration(F, Graph, Graph)
  return ΔF(g)
end

"""Given a structured decomposition d: FG → C and a sheaf F: C → S^{op} w.r.t to the decompositon topology, 
we can make a structured decomposition valued in S^{op} which has the relevant data by composition:
F ∘ d : FG → C → S^{op}.
This is done by first lifting the sheaf to a functor 𝐃_f: 𝐃C → 𝐃(S^{op}) between categories of structured decompositions. 
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

end
