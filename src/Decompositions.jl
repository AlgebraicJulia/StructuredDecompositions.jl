module Decompositions

export StructuredDecomposition, StrDecomp, CoStrDecomp, ∫, DecompType, Decomposition, CoDecomposition, 𝐃, bags, adhesions, adhesionSpans, op_graph, codecomp, GrCpt, Vertex, Edge

using PartialFunctions
using MLStyle

using Catlab
using Catlab.CategoricalAlgebra
using Catlab.Graphs
using Catlab.ACSetInterface
using Catlab.CategoricalAlgebra.Diagrams
#=
using Catlab.Theories
using Catlab.CategoricalAlgebra.FinSets
using Catlab.CategoricalAlgebra.Diagrams
using Catlab.CategoricalAlgebra.FreeDiagrams
using ..FinCats: FreeCatGraph, FinDomFunctor, collect_ob, collect_hom, Categories.op, graph
using Catlab.Programs
using Catlab.Graphics
=#

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
  domain       ::C
  diagram      ::D               
  decomp_type  ::DecompType
end

StrDecomp(the_decomp_shape, the_domain, the_diagram) = StrDecomp(the_decomp_shape, the_domain, the_diagram, Decomposition)

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

@data Indexing begin
  Indexed 
  NotIndexed
end

function get(c::StrDcmpCpt, d::StructuredDecomposition, ind::Indexing)
  # either apply the object- or the morphism component of the diagram of d
  evalDiagr(t::MapType, x) = @match t begin 
    ObMap  => ob_map( d.diagram, x) 
    HomMap => hom_map(d.diagram, x)
  end
  el = elements(d.decomp_shape) 
  get_cat_cpt_of_flavor(sc::ShapeCpt) = getFromDom(sc, d, el)

  map_ind(f, x) = @match ind begin
    Indexed    => Dict(zip( map(f, x), x) )
    NotIndexed =>           map(f, x)
  end
  # now just do the actual computation
  @match c begin 
    Bag          => map_ind(evalDiagr $ ObMap     , get_cat_cpt_of_flavor(ShapeVertex) ) 
    AdhesionApex => map_ind(evalDiagr $ ObMap     ,get_cat_cpt_of_flavor(ShapeEdge)   ) 
    AdhesionSpan => map_ind(map $ (evalDiagr $ HomMap), get_cat_cpt_of_flavor(ShapeSpan)   ) 
  end
end

bags(d)          = get(Bag,          d, NotIndexed)
adhesions(d)     = get(AdhesionApex, d, NotIndexed)
adhesionSpans(d) = get(AdhesionSpan, d, NotIndexed)

#=
@data GrCpt begin
  Vertex
  Edge
end

function (d::StrDecomp)(c::GrCpt, x)
  str_dcmp_cpt = @match c begin
    Vertex => Bag
    Edge   => AdhesionSpan
  end
  d_dict = get(str_dcmp_cpt, d, Indexed)
  haskey(d_dict, x) ? d_dict[x] : error(string(x) * " is not a " * string(c) )
end
=#
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
function 𝐃(f, d ::StructuredDecomposition, t::DecompType = Decomposition)::StructuredDecomposition 
  flip = @match t begin 
    Decomposition    => x -> x #nothing to do
    CoDecomposition => FinCat ∘ op_graph ∘ graph     #work with (∫G)^{op}
  end 
  δ   = d.diagram 
  X   = dom(δ)
  #Q is the composite Q = F ∘ d : FG → C → S^{op}
  Q   = FinDomFunctor(
          Dict(x => f(ob_map(δ,x))   for x ∈ ob_generators(X) ), #the ob  map Q₀
          Dict(g => f(hom_map(δ, g)) for g ∈ hom_generators(X)), #the hom map Q₁
          flip(X)
        )
  return StrDecomp(d.decomp_shape, flip(X), Q, t) 
end

end
