module Decompositions

export Decomposition, CoDecomposition, 𝐃, adhesionSpans, ∫

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

"""
Structured decompositions are diagrams.
"""
abstract type StructuredDecomposition{G, C, D} <: Diagram{id, C, D} end
export StructuredDecomposition

@data DecompType begin
  Decomposition 
  CoDecomposition
end
export DecompType

"""    Structrured decompositions

Think of these are graphs whose vertices are labeled by the objects of some category and whose edges are labeled by SPANS in this category
"""
struct StrDecomp{G, C, D} <: StructuredDecomposition{G, C, D}  
  decomp_shape ::G 
  diagram      ::D             
  decomp_type  ::DecompType
  domain       ::C
end
export StrDecomp

"""

One can construct a structured decomposition by simply providing the graph representing the shape of the decompostion and the relevant diagram. This constructor will default to setting the Decomposition Type to `Decomposition` (i.e. we default to viewing C-valued structured decompositions as diagrams into Span C)
"""
function StrDecomp(the_decomp_shape, the_diagram)
    StrDecomp(the_decomp_shape, the_diagram, Decomposition)
end

"""

If we want to explicitly specify the decomposition type of a structured decomposition, then this can be done by explicitly passing a some t::DecompType as an argument to the constructor. DecompType has two values: Decomposition and CoDecomposition. These are used to distinguish whether we think of a C-valued structured decomposition d: FG → Span C as a diagram into Span C or whether we think of it as a diagram into Cospan C of the form d: FG → Cospan C^op. 
"""
function StrDecomp(decomp_shape, diagram, decomp_type)
    d  = StrDecomp(decomp_shape, diagram, decomp_type, dom(diagram))
    dc = s -> @match decomp_type begin
        Decomposition => dom(s[1]) == dom(s[2])
        CoDecomposition => codom(s[1]) == codom(s[2])
    end
    all(dc, adhesionSpans(d)) ? d : throw(StrDecompError(d)) 
end
#construct a structured decomposition and check whether the decomposition shape actually makes sense. 
# TODO: check that the domain is also correct...

struct StrDecompError <: Exception
    d::StrDecomp
end

function Base.showerror(io::IO, e::StrDecompError)
    print(io, "$(str(d)) is not a $(str(the_decomp_type))")
end

ob_map(d::StructuredDecomposition, x)  = ob_map(d.diagram, x)
hom_map(d::StructuredDecomposition, f) = hom_map(d.diagram, f)

colimit(d::StructuredDecomposition) = colimit(FreeDiagram(d.diagram))
limit(d::StructuredDecomposition) = limit(FreeDiagram(d.diagram))

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

"""
Get the points in el::Elements corresponding to either vertices or edges
"""
function get_components(el::Elements, j::Int)
    filter(part -> el[part, :πₑ] == j, parts(el, :El))
end

"""    get_cat_components(d;:StructuredDecomposition, el::Elements, j::Int)

Get the points in el::Elements into actual objects of the category ∫G, where G is the decomposition shape.
"""
function get_cat_components(d::StructuredDecomposition, el::Elements, j::Int)
    obgen = ob_generators(d.domain)
    map(get_components(el, j)) do component
        obgen[component]
    end
end

"""    getFromDom(c::ShapeCpt, d::StructuredDecomposition, el::Elements)

For a given shape component, get the corresponding components in the category of elements 
"""
function getFromDom(c::ShapeCpt, d::StructuredDecomposition, el::Elements = elements(d.decomp_shape))
  @match c begin
    ShapeVertex => get_cat_components(d,el,1)
    ShapeEdge   => get_cat_components(d,el,2)
    ShapeSpan   => el[parts(el, :Arr) .∈ Ref(get_components(el, 2)), :src]
    # get all arrow components in the edge-objects for the given category of elements
  end
end

# TODO can we source this from an ACSet?
@data MapType begin
    ObMap
    HomMap
end

@data StrDcmpCpt begin
    Bag
    AdhesionApex
    AdhesionSpan
end

"""    map_index(f::Function, x, is_indexing::Bool=true)

If indexing, return a list of tuples (x, f.x), where x 

"""
function map_index(f::Function, x, is_indexing::Bool=true)
    is_indexing ? collect(zip(x, map(f, x))) : map(f, x)
end

"""    get(c::StrDcmpCpt, d::StructuredDecomposition, is_indexing::Bool)

Given an SD and type of component, 

`c` may be Bag, AdhesionApex, or AdhesionFoot
"""
function get(c::StrDcmpCpt, d::StructuredDecomposition, is_indexing::Bool)
    # get the category of elements from the decomposition shape
    el = elements(d.decomp_shape) 
    # cache a function for getting the component from SD and its elements
    get_cat_cpt(sc::ShapeCpt) = getFromDom(sc, d, el)
    # now just do the actual computation
    @match c begin 
        Bag          => map_index(x -> ob_map(d.diagram, x), get_cat_cpt(ShapeVertex), is_indexing) 
        AdhesionApex => map_index(x -> ob_map(d.diagram, x), get_cat_cpt(ShapeEdge), is_indexing) 
        AdhesionSpan => map_index(x -> hom_map.(Ref(d.diagram), x), get_cat_cpt(ShapeSpan), is_indexing) 
    end
end

"""    bags(d, index)
get a vector of indexed bags; i.e. a vector consisting of pairs (x, dx) where x ∈ FG and dx is the value mapped to x under the decompositon d
"""
bags(d::StructuredDecomposition, is_indexing) = get(Bag, d, is_indexing)

"""    bags(d)
get a vector of the bags of a decomposition
"""
bags(d::StructuredDecomposition) = bags(d, false)

"""
get a vector of indexed adhesions; i.e. a vector consisting of pairs (e, de) where e is an edge in ∫G and de is the value mapped to e under the decompositon d
"""
adhesions(d::StructuredDecomposition, is_indexing) = get(AdhesionApex, d, is_indexing)

"""
get a vector of the adhesions of a decomposition
"""
adhesions(d) = adhesions(d, false)     

"""  adhesionSpans(d, ind)

get a vector of indexed adhesion spans; i.e. a vector consisting of pairs (x₁ <- e -> x₂, dx₁ <- de -> dx₂) where x₁ <- e -> x₂ is span in ∫G and dx₁ <- de -> dx₂ is what its image under the decompositon d
"""
adhesionSpans(d, is_indexing) = get(AdhesionSpan, d, is_indexing)

"""  adhesionSpans(d)

get a vector of the adhesion spans of a decomposition
"""
adhesionSpans(d) = adhesionSpans(d, false)

export bags, adhesions, adhesionSpans

# TODO I'd like to develop this design idea
function FGraphTo∫(is_covariant::Bool)
    if is_covariant
        FinFunctor(Dict(:V => :El, :E => :Arr), Dict(:src => :src, :tgt => :tgt), SchGraph, SchElements)
    else
        FinFunctor(Dict(:V => :V, :E => :E), Dict(:src => :tgt, :tgt => :src), SchGraph, SchGraph)
        # TODO why?
    end
end

# TODO document this
function elements_graph(el::Elements)
    ΔF = DeltaMigration(FGraphTo∫(true))
    return migrate(Graph, el, ΔF)
end

#reverse direction of the edges
function op_graph(g::Graph)::Graph
    ΔF = DeltaMigration(FGraphTo∫(false))
    return migrate(Graph, g, ΔF)
end
# TODO upstream?

"""    ∫(G)
Syntactic sugar for constructing the category of elements of a graph. 
Note that ∫(G) has type Category whereas elements(G) has type Elements
"""
∫(G::T) where {T <: ACSet} = ∫(elements(G))
∫(G::Elements) = FinCat(elements_graph(G))

export ∫

function flip(d::StructuredDecomposition, t::DecompType)
    t == d.decomp_type ? identity : FinCat ∘ op_graph ∘ graph
end

"""    𝐃(f, d::StructuredDecomposition, t::DecompType)::StructuredDecomposition

The construction of categories of structured decompostions is functorial; 
it consists of a functor `𝐃: Cat_{pullback} → Cat` taking any category `C` with pullbacks to the category 
𝐃C of C-values structured decompositions. The functoriality of this construction allows us to lift any functor F : C → E to a functor 𝐃f : 𝐃 C → 𝐃 E which maps C-valued structured decompositions to E-valued structured decompositions. 

When we think of the functor F as a computational problem (taking inputs in C to solution spaces in E), then 𝐃f should be thought of as lifting the global comuputation of F to local computation on the constituent parts of C-valued decompositions. 

In particular, given a structured decomposition d: FG → C and a sheaf F: C → FinSet^{op} w.r.t to the decompositon topology, we can make a structured decomposition valued in FinSet^{op} by lifting the sheaf to a functor 𝐃_f: 𝐃C → 𝐃(S^{op}) between categories of structured decompositions. 
"""
function 𝐃(f, d::StructuredDecomposition, t::DecompType = d.decomp_type)::StructuredDecomposition 
  δ   = d.diagram 
  X   = dom(δ)
  ϕ = flip(d,t)
  #Q is the composite Q = F ∘ d : FG → C → S^{op}
  Q   = FinDomFunctor(
          Dict(x => f(ob_map(δ,x))   for x ∈ ob_generators(X)), #the ob  map Q₀
          Dict(g => f(hom_map(δ, g)) for g ∈ hom_generators(X)), #the hom map Q₁
          ϕ(X))
  StrDecomp(d.decomp_shape, Q, t) 
end

end
