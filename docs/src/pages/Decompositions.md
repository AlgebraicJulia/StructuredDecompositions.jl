# [Decompositions] (@id Decompositions)

Structured decompositions are diagrams. They can be thought of as graphs whose vertices are labeled by the objects of some category and whose edges are labeld by spans in this category. 

```julia
abstract type StructuredDecomposition{G, C, D} <: Diagram{id, C, D} end

struct StrDecomp{G, C, D} <: StructuredDecomposition{G, C, D}  
  decomp_shape ::G 
  diagram      ::D             
  decomp_type  ::DecompType
  domain       ::C  
end
```

A structured decomposition can be constructed by providing the graph representing the shap of the decomposition and the relevant diagram. The default constructor will set the decomposition type to Decomposition (viewing C-valued structured decompositions as diagrams into Span C).

```julia
StrDecomp(the_decomp_shape, the_diagram)=StrDecomp(the_decomp_shape, the_diagram, Decomposition)
```

The function StrDecomp will construct a structured decomposition and check whether the decomposition shape makes sense.

```julia
StrDecomp(the_decomp_shape, the_diagram, the_decomp_type)
```

The functions colimit and limit when called on a structured decomposition will take the colimit and limit of the diagram, respectively.

```julia
function colimit(d::StructuredDecomposition) colimit(FreeDiagram(d.diagram)) end
function limit(d::StructuredDecomposition) limit(FreeDiagram(d.diagram)) end
```

The construction of categories of structured decompositions consists of a functor 𝐃:Cat_{pullback} → Cat taking any category C to the category 𝐃C of C-valued structured decompositions. This allows the lifting of any functor F: C → E to a functor 𝐃f : 𝐃C → 𝐃E which maps C-valued structured decompositions to E-valued structured decompostions.

```julia
function 𝐃(f, d ::StructuredDecomposition, t::DecompType = d.decomp_type)::StructuredDecomposition
```
