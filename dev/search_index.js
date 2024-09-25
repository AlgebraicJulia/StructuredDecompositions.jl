var documenterSearchIndex = {"docs":
[{"location":"api/#Library-Reference","page":"Library Reference","title":"Library Reference","text":"","category":"section"},{"location":"api/","page":"Library Reference","title":"Library Reference","text":"Modules = [StructuredDecompositions, StructuredDecompositions.Decompositions, StructuredDecompositions.FunctorUtils, StructuredDecompositions.DecidingSheaves]","category":"page"},{"location":"api/#StructuredDecompositions.Decompositions.StrDecomp","page":"Library Reference","title":"StructuredDecompositions.Decompositions.StrDecomp","text":"Structrured decomposition struct     – think of these are graphs whose vertices are labeled by the objects of some category      and whose edges are labeled by SPANS in this category\n\n\n\n\n\n","category":"type"},{"location":"api/#StructuredDecompositions.Decompositions.StrDecomp-Tuple{Any, Any}","page":"Library Reference","title":"StructuredDecompositions.Decompositions.StrDecomp","text":"One can construct a structured decomposition by simply providing the graph representing the shape of the decompostion and the relevant diagram. This constructor will default to setting the Decompsotion Type to Decomposition (i.e. we default to viewing C-valued structured decompositions as diagrams into Span C)\n\n\n\n\n\n","category":"method"},{"location":"api/#StructuredDecompositions.Decompositions.StructuredDecomposition","page":"Library Reference","title":"StructuredDecompositions.Decompositions.StructuredDecomposition","text":"Structured decompositions are diagrams.\n\n\n\n\n\n","category":"type"},{"location":"api/#StructuredDecompositions.Decompositions.adhesionSpans-Tuple{Any, Any}","page":"Library Reference","title":"StructuredDecompositions.Decompositions.adhesionSpans","text":"get a vector of indexed adhesion spans; i.e. a vector consisting of pairs (x₁ <- e -> x₂, dx₁ <- de -> dx₂) where x₁ <- e -> x₂ is span in ∫G and dx₁ <- de -> dx₂ is what its image under the decompositon d\n\n\n\n\n\n","category":"method"},{"location":"api/#StructuredDecompositions.Decompositions.adhesionSpans-Tuple{Any}","page":"Library Reference","title":"StructuredDecompositions.Decompositions.adhesionSpans","text":"get a vector of the adhesion spans of a decomposition\n\n\n\n\n\n","category":"method"},{"location":"api/#StructuredDecompositions.Decompositions.adhesions-Tuple{Any, Any}","page":"Library Reference","title":"StructuredDecompositions.Decompositions.adhesions","text":"get a vector of indexed adhesions; i.e. a vector consisting of pairs (e, de) where e is an edge in ∫G and de is the value mapped to e under the decompositon d\n\n\n\n\n\n","category":"method"},{"location":"api/#StructuredDecompositions.Decompositions.adhesions-Tuple{Any}","page":"Library Reference","title":"StructuredDecompositions.Decompositions.adhesions","text":"get a vector of the adhesions of a decomposition\n\n\n\n\n\n","category":"method"},{"location":"api/#StructuredDecompositions.Decompositions.bags-Tuple{Any, Any}","page":"Library Reference","title":"StructuredDecompositions.Decompositions.bags","text":"get a vector of indexed bags; i.e. a vector consisting of pairs (x, dx) where x ∈ FG and dx is the value mapped to x under the decompositon d\n\n\n\n\n\n","category":"method"},{"location":"api/#StructuredDecompositions.Decompositions.bags-Tuple{Any}","page":"Library Reference","title":"StructuredDecompositions.Decompositions.bags","text":"get a vector of the bags of a decomposition\n\n\n\n\n\n","category":"method"},{"location":"api/#StructuredDecompositions.Decompositions.∫-Tuple{T} where T<:ACSets.ACSetInterface.ACSet","page":"Library Reference","title":"StructuredDecompositions.Decompositions.∫","text":"Syntactic sugar for costrucitng the category of elements of a graph.  Note that ∫(G) has type Category whereas elements(G) has type Elements\n\n\n\n\n\n","category":"method"},{"location":"api/#StructuredDecompositions.Decompositions.𝐃","page":"Library Reference","title":"StructuredDecompositions.Decompositions.𝐃","text":"The construction of categories of structured decompostions is functorial;  it consists of a functor 𝐃: Cat{pullback} → Cat taking any category C with pullbacks to the category  𝐃C of C-values structured decompositions. The functoriality of this construction allows us to lift any functor  F : C → E to a functor 𝐃f : 𝐃 C → 𝐃 E which maps C-valued structured decompositions to E-valued structured decompositions.  When we think of the functor F as a computational problem (taking inputs in C to solution spaces in E), then 𝐃f should be thought of as lifting the global comuputation of F to local computation on the constituent parts of C-valued decompositions.  In particular, given a structured decomposition d: FG → C and a sheaf F: C → FinSet^{op} w.r.t to the decompositon topology,  we can make a structured decomposition valued in FinSet^{op} by lifting the sheaf to a functor 𝐃f: 𝐃C → 𝐃(S^{op}) between  categories of structured decompositions. \n\n\n\n\n\n","category":"function"},{"location":"api/#StructuredDecompositions.DecidingSheaves.adhesion_filter-Tuple{Tuple, StructuredDecompositions.Decompositions.StructuredDecomposition}","page":"Library Reference","title":"StructuredDecompositions.DecidingSheaves.adhesion_filter","text":"Filtering algorithm.  Note: we are assuming that we only know how to work with FinSet(Int) !\n\nINPUT: a Finset^{op}-valued structured decomposition d : FG → Span Finset^{op}            (which is expected to be in co-decomposition form;              i.e. as a diagram d : FG → Cospan Finset )        and an indexed span ( (ℓ, r), ( d(ℓ), d(r) ) ) in d            (i.e a pair consisting of span (ℓ, r) in ∫G and its image under d)\n\nOUTPUT: a structured decomposition obtained by replacing the span de in d          by the span obtained by projecting the pullback of de (i.e. taking images)\n\n\n\n\n\n","category":"method"},{"location":"api/#StructuredDecompositions.DecidingSheaves.decide_sheaf_tree_shape","page":"Library Reference","title":"StructuredDecompositions.DecidingSheaves.decide_sheaf_tree_shape","text":"Solve the decision problem encoded by a sheaf.  The algorithm is as follows:    compute on each bag (optionally, if the decomposition of the solution space                         is already known, then it can be passed as an argument),   compute composites on edges,    project back down to bags   answer (providing a witness)     \"no\" if there is an empty bag;      \"yes\" otherwise.\n\n\n\n\n\n","category":"function"},{"location":"#StructuredDecompositions.jl","page":"StructuredDecompositions.jl","title":"StructuredDecompositions.jl","text":"","category":"section"},{"location":"","page":"StructuredDecompositions.jl","title":"StructuredDecompositions.jl","text":"CurrentModule = StructuredDecompositions","category":"page"},{"location":"","page":"StructuredDecompositions.jl","title":"StructuredDecompositions.jl","text":"Structured decompositions!","category":"page"}]
}