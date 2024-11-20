var documenterSearchIndex = {"docs":
[{"location":"api/#Library-Reference","page":"Library Reference","title":"Library Reference","text":"","category":"section"},{"location":"api/","page":"Library Reference","title":"Library Reference","text":"Modules = [StructuredDecompositions, StructuredDecompositions.Decompositions, StructuredDecompositions.FunctorUtils, StructuredDecompositions.DecidingSheaves, StructuredDecompositions.JunctionTrees, StructuredDecompositions.NestedUWDs]","category":"page"},{"location":"api/#StructuredDecompositions.Decompositions.StrDecomp","page":"Library Reference","title":"StructuredDecompositions.Decompositions.StrDecomp","text":"Structured decomposition struct     – Consider these as graphs whose vertices are labeled by the objects of some category      and whose edges are labeled by spans in this category\n\n\n\n\n\n","category":"type"},{"location":"api/#StructuredDecompositions.Decompositions.StrDecomp-Tuple{Any, Any}","page":"Library Reference","title":"StructuredDecompositions.Decompositions.StrDecomp","text":"One can construct a structured decomposition by simply providing the graph representing the shape of the decompostion and the relevant diagram. This constructor will default to setting the Decompsotion Type to Decomposition (i.e. we default to viewing C-valued structured decompositions as diagrams into Span C)\n\n\n\n\n\n","category":"method"},{"location":"api/#StructuredDecompositions.Decompositions.StructuredDecomposition","page":"Library Reference","title":"StructuredDecompositions.Decompositions.StructuredDecomposition","text":"Structured decompositions are diagrams.\n\n\n\n\n\n","category":"type"},{"location":"api/#StructuredDecompositions.Decompositions.adhesionSpans-Tuple{Any, Any}","page":"Library Reference","title":"StructuredDecompositions.Decompositions.adhesionSpans","text":"Get a vector of indexed adhesion spans; i.e. a vector consisting of pairs (x₁ <- e -> x₂, dx₁ <- de -> dx₂) where x₁ <- e -> x₂ is span in ∫G and dx₁ <- de -> dx₂ is what its image under the decompositon d\n\n\n\n\n\n","category":"method"},{"location":"api/#StructuredDecompositions.Decompositions.adhesionSpans-Tuple{Any}","page":"Library Reference","title":"StructuredDecompositions.Decompositions.adhesionSpans","text":"Get a vector of the adhesion spans of a decomposition\n\n\n\n\n\n","category":"method"},{"location":"api/#StructuredDecompositions.Decompositions.adhesions-Tuple{Any, Any}","page":"Library Reference","title":"StructuredDecompositions.Decompositions.adhesions","text":"Get a vector of indexed adhesions; i.e. a vector consisting of pairs (e, de) where e is an edge in ∫G and de is the value mapped to e under the decompositon d\n\n\n\n\n\n","category":"method"},{"location":"api/#StructuredDecompositions.Decompositions.adhesions-Tuple{Any}","page":"Library Reference","title":"StructuredDecompositions.Decompositions.adhesions","text":"Get a vector of the adhesions of a decomposition\n\n\n\n\n\n","category":"method"},{"location":"api/#StructuredDecompositions.Decompositions.bags-Tuple{Any, Any}","page":"Library Reference","title":"StructuredDecompositions.Decompositions.bags","text":"Get a vector of indexed bags; i.e. a vector consisting of pairs (x, dx) where x ∈ FG and dx is the value mapped to x under the decompositon d\n\n\n\n\n\n","category":"method"},{"location":"api/#StructuredDecompositions.Decompositions.bags-Tuple{Any}","page":"Library Reference","title":"StructuredDecompositions.Decompositions.bags","text":"Get a vector of the bags of a decomposition\n\n\n\n\n\n","category":"method"},{"location":"api/#StructuredDecompositions.Decompositions.∫-Tuple{T} where T<:ACSets.ACSetInterface.ACSet","page":"Library Reference","title":"StructuredDecompositions.Decompositions.∫","text":"Syntactic sugar for costructing the category of elements of a graph.  Note that ∫(G) has type Category whereas elements(G) has type Elements\n\n\n\n\n\n","category":"method"},{"location":"api/#StructuredDecompositions.Decompositions.𝐃","page":"Library Reference","title":"StructuredDecompositions.Decompositions.𝐃","text":"The construction of categories of structured decompostions is functorial;  it consists of a functor 𝐃: Cat{pullback} → Cat taking any category C with pullbacks to the category  𝐃C of C-values structured decompositions. The functoriality of this construction allows us to lift any functor  F : C → E to a functor 𝐃f : 𝐃 C → 𝐃 E which maps C-valued structured decompositions to E-valued structured decompositions.  When we think of the functor F as a computational problem (taking inputs in C to solution spaces in E), then 𝐃f should be thought of as lifting the global comuputation of F to local computation on the constituent parts of C-valued decompositions.  In particular, given a structured decomposition d: FG → C and a sheaf F: C → FinSet^{op} w.r.t to the decompositon topology,  we can make a structured decomposition valued in FinSet^{op} by lifting the sheaf to a functor 𝐃f: 𝐃C → 𝐃(S^{op}) between  categories of structured decompositions. \n\n\n\n\n\n","category":"function"},{"location":"api/#StructuredDecompositions.DecidingSheaves.decide_sheaf_tree_shape","page":"Library Reference","title":"StructuredDecompositions.DecidingSheaves.decide_sheaf_tree_shape","text":"Solve the decision problem encoded by a sheaf.  The algorithm is as follows:    compute on each bag    (optionally, if the decomposition of the solution space   is already known, then it can be passed as an argument),   compute composites on edges,    project back down to bags   answer (providing a witness)     \"no\" if there is an empty bag;      \"yes\" otherwise.\n\n\n\n\n\n","category":"function"},{"location":"api/#StructuredDecompositions.DecidingSheaves.old_adhesion_filter-Tuple{Tuple, StructuredDecompositions.Decompositions.StructuredDecomposition}","page":"Library Reference","title":"StructuredDecompositions.DecidingSheaves.old_adhesion_filter","text":"Filtering algorithm.  Note: we are assuming that we only know how to work with FinSet(Int) !\n\nINPUT: a Finset^{op}-valued structured decomposition d : FG → Span Finset^{op}            (which is expected to be in co-decomposition form;            i.e. as a diagram d : FG → Cospan Finset )           and an indexed span ( (ℓ, r), ( d(ℓ), d(r) ) ) in d            (i.e a pair consisting of span (ℓ, r) in ∫G and its image under d)\n\nOUTPUT: a structured decomposition obtained by replacing the span de in d          by the span obtained by projecting the pullback of de (i.e. taking images)\n\n\n\n\n\n","category":"method"},{"location":"api/#StructuredDecompositions.DecidingSheaves.old_decide_sheaf_tree_shape","page":"Library Reference","title":"StructuredDecompositions.DecidingSheaves.old_decide_sheaf_tree_shape","text":"Solve the decision problem encoded by a sheaf.  The algorithm is as follows:    compute on each bag (optionally, if the decomposition of the solution space                         is already known, then it can be passed as an argument),   compute composites on edges,    project back down to bags   answer (providing a witness)     \"no\" if there is an empty bag;      \"yes\" otherwise.\n\n\n\n\n\n","category":"function"},{"location":"api/#StructuredDecompositions.JunctionTrees.AMDJL_AMD","page":"Library Reference","title":"StructuredDecompositions.JunctionTrees.AMDJL_AMD","text":"AMDJL_AMD <: EliminationAlgorithm\n\nThe approximate minimum degree algorithm. Uses AMD.jl.\n\n\n\n\n\n","category":"type"},{"location":"api/#StructuredDecompositions.JunctionTrees.CuthillMcKeeJL_RCM","page":"Library Reference","title":"StructuredDecompositions.JunctionTrees.CuthillMcKeeJL_RCM","text":"CuthillMcKeeJL_RCM <: EliminationAlgorithm\n\nThe reverse Cuthill-McKee algorithm. Uses CuthillMckee.jl.\n\n\n\n\n\n","category":"type"},{"location":"api/#StructuredDecompositions.JunctionTrees.EliminationAlgorithm","page":"Library Reference","title":"StructuredDecompositions.JunctionTrees.EliminationAlgorithm","text":"EliminationAlgorithm\n\nA graph elimination algorithm. The options are\n\nCuthillMcKeeJL_RCM\nAMDJL_AMD\nMetisJL_ND\nMCS\n\n\n\n\n\n","category":"type"},{"location":"api/#StructuredDecompositions.JunctionTrees.FundamentalSupernode","page":"Library Reference","title":"StructuredDecompositions.JunctionTrees.FundamentalSupernode","text":"FundamentalSupernode <: Supernode\n\nA fundamental supernode.\n\n\n\n\n\n","category":"type"},{"location":"api/#StructuredDecompositions.JunctionTrees.MCS","page":"Library Reference","title":"StructuredDecompositions.JunctionTrees.MCS","text":"MCS <: EliminationAlgorithm\n\nThe maximum cardinality search algorithm.\n\n\n\n\n\n","category":"type"},{"location":"api/#StructuredDecompositions.JunctionTrees.MaximalSupernode","page":"Library Reference","title":"StructuredDecompositions.JunctionTrees.MaximalSupernode","text":"MaximalSupernode <: Supernode\n\nA maximal supernode.\n\n\n\n\n\n","category":"type"},{"location":"api/#StructuredDecompositions.JunctionTrees.MetisJL_ND","page":"Library Reference","title":"StructuredDecompositions.JunctionTrees.MetisJL_ND","text":"MetisJL_ND <: EliminationAlgorithm\n\nThe nested dissection heuristic. Uses Metis.jl.\n\n\n\n\n\n","category":"type"},{"location":"api/#StructuredDecompositions.JunctionTrees.Node","page":"Library Reference","title":"StructuredDecompositions.JunctionTrees.Node","text":"Node <: Supernode\n\nA single-vertex supernode.\n\n\n\n\n\n","category":"type"},{"location":"api/#StructuredDecompositions.JunctionTrees.Supernode","page":"Library Reference","title":"StructuredDecompositions.JunctionTrees.Supernode","text":"Supernode\n\nA type of supernode. The options are\n\nNode\nMaximalSupernode\nFundamentalSupernode\n\n\n\n\n\n","category":"type"},{"location":"api/#StructuredDecompositions.NestedUWDs.NestedUWD","page":"Library Reference","title":"StructuredDecompositions.NestedUWDs.NestedUWD","text":"NestedUWD(\n    diagram::UndirectedWiringDiagram,\n    [, algorithm::Union{Order, EliminationAlgorithm}]\n    [, supernode::Supernode])\n\nConstruct a nested undirected wiring diagram.\n\n\n\n\n\n","category":"type"},{"location":"api/#StructuredDecompositions.NestedUWDs.NestedUWD-2","page":"Library Reference","title":"StructuredDecompositions.NestedUWDs.NestedUWD","text":"NestedUWD{T, B, V}\n\nAn undirected wiring diagram, represented as a nested collected of undirected wiring diagrams.\n\n\n\n\n\n","category":"type"},{"location":"api/#StructuredDecompositions.NestedUWDs.evalschedule-Union{Tuple{T}, Tuple{Any, StructuredDecompositions.NestedUWDs.NestedUWD, AbstractVector{T}}, Tuple{Any, StructuredDecompositions.NestedUWDs.NestedUWD, AbstractVector{T}, AbstractVector}} where T","page":"Library Reference","title":"StructuredDecompositions.NestedUWDs.evalschedule","text":"function evalschedule(\n    f,\n    nuwd::NestedUWD,\n    generators::Union{AbstractVector, AbstractDict}\n    [, operations::AbstractVector])\n\nEvaluate an undirected wiring diagrams given a set of generators for the boxes. The optional first argument f should be callable with the signature\n\n    f(diagram, generators)\n\nwhere diagram is an undirected wiring diagram, and generators is a vector. If f is not specified, then it defaults to oapply.\n\n\n\n\n\n","category":"method"},{"location":"api/#StructuredDecompositions.NestedUWDs.makeschedule-Union{Tuple{StructuredDecompositions.NestedUWDs.NestedUWD{<:Any, T}}, Tuple{T}} where T","page":"Library Reference","title":"StructuredDecompositions.NestedUWDs.makeschedule","text":"makeschedule(nuwd::NestedUWD)\n\nConstruct a directed wiring diagram that represents the nesting structure of a nested UWD.\n\n\n\n\n\n","category":"method"},{"location":"pages/decompositions/#Decompositions","page":"Decompositions","title":"Decompositions","text":"","category":"section"},{"location":"pages/decompositions/","page":"Decompositions","title":"Decompositions","text":"Structured decompositions are diagrams. They can be thought of as graphs whose vertices are labeled by the objects of some category and whose edges are labeld by spans in this category. ","category":"page"},{"location":"pages/decompositions/","page":"Decompositions","title":"Decompositions","text":"abstract type StructuredDecomposition{G, C, D} <: Diagram{id, C, D} end\n\nstruct StrDecomp{G, C, D} <: StructuredDecomposition{G, C, D}  \n  decomp_shape ::G \n  diagram      ::D             \n  decomp_type  ::DecompType\n  domain       ::C  \nend","category":"page"},{"location":"pages/decompositions/","page":"Decompositions","title":"Decompositions","text":"A structured decomposition can be constructed by providing the graph representing the shap of the decomposition and the relevant diagram. The default constructor will set the decomposition type to Decomposition (viewing C-valued structured decompositions as diagrams into Span C).","category":"page"},{"location":"pages/decompositions/","page":"Decompositions","title":"Decompositions","text":"StrDecomp(the_decomp_shape, the_diagram)=StrDecomp(the_decomp_shape, the_diagram, Decomposition)","category":"page"},{"location":"pages/decompositions/","page":"Decompositions","title":"Decompositions","text":"The function StrDecomp will construct a structured decomposition and check whether the decomposition shape makes sense.","category":"page"},{"location":"pages/decompositions/","page":"Decompositions","title":"Decompositions","text":"StrDecomp(the_decomp_shape, the_diagram, the_decomp_type)","category":"page"},{"location":"pages/decompositions/","page":"Decompositions","title":"Decompositions","text":"The functions colimit and limit when called on a structured decomposition will take the colimit and limit of the diagram, respectively.","category":"page"},{"location":"pages/decompositions/","page":"Decompositions","title":"Decompositions","text":"function colimit(d::StructuredDecomposition) colimit(FreeDiagram(d.diagram)) end\nfunction limit(d::StructuredDecomposition) limit(FreeDiagram(d.diagram)) end","category":"page"},{"location":"pages/decompositions/","page":"Decompositions","title":"Decompositions","text":"The construction of categories of structured decompositions consists of a functor 𝐃:Cat_{pullback} → Cat taking any category C to the category 𝐃C of C-valued structured decompositions. This allows the lifting of any functor F: C → E to a functor 𝐃f : 𝐃C → 𝐃E which maps C-valued structured decompositions to E-valued structured decompostions.","category":"page"},{"location":"pages/decompositions/","page":"Decompositions","title":"Decompositions","text":"function 𝐃(f, d ::StructuredDecomposition, t::DecompType = d.decomp_type)::StructuredDecomposition","category":"page"},{"location":"pages/junction_trees/#JunctionTrees","page":"Junction Trees","title":"Junction Trees","text":"","category":"section"},{"location":"pages/nested_uwds/#NestedUWDs","page":"Nested UWDs","title":"Nested UWDs","text":"","category":"section"},{"location":"pages/functor_utils/#FunctorUtils","page":"Functor Utils","title":"Functor Utils","text":"","category":"section"},{"location":"pages/functor_utils/","page":"Functor Utils","title":"Functor Utils","text":"Functor Utils only includes 4 functions and builds off of Decompositions.","category":"page"},{"location":"pages/functor_utils/","page":"Functor Utils","title":"Functor Utils","text":"We first define the forgetful functors vs.","category":"page"},{"location":"pages/functor_utils/","page":"Functor Utils","title":"Functor Utils","text":"function vs(X::Graph) FinSet(length(vertices(X))) end\nfunction vs(f::ACSetTransformation) components(f)[1] end","category":"page"},{"location":"pages/functor_utils/","page":"Functor Utils","title":"Functor Utils","text":"We also define the functor skeleton taking set to the skeleton of the set.","category":"page"},{"location":"pages/functor_utils/","page":"Functor Utils","title":"Functor Utils","text":"function skeleton(s::FinSet) FinSet(length(s)) end\nfunction skeleton(f::FinFunction)\n  (dd, cc) = (dom(f), codom(f))\n  ℓ = isempty(dd) ? Int[] : [findfirst(item -> item == f(x), collect(cc)) for x ∈ collect(dd)]\n  FinFunction(ℓ, skeleton(dd), skeleton(cc))\nend","category":"page"},{"location":"pages/decidingsheaves/#DecidingSheaves","page":"Deciding Sheaves","title":"Deciding Sheaves","text":"","category":"section"},{"location":"pages/decidingsheaves/","page":"Deciding Sheaves","title":"Deciding Sheaves","text":"There are two functions that are used here. ","category":"page"},{"location":"pages/decidingsheaves/","page":"Deciding Sheaves","title":"Deciding Sheaves","text":"The first one called adhesion_filter will take an input Finset^{op}-valued structured decomposition and return a new structured decompostion replacing the span de in d by the span obtained by projecting the pullback of de. ","category":"page"},{"location":"pages/decidingsheaves/","page":"Deciding Sheaves","title":"Deciding Sheaves","text":"function adhesion_filter(tup::Tuple, d::StructuredDecomposition)","category":"page"},{"location":"pages/decidingsheaves/","page":"Deciding Sheaves","title":"Deciding Sheaves","text":"The second function is called decidesheaftree_shape and solves a decision problem that is encoded by a sheaf. This function works by computing first on each of the bags, then computes composites on edges and projects back down to bags.","category":"page"},{"location":"pages/decidingsheaves/","page":"Deciding Sheaves","title":"Deciding Sheaves","text":"function decide_sheaf_tree_shape(f, d::StructuredDecomposition, solution_space_decomp::StructuredDecomposition = 𝐃(f, d, CoDecomposition))","category":"page"},{"location":"pages/decidingsheaves/#Graph-Coloring-Structure-and-Examples","page":"Deciding Sheaves","title":"Graph Coloring Structure and Examples","text":"","category":"section"},{"location":"pages/decidingsheaves/","page":"Deciding Sheaves","title":"Deciding Sheaves","text":"One of the many decision problems that can be solved using this system of functions is figuring out graph colorings. To put it simply, a graph coloring is an assignment of labels(colors) to each vertex of a graph so that no two adjacent vertices have the same label(color).","category":"page"},{"location":"pages/decidingsheaves/","page":"Deciding Sheaves","title":"Deciding Sheaves","text":"We first define the structure of a coloring.","category":"page"},{"location":"pages/decidingsheaves/","page":"Deciding Sheaves","title":"Deciding Sheaves","text":"#an H-coloring is a hom onto H\nstruct Coloring\n  n     #the target graph\n  func  #the function mapping opens to lists of homs from G to K_n\nend","category":"page"},{"location":"pages/decidingsheaves/","page":"Deciding Sheaves","title":"Deciding Sheaves","text":"We then define how we are going to create and test colorings.","category":"page"},{"location":"pages/decidingsheaves/","page":"Deciding Sheaves","title":"Deciding Sheaves","text":"#construct an n-coloring\nK(n) = complete_graph(Graph, n)\nColoring(n) = Coloring(n, g -> homomorphisms(g, K(n) ))\n#make it callable\n(c::Coloring)(X::Graph) = FinSet(c.func(X)) # given graph homos #f: G₁ → G₂ get morphism col(G₂) → col(G₁) by precomposition: take each λ₂ ∈ col(G₂) to hf ∈ col(G)\nfunction (c::Coloring)(f::ACSetTransformation)  \n  (G₁, G₂)   = (dom(f), codom(f)) \n  (cG₁, cG₂) = (c(G₁), c(G₂))\n  FinFunction( λ₂ -> compose(f,λ₂), cG₂, cG₁ ) #note the contravariance\nend\n\nskeletalColoring(n) = skeleton ∘ Coloring(n)\n\ncolorability_test(n, the_test_case) = is_homomorphic(ob(colimit(the_test_case)), K(n)) == decide_sheaf_tree_shape(skeletalColoring(n), the_test_case)[1]","category":"page"},{"location":"pages/decidingsheaves/","page":"Deciding Sheaves","title":"Deciding Sheaves","text":"We can consider the example of a ring with seven nodes as our graph.  We first seperate the nodes into bags with adhesions and what our adhesions look like.","category":"page"},{"location":"pages/decidingsheaves/","page":"Deciding Sheaves","title":"Deciding Sheaves","text":"#bag 1\nH₁ = @acset Graph begin\n    V = 3\n    E = 2\n    src = [1, 2]\n    tgt = [2, 3]\nend\n\n#adhesion 1,2\nH₁₂ = @acset Graph begin\n    V = 2\nend\n\n#bag 2\nH₂ = @acset Graph begin\n    V = 4\n    E = 3\n    src = [1, 2, 3]\n    tgt = [2, 3, 4]\nend\n\nGₛ = @acset Graph begin\n    V = 2\n    E = 1\n    src = [1]\n    tgt = [2]\nend","category":"page"},{"location":"pages/decidingsheaves/","page":"Deciding Sheaves","title":"Deciding Sheaves","text":"We can then construct our structured decomposition through transformations of these bags and adhesions.","category":"page"},{"location":"pages/decidingsheaves/","page":"Deciding Sheaves","title":"Deciding Sheaves","text":"#transformations\nΓₛ⁰ = Dict(1 => H₁, 2 => H₂, 3 => H₁₂)\nΓₛ = FinDomFunctor(\n    Γₛ⁰,\n    Dict(\n      1 => ACSetTransformation(Γₛ⁰[3], Γₛ⁰[1], V=[1, 3]),\n      2 => ACSetTransformation(Γₛ⁰[3], Γₛ⁰[2], V=[4, 1]),\n    ),\n    ∫(Gₛ)\n)\n\nmy_decomp1  = StrDecomp(Gₛ, Γₛ)","category":"page"},{"location":"#StructuredDecompositions.jl","page":"StructuredDecompositions.jl","title":"StructuredDecompositions.jl","text":"","category":"section"},{"location":"","page":"StructuredDecompositions.jl","title":"StructuredDecompositions.jl","text":"CurrentModule = StructuredDecompositions","category":"page"},{"location":"","page":"StructuredDecompositions.jl","title":"StructuredDecompositions.jl","text":"Structural graph theorists and algorithmicists alike know that it's usually a smart idea to decompose graphs into smaller and simpler parts before trying answer difficult computational questions. Tree decompositions are one of the best-known ways of systematically chopping graphs up and, as such they have been key tools for establishing deep results in many areas of discrete mathematics including graph minor theory and algorithmic meta-theorems. ","category":"page"},{"location":"#Generality","page":"StructuredDecompositions.jl","title":"Generality","text":"","category":"section"},{"location":"","page":"StructuredDecompositions.jl","title":"StructuredDecompositions.jl","text":"This project consists of an implementation of the category-theoretic notion of structured decompositions. These provide a formalism for decomposing arbitrary mathematical objects (not just graphs) and morphisms between them into smaller constituent parts. Since the definition of a structured decompositions is functorial, one can easily lift computational problems (defined as functors mapping inputs to solution spaces) to functors between from decompositions of the inputs to decompositions of solution spaces.","category":"page"},{"location":"#What-is-in-this-package?","page":"StructuredDecompositions.jl","title":"What is in this package?","text":"","category":"section"},{"location":"","page":"StructuredDecompositions.jl","title":"StructuredDecompositions.jl","text":"This package allows one to leverage insights to solve decision problems that are encoded as sheaves efficiently (i.e. in fixed-parameter-tractable time parameterized by the width of the decompositions). Currently, this packages includes many general tools that can be used to decompose arbitrary mathematical objects. One of the many applications of this package, exampled in the Deciding Sheaves module, is graph colorings.","category":"page"},{"location":"","page":"StructuredDecompositions.jl","title":"StructuredDecompositions.jl","text":"Pages = [\n    \"pages/decompostions.md\",\n    \"pages/decidingsheaves.md\",\n    \"pages/junction_trees.md\",\n    \"pages/nested_uwds.md\",\n    \"pages/functor_utils.md\",\n    ]\nDepth = 2","category":"page"}]
}
