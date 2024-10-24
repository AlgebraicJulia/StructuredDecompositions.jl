module TestDecidingSheaves

using Test
using PartialFunctions
using MLStyle

using StructuredDecompositions.Decompositions
using StructuredDecompositions.DecidingSheaves
using StructuredDecompositions.FunctorUtils

using Catlab.Graphs
using Catlab.ACSetInterface
using Catlab.CategoricalAlgebra
using Catlab.Graphics

############################
#     EXAMPLE INSTANCE str decomp
############################

# Bag 1
# graph: 1 -> 2 -> 3
H₁ = @acset Graph begin
  V = 3
  E = 2
  src = [1, 2]
  tgt = [2, 3]
end
# Adhesion 1,2
# graph: 1 2
H₁₂ = @acset Graph begin
  V = 2
end
# Bag 2
# graph: 1 -> 2 -> 3 -> 4
H₂ = @acset Graph begin
  V = 4
  E = 3
  src = [1, 2, 3]
  tgt = [2, 3, 4]
end
# Decomposition Shape
# graph: 1 -> 2
Gₛ = @acset Graph begin
  V = 2
  E = 1
  src = [1]
  tgt = [2]
end
# ∫(Gₛ) produces a finitely-presented category
#  1:2 ⇉ 1:3
#  -- it accepts an ACSet
#  -- produces its Elements
#  -- produces an elements_graph

decomp_shape = ∫(Gₛ)

#= build a functor from ∫G --> FinSet
The first bag goes to the ACSetTransformation(Adhesion to First Bag). The first vertex of the adhesion goes to the first vertex of the bag and the second vertex goes to the third vertex of the bag.
The second bag goes to the ACSetTransformation into the Second Bag. The first vertex of the adhesion goes to the fourth vertex of the bag and the second vertex goes to the first vertex of the bag.
Quotienting the resulting decomposition returns a five-vertex cycle graph.
=#
Γₛ⁰ = Dict(1 => H₁, 2 => H₂, 3 => H₁₂)
Γₛ = FinDomFunctor(
  Γₛ⁰,
  Dict(
    1 => ACSetTransformation(Γₛ⁰[3], Γₛ⁰[1], V=[1, 3]),
    2 => ACSetTransformation(Γₛ⁰[3], Γₛ⁰[2], V=[4, 1]),
  ),
  decomp_shape
)

my_decomp   = StrDecomp(Gₛ, Γₛ)

"""
An example: graph colorings
"""
#an H-coloring is a hom onto H
struct Coloring
  n     #the target graph
  func  #the function mappgin opens to lists of homs from G to K_n
end

#construct an n-coloring
K(n) = complete_graph(Graph, n)
Coloring(n) = Coloring(n, g -> homomorphisms(g, K(n) ))
#make it callable
(c::Coloring)(X::Graph) = FinSet(c.func(X)) # given graph homos #f: G₁ → G₂ get morphism col(G₂) → col(G₁) by precomposition: take each λ₂ ∈ col(G₂) to hf ∈ col(G)
function (c::Coloring)(f::ACSetTransformation)  
  (G₁, G₂)   = (dom(f), codom(f)) 
  (cG₁, cG₂) = (c(G₁), c(G₂))
  FinFunction( λ₂ -> compose(f,λ₂), cG₂, cG₁ ) #note the contravariance
end

skeletalColoring(n) = skeleton ∘ Coloring(n)

function colorability_test(n, the_test_case)
    is_homomorphic(ob(colimit(the_test_case)), K(n)) == decide_sheaf_tree_shape(skeletalColoring(n), the_test_case)[1]
end

is_homomorphic(ob(colimit(my_decomp)), K(2))

@test decide_sheaf_tree_shape(skeletalColoring(2), my_decomp)[1] == false

@test all(colorability_test(n, my_decomp) for n ∈ range(1,10))

# test functionality

using StructuredDecompositions.Decompositions: get_components, get_cat_components

elH1 = elements(H₁)
elH2 = elements(H₂)
elH12 = elements(H₁₂)

@test get_components(elH1, 1) == [1,2,3] # there are three vertex-objects
@test get_components(elH1, 2) == [4,5] # there are two edge-objects
@test get_components(elH2, 1) == [1,2,3,4] # there are four vertex-objects
@test get_components(elH2, 2) == [5,6,7] # there are three edge-objects
@test get_components(elH12, 1) == [1,2] # there are two vertex-objects
@test get_components(elH12, 2) == [] # there are no edge-objects

# @test get_cat_components(ssd, elH1, 1) == [1,2,3]
# @test_broken get_cat_components(ssd, elH2, 1) == [1,2]
# @test get_cat_components(ssd, elH12, 1) == [1,2]
# @test get_cat_components(ssd, elH12, 2) == []

# using StructuredDecompositions.Decompositions: ShapeVertex, ShapeEdge, ShapeSpan, getFromDom

# @test getFromDom(ShapeVertex, ssd, elH1) == [1,2,3]
# @test_broken getFromDom(ShapeEdge, ssd, elH1) == [4]
# @test getFromDom(ShapeSpan, ssd, elH1) == [5]

# @test bags(ssd) == [FinSet(2), FinSet(2)]
# @test adhesions(ssd) == [FinSet(4)]

# @test get(Bag, ssd, false) == [FinSet(2), FinSet(2)]
# @test get(AdhesionApex, ssd, false) == [FinSet(4)]

end
