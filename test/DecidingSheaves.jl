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

K(n)=complete_graph(Graph, n)
Gs = Dict([i => ∫(K(i)) for i in 1:7])

# we can see that
# ∫(K(1)) = *
# ∫(K(2)) = 4 -> 2 <- 3 -> 1 <- 4

############################
#     EXAMPLE INSTANCE str decomp
############################

# bag 1
H₁ = @acset Graph begin
  V = 3
  E = 2
  src = [1, 2]
  tgt = [2, 3]
end

# adhesion 1,2
H₁₂ = @acset Graph begin
  V = 2
end

# bag 2
H₂ = @acset Graph begin
  V = 4
  E = 3
  src = [1, 2, 3]
  tgt = [2, 3, 4]
end

# the shape of the decomposition
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

# build a functor from ∫G --> FinSet
Γₛ⁰ = Dict(1 => H₁, 2 => H₂, 3 => H₁₂)
Γₛ = FinDomFunctor(
  Γₛ⁰,
  Dict(
    1 => ACSetTransformation(Γₛ⁰[3], Γₛ⁰[1], V=[1, 3]),
    2 => ACSetTransformation(Γₛ⁰[3], Γₛ⁰[2], V=[4, 1]),
  ),
  ∫(Gₛ)
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

colorability_test(n, the_test_case) = is_homomorphic(ob(colimit(the_test_case)), K(n)) == decide_sheaf_tree_shape(skeletalColoring(n), the_test_case)[1]

is_homomorphic(ob(colimit(my_decomp)), K(2))

@test decide_sheaf_tree_shape(skeletalColoring(2), my_decomp)[1] == false

@test all(colorability_test(n, my_decomp) for n ∈ range(1,10))


end
