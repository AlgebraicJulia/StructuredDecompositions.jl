module TestDecidingSheaves

using Test
using PartialFunctions
using MLStyle

using ..Decompositions
using ..DecidingSheaves
using ..FunctorUtils

using Catlab.Graphs
using Catlab.ACSetInterface
using Catlab.CategoricalAlgebra
using Catlab.Graphics

############################
#     EXAMPLE INSTANCE str decomp
############################

H₁ = @acset Graph begin
  V = 3
  E = 2
  src = [1, 2]
  tgt = [2, 3]
end

#to_graphviz(H₁)

#adhesion 1,2
H₁₂ = @acset Graph begin
  V = 2
end

#bag 2
H₂ = @acset Graph begin
  V = 4
  E = 3
  src = [1, 2, 3]
  tgt = [2, 3, 4]
end


Gₛ = @acset Graph begin
  V = 2
  E = 1
  src = [1]
  tgt = [2]
end

Γₛ⁰ = Dict(1 => H₁, 2 => H₂, 3 => H₁₂)
Γₛ = FinDomFunctor(
  Γₛ⁰,
  Dict(
    1 => ACSetTransformation(Γₛ⁰[3], Γₛ⁰[1], V=[1, 3]),
    2 => ACSetTransformation(Γₛ⁰[3], Γₛ⁰[2], V=[4, 1]),
  ),
  ∫(Gₛ)
)
smallSD = StrDecomp(Gₛ, ∫(Gₛ), Γₛ)
as_SD   = adhesionSpans(smallSD)
length(as_SD)
colim_SD = pushout(as_SD[1])
to_graphviz(ob(colim_SD))

"""
An example: graph colorings
"""
#an H-coloring is a hom onto H
K₂ = @acset Graph begin
  V = 2
  E = 1
  src = [1]
  tgt = [2]
end

struct Coloring
  n     #the target graph
  func  #the function mappgin opens to lists of homs from G to K_n
end

#construct an n-coloring
Coloring(n) = Coloring(n, g -> homomorphisms(g, complete_graph(Graph, n)) )
#make it callable
(c::Coloring)(X::Graph) = FinSet(c.func(X)) # given graph homos #f: G₁ → G₂ get morphism col(G₂) → col(G₁) by precomposition: take each λ₂ ∈ col(G₂) to hf ∈ col(G)
function (c::Coloring)(f::ACSetTransformation)  
  (G₁, G₂)   = (dom(f), codom(f)) 
  (cG₁, cG₂) = (c(G₁), c(G₂))
  FinFunction( λ₂ -> compose(f,λ₂), cG₂, cG₁ ) #note the contravariance
end

skeletalColoring(n) = skeleton ∘ Coloring(n)
#is_homomorphic(ob(colim_SD), complete_graph(Graph, 2))

#skeletal_coloring(n) = skeleton ∘ (Coloring(n))

𝐃_col = (𝐃 $ skeleton) ∘ (x -> 𝐃(Coloring(3), x, CoDecomposition))
#Now you can use this functor to conert a structured decomposition of graphs into a structured decomposition of the solution spaces on those graphs. 
#coloring_decomp = 𝐃_col(smallSD)
three_d = 𝐃_col(smallSD)
adhesion_filter(adhesionSpans(smallSD, true)[1], three_d)
#decide_sheaf_tree_shape(skeletalColoring(2), smallSD)

end