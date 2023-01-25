module TestDecidingSheaves

using Test
using PartialFunctions
using MLStyle

using ..Decompositions
using ..DecidingSheaves

using Catlab.Graphs
using Catlab.ACSetInterface
using Catlab.CategoricalAlgebra

############################
#     EXAMPLE INSTANCES
############################
"""
An example: graph colorings
"""
#an H-coloring is a hom onto H
codom_first_homs(target,source) = homomorphisms(source,target)

struct Coloring <: Sheaf
  n     #n-coloring
  func  #the function mappgin opens to lists of homs from G to K_n
end

#construct an n-coloring
Coloring(n) = Coloring(n, codom_first_homs $ (complete_graph(Graph, n)) )
#make it callable
(c::Coloring)(X::Graph) = FinSet(c.func(X))
# given graph homos f: G → H get morphism col(H) → col(G) by precomposition: take each h ∈ col(H) to hf ∈ col(G)
(c::Coloring)(f::ACSetTransformation) = FinFunction(Dict(ycol => compose(f, ycol) for ycol ∈ c(codom(f)) ) )

#(c::Coloring)(f::ACSetTransformation) = ycol -> compose(f, ycol)
#FinFunction(c::Coloring, f::ACSetTransformation) = FinFunction(c(f), c(codom(f)))

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


#lift the sheaf to a functor between categories of sructured decompositions 
𝐃_f = 𝐃 $ Coloring(3)
#Now you can use this functor to conert a structured decomposition of graphs into a structured decomposition of the solution spaces on those graphs. 
three_d = 𝐃_f(smallSD)

as = adhesionSpans(three_d)
#trial =  map(FinFunction, as[1])
sp₁ = as[1] #map(FinFunction , as[1])
#asd = FinFunction(sp₁[1]) 
dom(sp₁[1]) == dom(sp₁[2])
pushout(sp₁[1], sp₁[2])

#pullback(trial)
#pullback(s₁[1])
#trying out pullbacks

#=
s₁ = FinSet(5)
s₂ = FinSet(4)
s  = FinSet(3)
f₁ = FinFunction([1,1,2,2,3], s₁, s)
f₂ = FinFunction([2,3,1,2], s₂, s)
#σ  = Span(f₁, f₂) 
ℓ  = pullback([f₁,f₂])
=#

#the functor
#=
D₀ = Dict(1 => FinSet(3), 2 => FinSet(3), 3 => FinSet(4), 4 => FinSet(1), 5 => FinSet(2))
D₁ = Dict(
      1 => FinFunction([3], D₀[4], D₀[1]), 
      2 => FinFunction([1,2], D₀[5], D₀[2]),
      3 => FinFunction([4], D₀[4], D₀[2]),
      4 => FinFunction([1,5], D₀[5], D₀[3])
    )
D = FinDomFunctor(D₀, D₁, ∫(G))
=#
#the structured decomposition
############################# 

end