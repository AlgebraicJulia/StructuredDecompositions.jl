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

struct Coloring <: Sheaf
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

#skeletal_coloring(n) = skeleton ∘ (Coloring(n))

𝐃_col = (𝐃 $ skeleton) ∘ (x -> 𝐃(Coloring(3), x, CoDecomposition))
#Now you can use this functor to conert a structured decomposition of graphs into a structured decomposition of the solution spaces on those graphs. 
#coloring_decomp = 𝐃_col(smallSD)
three_d = 𝐃_col(smallSD)
to_graphviz(graph(three_d.domain))
shape(three_d) == three_d.domain
diagram(three_d) == three_d.diagram
#graph(three_d.domain)
#three_d = 𝐃(skeleton, coloring_decomp)
bags(three_d)
adhesions(three_d)
as = adhesionSpans(three_d)
#trial =  map(FinFunction, as[1])

sp = as[1] 
#asd = FinFunction(sp₁[1]) 
pullback(sp[1], sp[2])

end



#=
f₁ = ACSetTransformation(H₁₂, H₁, V=[1, 3])
f₂ = ACSetTransformation(H₁₂, H₂, V=[4, 1])
fs = [f₁, f₂]
cfs = map(Coloring(3), fs)
codom(cfs[1]) == codom(cfs[2])
scfs = map(skeleton, cfs)
codom(scfs[1]) == codom(scfs[2])
dom(scfs[1])
pullback(scfs[1], scfs[2])
=#


#=
sc3 = skeletal_coloring(3)

cH₁ = collect(Coloring(3)(H₁))
q = cH₁[1]
dom(q)
f = ACSetTransformation(H₁₂, H₁, V=[1, 3])
compose(f,q)
vertices(H₁₂)
cf = Coloring(3)(f)
scf = skeleton(cf)
=#
#K = Coloring(3)
#coloring_map = K(ACSetTransformation(Γₛ⁰[3], Γₛ⁰[1], V=[1, 3]))
#all(κ -> coloring_map(κ) ∈ K(Γₛ⁰[3]), K(Γₛ⁰[1]))

#skel_map = skeleton(coloring_map)
#all(κ -> skel_map(κ) ∈ skeleton(K(Γₛ⁰[3])), skeleton(K(Γₛ⁰[1])))

#(c::Coloring)(f::ACSetTransformation) = ycol -> compose(f, ycol)
#FinFunction(c::Coloring, f::ACSetTransformation) = FinFunction(c(f), c(codom(f)))

#=
adhesionSpans(smallSD)
littlespan = adhesionSpans(smallSD)[1]
colspan = map(Coloring(3), littlespan)
codom(colspan[1]) == codom(colspan[2])
=#
#limit(colspan)

#=
c₁ = Coloring(3)(H₁) 
c₁₂ = Coloring(3)(H₁₂) 
ff = Coloring(3)(ACSetTransformation(H₁₂, H₁, V=[1, 3]))
all(x -> ff(x) ∈  c₁₂, c₁)

s₁  = FinSet(c₁)
s₁₂ = FinSet(c₁₂)
=#
#lift the sheaf to a functor between categories of sructured decompositions 

#=
K₂ = @acset Graph begin
  V = 2
  E = 1
  src = [1]
  tgt = [2]
end

spanCat = ∫(K₂) #the category x <- e -> y

function mk_span_functor(ℓ₁, ℓ₂)
  FinDomFunctor(
    Dict(1 => dom(ℓ₁), 2 => dom(ℓ₂), 3 => codom(ℓ₁)),
    Dict(1 => ℓ₁, 2 => ℓ₂),
    spanCat 
  )
end

=#






#=
pb = pullback( sp[1], sp[2] )

ob(pb)
π₁, π₂ = legs(pb)

map(g -> homomorphisms(g, K₂), bags(smallSD))
=#




#=
function myfunc(a,b) 
  x -> @match x begin
    1 => a
    2 => b
  end
end 

three = FinSet(3)
two = FinSet(2)
four = FinSet(4)
f₁ = FinFunction(myfunc(2,3), two, three)
f₂ = FinFunction(myfunc(2,4), two, four)
pushout(f₁, f₂)

s = FinSet(["a", "b", "d"])

m = FinSet([1,2])

t = FinSet(["c", "b", "d"])


ms = FinFunction(myfunc("b", "d")  , m, s)
mt = FinFunction(myfunc("b", "d") , m, t)

length(s)
dom(ms)
skeleton(ms)

pushout(skeleton(ms), skeleton(mt))

=#

#=

#ms = FinFunction(Dict(1 => "b", 2 => "d"), sl)
ms = FinFunction(α, m, s)
dom(ms) == m

codom(ms)
dom(ms) == m
mt = FinFunction(Dict(1 => "b'", 2 => "d'"))
dom(ms)
m
dom(ms) == dom(mt)


FinSet([1,2]) == FinSet([1,2])
pushout([ms, mt])

homomorphisms(H₁₂, H₁)
homomorphisms(H₁₂, H₁)
FinSet(homomorphisms(H₁₂, H₁)) == FinSet(homomorphisms(H₁₂, H₁))
#pullback(trial)
#pullback(s₁[1])
#trying out pullbacks


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
#####################