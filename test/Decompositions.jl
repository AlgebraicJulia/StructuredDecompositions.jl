module TestDecompositions

using Test
using PartialFunctions

using ..Decompositions 
using Catlab.Graphics
using Catlab.Graphs
using Catlab.ACSetInterface
using Catlab.CategoricalAlgebra

#using Catlab.Graphi

#Define the instance#######################
#bag 1
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

#adhesion 2,3
H₂₃ = @acset Graph begin
  V = 1
end

#bag 3
H₃ = @acset Graph begin
  V = 2
  E = 1
  src = [1]
  tgt = [2]
end

# Make the decomp ###########
#The shape
G = @acset Graph begin
  V = 3
  E = 2
  src = [1, 2]
  tgt = [2, 3]
end

#to_graphviz( graph( elements(G)) )

#the functor
Γ₀ = Dict(1 => H₁, 2 => H₂, 3 => H₃, 4 => H₁₂, 5 => H₂₃)
Γ = FinDomFunctor(
  Γ₀,
  Dict(
    1 => ACSetTransformation(Γ₀[4], Γ₀[1], V=[1, 3]),
    2 => ACSetTransformation(Γ₀[5], Γ₀[2], V=[1]   ),
    3 => ACSetTransformation(Γ₀[4], Γ₀[2], V=[4, 1]),
    4 => ACSetTransformation(Γ₀[5], Γ₀[3], V=[1]   )
  ),
  ∫(G)
)
#the decomposition
bigdecomp = StrDecomp(G, ∫(G), Γ)
f = ACSetTransformation(Γ₀[4], Γ₀[1], V=[1, 3])
ob_generators(bigdecomp.domain) == ob_generators((FinCat ∘ op_graph ∘ graph)(bigdecomp.domain))
hom_generators(bigdecomp.domain) == hom_generators((FinCat ∘ op_graph ∘ graph)(bigdecomp.domain))
codom(f)
FinSet(length(vertices(dom(f)))) == dom(components(f)[1])

@test H₁ ∈ bags(bigdecomp) && H₂ ∈ bags(bigdecomp) && !(H₁₂ ∈ bags(bigdecomp))

#forgetful functor vs: Gr → Set taking G to VG
function vs(X::Graph) FinSet(length(vertices(X))) end
function vs(f::ACSetTransformation) components(f)[1] end

𝐃ᵥ = 𝐃 $ vs
bigdecomp_to_sets = 𝐃ᵥ(bigdecomp)
@test all( 
          s -> dom(s[1]) == dom(s[2]), 
          adhesionSpans(bigdecomp_to_sets)
        )

end