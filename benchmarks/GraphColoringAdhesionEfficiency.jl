# This file will benchmark the decrease in efficiency of the adhesion filter and deciding sheaves algorithm when the number of adhesions is increased.
# An increase in the number of adhesions brings the algorithm more in line with brute force.
# This is something that can be seen by looking at the complexity but is worth at least noting and giving an example of.

import Graphs as GraphsPkg

using BenchmarkTools
using PkgBenchmark
using Profile

using StructuredDecompositions.Decompositions
using StructuredDecompositions.DecidingSheaves
using StructuredDecompositions.FunctorUtils

using Catlab.Graphs
using Catlab.ACSetInterface
using Catlab.CategoricalAlgebra
using Catlab.Graphics

# The following structures and functions were pulled from DecidingSheaves.jl 

struct Coloring
    n     
    func
end

K(n) = complete_graph(Graph, n)
Coloring(n) = Coloring(n, g -> homomorphisms(g, K(n)))
(c::Coloring)(X::Graph) = FinSet(c.func(X))
function (c::Coloring)(f::ACSetTransformation)  
    (G₁, G₂)   = (dom(f), codom(f)) 
    (cG₁, cG₂) = (c(G₁), c(G₂))
    FinFunction( λ₂ -> compose(f,λ₂), cG₂, cG₁ )
end
  
skeletalColoring(n) = skeleton ∘ Coloring(n)
  
function colorability_test(n, the_test_case)
  hom = is_homomorphic(ob(colimit(the_test_case)), K(n))
  dec = decide_sheaf_tree_shape(skeletalColoring(n), the_test_case)[1]
  if hom == dec
    return hom
  else
    error("is_homomorphic != decide_sheaf_tree_shape")
  end
end

# my_decomp1 is a graph of 4 diamonds composed of two triangles. The diamonds are in a line, connected at end points
# my_decomp2 is similarly 4 diamonds each composed of 2 triangles. Each diamond is connected to the next at exactly 2 points opposite of where it is connected to the other.


# my_decomp1 is a graph of 4 diamonds composed of two triangles. The diamonds are in a line, connected at end points
#bag 1
H₁ = @acset Graph begin
  V = 4
  E = 5
  src = [1, 1, 2, 2, 3]
  tgt = [2, 3, 3, 4, 4]
end
  
#adhesion 1,2
H₁₂ = @acset Graph begin
  V = 1
end
  
#bag 2
H₂ = @acset Graph begin
  V = 4
  E = 5
  src = [1, 1, 2, 2, 3]
  tgt = [2, 3, 3, 4, 4]
end
  
#adhesion 3,4
H₂₃ = @acset Graph begin
  V = 1
end

#bag 3
H₃ = @acset Graph begin
  V = 4
  E = 5
  src = [1, 1, 2, 2, 3]
  tgt = [2, 3, 3, 4, 4]
end

#adhesion 3,4
H₃₄ = @acset Graph begin
  V = 1
end

#bag 4
H₄ = @acset Graph begin
  V = 4
  E = 5
  src = [1, 1, 2, 2, 3]
  tgt = [2, 3, 3, 4, 4]
end

Gₛ = @acset Graph begin
  V = 4
  E = 3
  src = [1, 1, 1]
  tgt = [2, 3, 4]
end

#transformations
Γₛ⁰ = Dict(1 => H₁, 2 => H₂, 3 => H₃, 4 => H₄, 5 => H₁₂, 6 => H₂₃, 7 => H₃₄)
Γₛ = FinDomFunctor(
  Γₛ⁰,
  Dict(
    1 => ACSetTransformation(Γₛ⁰[5], Γₛ⁰[1], V=[4]),
    2 => ACSetTransformation(Γₛ⁰[5], Γₛ⁰[2], V=[1]),
    3 => ACSetTransformation(Γₛ⁰[6], Γₛ⁰[2], V=[4]),
    4 => ACSetTransformation(Γₛ⁰[6], Γₛ⁰[3], V=[1]),
    5 => ACSetTransformation(Γₛ⁰[7], Γₛ⁰[3], V=[4]),
    6 => ACSetTransformation(Γₛ⁰[7], Γₛ⁰[4], V=[1])
  ),
  ∫(Gₛ)
)
  
my_decomp1  = StrDecomp(Gₛ, Γₛ)

# my_decomp2 is similarly 4 diamonds each composed of 2 triangles. Each diamond is connected to the next at exactly 2 points opposite of where it is connected to the other.
#bag 1
H₁ = @acset Graph begin
  V = 4
  E = 5
  src = [1, 1, 2, 2, 3]
  tgt = [2, 3, 3, 4, 4]
end
  
#adhesion 1,2
H₁₂ = @acset Graph begin
  V = 2
end
  
#bag 2
H₂ = @acset Graph begin
  V = 4
  E = 5
  src = [1, 1, 2, 2, 3]
  tgt = [2, 3, 3, 4, 4]
end
    
#adhesion 3,4
H₂₃ = @acset Graph begin
  V = 2
end
  
#bag 3
H₃ = @acset Graph begin
  V = 4
  E = 5
  src = [1, 1, 2, 2, 3]
  tgt = [2, 3, 3, 4, 4]
end
  
#adhesion 3,4
H₃₄ = @acset Graph begin
  V = 2
end
  
#bag 4
H₄ = @acset Graph begin
  V = 4
  E = 5
  src = [1, 1, 2, 2, 3]
  tgt = [2, 3, 3, 4, 4]
end

Gₛ = @acset Graph begin
  V = 4
  E = 3
  src = [1, 2, 3]
  tgt = [2, 3, 4]
end

Γₛ⁰ = Dict(1 => H₁, 2 => H₂, 3 => H₃, 4 => H₄, 5 => H₁₂, 6 => H₂₃, 7 => H₃₄)
Γₛ = FinDomFunctor(
  Γₛ⁰,
  Dict(
    1 => ACSetTransformation(Γₛ⁰[5], Γₛ⁰[1], V=[3, 4]),
    2 => ACSetTransformation(Γₛ⁰[5], Γₛ⁰[2], V=[1, 2]),
    3 => ACSetTransformation(Γₛ⁰[6], Γₛ⁰[2], V=[3, 4]),
    4 => ACSetTransformation(Γₛ⁰[6], Γₛ⁰[3], V=[1, 2]),
    5 => ACSetTransformation(Γₛ⁰[7], Γₛ⁰[3], V=[3, 4]),
    6 => ACSetTransformation(Γₛ⁰[7], Γₛ⁰[4], V=[1, 2])
  ),
  ∫(Gₛ)
)

my_decomp2 = StrDecomp(Gₛ, Γₛ)

decide_sheaf_tree_shape(skeletalColoring(2), my_decomp1)[1] == false
@benchmark decide_sheaf_tree_shape(skeletalColoring(2), my_decomp1)[1]
  
decide_sheaf_tree_shape(skeletalColoring(3), my_decomp1)[1] == true
@benchmark decide_sheaf_tree_shape(skeletalColoring(3), my_decomp1)[1]
  
decide_sheaf_tree_shape(skeletalColoring(4), my_decomp1)[1] == true
@benchmark decide_sheaf_tree_shape(skeletalColoring(4), my_decomp1)[1]

decide_sheaf_tree_shape(skeletalColoring(2), my_decomp2)[1] == false
@benchmark decide_sheaf_tree_shape(skeletalColoring(2), my_decomp2)[1]
  
decide_sheaf_tree_shape(skeletalColoring(3), my_decomp2)[1] == true
@benchmark decide_sheaf_tree_shape(skeletalColoring(3), my_decomp2)[1]
  
decide_sheaf_tree_shape(skeletalColoring(4), my_decomp2)[1] == true
@benchmark decide_sheaf_tree_shape(skeletalColoring(4), my_decomp2)[1]