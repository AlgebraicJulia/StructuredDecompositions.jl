# This file examines the relationship between run time, bag size, and number of bags for a fixed 40 node graph.
# We will compare the results against running the graph coloring algorithm on Graphs.jl

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

# We will start with a small number of vertices in a small number of bags and consider the case below
# Adding vertices to the bags at a constant rate


# We first consider the 2 bag, 5 vertices graph

#bag 1
H₁ = @acset Graph begin
  V = 3
  E = 2
  src = [1, 2]
  tgt = [2, 3]
end

#adhesion 1,2
H₁₂ = @acset Graph begin
  V = 1
end

#bag 2
H₂ = @acset Graph begin
  V = 3
  E = 2
  src = [1, 2]
  tgt = [2, 3]
end

Gₛ = @acset Graph begin
  V = 2
  E = 1
  src = [1]
  tgt = [1]
end

#transformations
Γₛ⁰ = Dict(1 => H₁, 2 => H₂, 3 => H₁₂)
Γₛ = FinDomFunctor(
  Γₛ⁰,
  Dict(
   1 => ACSetTransformation(Γₛ⁰[3], Γₛ⁰[1], V=[2]),
   2 => ACSetTransformation(Γₛ⁰[3], Γₛ⁰[2], V=[2]),
 ),
 ∫(Gₛ)
)
  
my_decomp1  = StrDecomp(Gₛ, Γₛ)

# We now consider the 2 bag, 8 vertices graph

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

Gₛ = @acset Graph begin
  V = 2
  E = 1
  src = [1]
  tgt = [1]
end

#transformations
Γₛ⁰ = Dict(1 => H₁, 2 => H₂, 3 => H₁₂)
Γₛ = FinDomFunctor(
  Γₛ⁰,
  Dict(
   1 => ACSetTransformation(Γₛ⁰[3], Γₛ⁰[1], V=[1]),
   2 => ACSetTransformation(Γₛ⁰[3], Γₛ⁰[2], V=[1]),
 ),
 ∫(Gₛ)
)
  
my_decomp2  = StrDecomp(Gₛ, Γₛ)

# We now consider the 2 bag, 12 vertices graph

#bag 1
H₁ = @acset Graph begin
  V = 6
  E = 8
  src = [1, 1, 2, 2, 3, 4, 4, 5]
  tgt = [2, 3, 3, 4, 4, 5, 6, 6]
end

#adhesion 1,2
H₁₂ = @acset Graph begin
  V = 1
end

#bag 2
H₂ = @acset Graph begin
  V = 6
  E = 8
  src = [1, 1, 2, 2, 3, 4, 4, 5]
  tgt = [2, 3, 3, 4, 4, 5, 6, 6]
end

Gₛ = @acset Graph begin
  V = 2
  E = 1
  src = [1]
  tgt = [1]
end

#transformations
Γₛ⁰ = Dict(1 => H₁, 2 => H₂, 3 => H₁₂)
Γₛ = FinDomFunctor(
  Γₛ⁰,
  Dict(
   1 => ACSetTransformation(Γₛ⁰[3], Γₛ⁰[1], V=[1]),
   2 => ACSetTransformation(Γₛ⁰[3], Γₛ⁰[2], V=[1]),
 ),
 ∫(Gₛ)
)
  
my_decomp3  = StrDecomp(Gₛ, Γₛ)

# We now consider the 2 bag, 16 vertices graph

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

Gₛ = @acset Graph begin
  V = 2
  E = 1
  src = [1]
  tgt = [1]
end
  
#transformations
Γₛ⁰ = Dict(1 => H₁, 2 => H₂, 3 => H₁₂)
Γₛ = FinDomFunctor(
  Γₛ⁰,
  Dict(
   1 => ACSetTransformation(Γₛ⁰[3], Γₛ⁰[1], V=[1]),
   2 => ACSetTransformation(Γₛ⁰[3], Γₛ⁰[2], V=[1]),
  ),
 ∫(Gₛ)
)
  
my_decomp2  = StrDecomp(Gₛ, Γₛ)

# We now consider the 2 bag, 12 vertices graph

#bag 1
H₁ = @acset Graph begin
  V = 8
  E = 10
  src = [1, 1, 2, 2, 3, 4, 4, 5, 5, 6]
  tgt = [2, 3, 3, 4, 4, 5, 6, 6, 7, 7]
end
  
#adhesion 1,2
H₁₂ = @acset Graph begin
  V = 1
end
  
#bag 2
H₂ = @acset Graph begin
  V = 8
  E = 10
  src = [1, 1, 2, 2, 3, 4, 4, 5, 5, 6]
  tgt = [2, 3, 3, 4, 4, 5, 6, 6, 7, 7]
end
  
Gₛ = @acset Graph begin
  V = 2
  E = 1
  src = [1]
  tgt = [1]
end
  
#transformations
Γₛ⁰ = Dict(1 => H₁, 2 => H₂, 3 => H₁₂)
Γₛ = FinDomFunctor(
  Γₛ⁰,
  Dict(
   1 => ACSetTransformation(Γₛ⁰[3], Γₛ⁰[1], V=[1]),
   2 => ACSetTransformation(Γₛ⁰[3], Γₛ⁰[2], V=[1]),
  ),
  ∫(Gₛ)
)
  
my_decomp4  = StrDecomp(Gₛ, Γₛ)

decide_sheaf_tree_shape(skeletalColoring(2), my_decomp1)[1] == true
decide_sheaf_tree_shape(skeletalColoring(2), my_decomp2)[1] == false
decide_sheaf_tree_shape(skeletalColoring(2), my_decomp3)[1] == false
decide_sheaf_tree_shape(skeletalColoring(2), my_decomp4)[1] == false
@benchmark decide_sheaf_tree_shape(skeletalColoring(2), my_decomp1)[1]
@benchmark decide_sheaf_tree_shape(skeletalColoring(2), my_decomp2)[1]
@benchmark decide_sheaf_tree_shape(skeletalColoring(2), my_decomp3)[1]
@benchmark decide_sheaf_tree_shape(skeletalColoring(2), my_decomp4)[1]
@benchmark colorability_test(2, my_decomp1)
@benchmark colorability_test(2, my_decomp2)
@benchmark colorability_test(2, my_decomp3)
@benchmark colorability_test(2, my_decomp4)

decide_sheaf_tree_shape(skeletalColoring(3), my_decomp1)[1] == true
decide_sheaf_tree_shape(skeletalColoring(3), my_decomp2)[1] == true
decide_sheaf_tree_shape(skeletalColoring(3), my_decomp3)[1] == true
decide_sheaf_tree_shape(skeletalColoring(3), my_decomp4)[1] == true
@benchmark decide_sheaf_tree_shape(skeletalColoring(3), my_decomp1)[1]
@benchmark decide_sheaf_tree_shape(skeletalColoring(3), my_decomp2)[1]
@benchmark decide_sheaf_tree_shape(skeletalColoring(3), my_decomp3)[1]
@benchmark decide_sheaf_tree_shape(skeletalColoring(3), my_decomp4)[1]
@benchmark colorability_test(3, my_decomp1)
@benchmark colorability_test(3, my_decomp2)
@benchmark colorability_test(3, my_decomp3)
@benchmark colorability_test(3, my_decomp4)

decide_sheaf_tree_shape(skeletalColoring(4), my_decomp1)[1] == true
decide_sheaf_tree_shape(skeletalColoring(4), my_decomp2)[1] == true
decide_sheaf_tree_shape(skeletalColoring(4), my_decomp3)[1] == true
decide_sheaf_tree_shape(skeletalColoring(4), my_decomp4)[1] == true
@benchmark decide_sheaf_tree_shape(skeletalColoring(4), my_decomp1)[1]
@benchmark decide_sheaf_tree_shape(skeletalColoring(4), my_decomp2)[1]
@benchmark decide_sheaf_tree_shape(skeletalColoring(4), my_decomp3)[1]
@benchmark decide_sheaf_tree_shape(skeletalColoring(4), my_decomp4)[1]
@benchmark colorability_test(4, my_decomp1)
@benchmark colorability_test(4, my_decomp2)
@benchmark colorability_test(4, my_decomp3)
@benchmark colorability_test(4, my_decomp4)