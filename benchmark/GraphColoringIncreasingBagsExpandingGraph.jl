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
# Adding more bags with the same number of vertices at a constant rate
# For the purpose of this benchmark we will only consider bags of size 4 and size 8

# 2 bags size 4

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
  
bags2size4  = StrDecomp(Gₛ, Γₛ)

# 3 bags size 4

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

#adhesion 2,3
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

Gₛ = @acset Graph begin
  V = 3
  E = 2
  src = [1, 1]
  tgt = [2, 3]
end
    
#transformations
Γₛ⁰ = Dict(1 => H₁, 2 => H₂, 3 => H₃, 4 => H₁₂, 5 => H₂₃)
Γₛ = FinDomFunctor(
  Γₛ⁰,
  Dict(
   1 => ACSetTransformation(Γₛ⁰[4], Γₛ⁰[1], V=[1]),
   2 => ACSetTransformation(Γₛ⁰[4], Γₛ⁰[2], V=[1]),
   3 => ACSetTransformation(Γₛ⁰[5], Γₛ⁰[2], V=[1]),
   4 => ACSetTransformation(Γₛ⁰[5], Γₛ⁰[3], V=[1]),
 ),
 ∫(Gₛ)
)
  
bags3size4  = StrDecomp(Gₛ, Γₛ)

# 4 bags size 4

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

#adhesion 2,3
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
   1 => ACSetTransformation(Γₛ⁰[5], Γₛ⁰[1], V=[1]),
   2 => ACSetTransformation(Γₛ⁰[5], Γₛ⁰[2], V=[1]),
   3 => ACSetTransformation(Γₛ⁰[6], Γₛ⁰[2], V=[1]),
   4 => ACSetTransformation(Γₛ⁰[6], Γₛ⁰[3], V=[1]),
   5 => ACSetTransformation(Γₛ⁰[7], Γₛ⁰[3], V=[1]),
   6 => ACSetTransformation(Γₛ⁰[7], Γₛ⁰[4], V=[1])
 ),
 ∫(Gₛ)
)
  
bags4size4  = StrDecomp(Gₛ, Γₛ)

# 8 bags size 4

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

#adhesion 2,3
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

#adhesion 4,5
H₄₅ = @acset Graph begin
  V = 1
end

#bag 5
H₅ = @acset Graph begin
  V = 4
  E = 5
  src = [1, 1, 2, 2, 3]
  tgt = [2, 3, 3, 4, 4]
end
  
#adhesion 5,6
H₅₆ = @acset Graph begin
  V = 1
end

#bag 6
H₆ = @acset Graph begin
  V = 4
  E = 5
  src = [1, 1, 2, 2, 3]
  tgt = [2, 3, 3, 4, 4]
end

#adhesion 6,7
H₆₇ = @acset Graph begin
  V = 1
end

#bag 7
H₇ = @acset Graph begin
  V = 4
  E = 5
  src = [1, 1, 2, 2, 3]
  tgt = [2, 3, 3, 4, 4]
end
  
#adhesion 7,8
H₇₈ = @acset Graph begin
  V = 1
end
  
#bag 8
H₈ = @acset Graph begin
  V = 4
  E = 5
  src = [1, 1, 2, 2, 3]
  tgt = [2, 3, 3, 4, 4]
end


Gₛ = @acset Graph begin
  V = 8
  E = 7
  src = [1, 1, 1, 1, 1, 1, 1]
  tgt = [2, 3, 4, 5, 6, 7, 8]
end
    
#transformations
Γₛ⁰ = Dict(1 => H₁, 2 => H₂, 3 => H₃, 4 => H₄, 5 => H₅, 6 => H₆, 7 => H₇, 8 => H₈, 9 => H₁₂, 10 => H₂₃, 11 => H₃₄,
          12 => H₄₅, 13 => H₅₆, 14 => H₆₇, 15 => H₇₈)
Γₛ = FinDomFunctor(
  Γₛ⁰,
  Dict(
   1 => ACSetTransformation(Γₛ⁰[9], Γₛ⁰[1], V=[1]),
   2 => ACSetTransformation(Γₛ⁰[9], Γₛ⁰[2], V=[1]),
   3 => ACSetTransformation(Γₛ⁰[10], Γₛ⁰[2], V=[1]),
   4 => ACSetTransformation(Γₛ⁰[10], Γₛ⁰[3], V=[1]),
   5 => ACSetTransformation(Γₛ⁰[11], Γₛ⁰[3], V=[1]),
   6 => ACSetTransformation(Γₛ⁰[11], Γₛ⁰[4], V=[1]),
   7 => ACSetTransformation(Γₛ⁰[12], Γₛ⁰[4], V=[1]),
   8 => ACSetTransformation(Γₛ⁰[12], Γₛ⁰[5], V=[1]),
   9 => ACSetTransformation(Γₛ⁰[13], Γₛ⁰[5], V=[1]),
   10 => ACSetTransformation(Γₛ⁰[13], Γₛ⁰[6], V=[1]),
   11 => ACSetTransformation(Γₛ⁰[14], Γₛ⁰[6], V=[1]),
   12 => ACSetTransformation(Γₛ⁰[14], Γₛ⁰[7], V=[1]),
   13 => ACSetTransformation(Γₛ⁰[15], Γₛ⁰[7], V=[1]),
   14 => ACSetTransformation(Γₛ⁰[15], Γₛ⁰[8], V=[1])
 ),
 ∫(Gₛ)
)
    
bags8size4  = StrDecomp(Gₛ, Γₛ)

# benchmarks on size 4 bags
decide_sheaf_tree_shape(skeletalColoring(2), bags2size4)[1] == false
decide_sheaf_tree_shape(skeletalColoring(2), bags3size4)[1] == false
decide_sheaf_tree_shape(skeletalColoring(2), bags4size4)[1] == false
decide_sheaf_tree_shape(skeletalColoring(2), bags8size4)[1] == false
@benchmark decide_sheaf_tree_shape(skeletalColoring(2), bags2size4)[1]
@benchmark decide_sheaf_tree_shape(skeletalColoring(2), bags3size4)[1]
@benchmark decide_sheaf_tree_shape(skeletalColoring(2), bags4size4)[1]
@benchmark decide_sheaf_tree_shape(skeletalColoring(2), bags8size4)[1]
@benchmark colorability_test(2, bags2size4)
@benchmark colorability_test(2, bags3size4)
@benchmark colorability_test(2, bags4size4)
@benchmark colorability_test(2, bags8size4)

decide_sheaf_tree_shape(skeletalColoring(3), bags2size4)[1] == true
decide_sheaf_tree_shape(skeletalColoring(3), bags3size4)[1] == true
decide_sheaf_tree_shape(skeletalColoring(3), bags4size4)[1] == true
decide_sheaf_tree_shape(skeletalColoring(3), bags8size4)[1] == true
@benchmark decide_sheaf_tree_shape(skeletalColoring(3), bags2size4)[1]
@benchmark decide_sheaf_tree_shape(skeletalColoring(3), bags3size4)[1]
@benchmark decide_sheaf_tree_shape(skeletalColoring(3), bags4size4)[1]
@benchmark decide_sheaf_tree_shape(skeletalColoring(3), bags8size4)[1]
@benchmark colorability_test(3, bags2size4)
@benchmark colorability_test(3, bags3size4)
@benchmark colorability_test(3, bags4size4)
@benchmark colorability_test(3, bags8size4)

decide_sheaf_tree_shape(skeletalColoring(4), bags2size4)[1] == true
decide_sheaf_tree_shape(skeletalColoring(4), bags3size4)[1] == true
decide_sheaf_tree_shape(skeletalColoring(4), bags4size4)[1] == true
decide_sheaf_tree_shape(skeletalColoring(4), bags8size4)[1] == true
@benchmark decide_sheaf_tree_shape(skeletalColoring(4), bags2size4)[1]
@benchmark decide_sheaf_tree_shape(skeletalColoring(4), bags3size4)[1]
@benchmark decide_sheaf_tree_shape(skeletalColoring(4), bags4size4)[1]
@benchmark decide_sheaf_tree_shape(skeletalColoring(4), bags8size4)[1]
@benchmark colorability_test(4, bags2size4)
@benchmark colorability_test(4, bags3size4)
@benchmark colorability_test(4, bags4size4)
@benchmark colorability_test(4, bags8size4)

# 2 bags size 8

#bag 1
H₁ = @acset Graph begin
  V = 8
  E = 11
  src = [1, 1, 2, 2, 3, 4, 4, 5, 5, 6, 7]
  tgt = [2, 3, 3, 4, 4, 5, 6, 6, 7, 7, 8]
end
  
#adhesion 1,2
H₁₂ = @acset Graph begin
  V = 1
end

#bag 2
H₂ = @acset Graph begin
  V = 8
  E = 11
  src = [1, 1, 2, 2, 3, 4, 4, 5, 5, 6, 7]
  tgt = [2, 3, 3, 4, 4, 5, 6, 6, 7, 7, 8]
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
  
bags2size8  = StrDecomp(Gₛ, Γₛ)

# 3 bags size 8

#bag 1
H₁ = @acset Graph begin
  V = 8
  E = 11
  src = [1, 1, 2, 2, 3, 4, 4, 5, 5, 6, 7]
  tgt = [2, 3, 3, 4, 4, 5, 6, 6, 7, 7, 8]
end
  
#adhesion 1,2
H₁₂ = @acset Graph begin
  V = 1
end

#bag 2
H₂ = @acset Graph begin
  V = 8
  E = 11
  src = [1, 1, 2, 2, 3, 4, 4, 5, 5, 6, 7]
  tgt = [2, 3, 3, 4, 4, 5, 6, 6, 7, 7, 8]
end

#adhesion 2,3
H₂₃ = @acset Graph begin
  V = 1
end

#bag 3
H₃ = @acset Graph begin
  V = 8
  E = 11
  src = [1, 1, 2, 2, 3, 4, 4, 5, 5, 6, 7]
  tgt = [2, 3, 3, 4, 4, 5, 6, 6, 7, 7, 8]
end

Gₛ = @acset Graph begin
  V = 3
  E = 2
  src = [1, 1]
  tgt = [2, 3]
end
    
#transformations
Γₛ⁰ = Dict(1 => H₁, 2 => H₂, 3 => H₃, 4 => H₁₂, 5 => H₂₃)
Γₛ = FinDomFunctor(
  Γₛ⁰,
  Dict(
   1 => ACSetTransformation(Γₛ⁰[4], Γₛ⁰[1], V=[1]),
   2 => ACSetTransformation(Γₛ⁰[4], Γₛ⁰[2], V=[1]),
   3 => ACSetTransformation(Γₛ⁰[5], Γₛ⁰[2], V=[1]),
   4 => ACSetTransformation(Γₛ⁰[5], Γₛ⁰[3], V=[1]),
 ),
 ∫(Gₛ)
)
  
bags3size8  = StrDecomp(Gₛ, Γₛ)

# 4 bags size 8

#bag 1
H₁ = @acset Graph begin
  V = 8
  E = 11
  src = [1, 1, 2, 2, 3, 4, 4, 5, 5, 6, 7]
  tgt = [2, 3, 3, 4, 4, 5, 6, 6, 7, 7, 8]
end
  
#adhesion 1,2
H₁₂ = @acset Graph begin
  V = 1
end

#bag 2
H₂ = @acset Graph begin
  V = 8
  E = 11
  src = [1, 1, 2, 2, 3, 4, 4, 5, 5, 6, 7]
  tgt = [2, 3, 3, 4, 4, 5, 6, 6, 7, 7, 8]
end

#adhesion 2,3
H₂₃ = @acset Graph begin
  V = 1
end

#bag 3
H₃ = @acset Graph begin
  V = 8
  E = 11
  src = [1, 1, 2, 2, 3, 4, 4, 5, 5, 6, 7]
  tgt = [2, 3, 3, 4, 4, 5, 6, 6, 7, 7, 8]
end

#adhesion 3,4
H₃₄ = @acset Graph begin
  V = 1
end

#bag 4
H₄ = @acset Graph begin
  V = 8
  E = 11
  src = [1, 1, 2, 2, 3, 4, 4, 5, 5, 6, 7]
  tgt = [2, 3, 3, 4, 4, 5, 6, 6, 7, 7, 8]
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
   1 => ACSetTransformation(Γₛ⁰[5], Γₛ⁰[1], V=[1]),
   2 => ACSetTransformation(Γₛ⁰[5], Γₛ⁰[2], V=[1]),
   3 => ACSetTransformation(Γₛ⁰[6], Γₛ⁰[2], V=[1]),
   4 => ACSetTransformation(Γₛ⁰[6], Γₛ⁰[3], V=[1]),
   5 => ACSetTransformation(Γₛ⁰[7], Γₛ⁰[3], V=[1]),
   6 => ACSetTransformation(Γₛ⁰[7], Γₛ⁰[4], V=[1])
 ),
 ∫(Gₛ)
)
  
bags4size8  = StrDecomp(Gₛ, Γₛ)

# 8 bags size 8

#bag 1
H₁ = @acset Graph begin
  V = 8
  E = 11
  src = [1, 1, 2, 2, 3, 4, 4, 5, 5, 6, 7]
  tgt = [2, 3, 3, 4, 4, 5, 6, 6, 7, 7, 8]
end
  
#adhesion 1,2
H₁₂ = @acset Graph begin
  V = 1
end

#bag 2
H₂ = @acset Graph begin
  V = 8
  E = 11
  src = [1, 1, 2, 2, 3, 4, 4, 5, 5, 6, 7]
  tgt = [2, 3, 3, 4, 4, 5, 6, 6, 7, 7, 8]
end

#adhesion 2,3
H₂₃ = @acset Graph begin
  V = 1
end

#bag 3
H₃ = @acset Graph begin
  V = 8
  E = 11
  src = [1, 1, 2, 2, 3, 4, 4, 5, 5, 6, 7]
  tgt = [2, 3, 3, 4, 4, 5, 6, 6, 7, 7, 8]
end

#adhesion 3,4
H₃₄ = @acset Graph begin
  V = 1
end

#bag 4
H₄ = @acset Graph begin
  V = 8
  E = 11
  src = [1, 1, 2, 2, 3, 4, 4, 5, 5, 6, 7]
  tgt = [2, 3, 3, 4, 4, 5, 6, 6, 7, 7, 8]
end

#adhesion 4,5
H₄₅ = @acset Graph begin
  V = 1
end

#bag 5
H₅ = @acset Graph begin
  V = 8
  E = 11
  src = [1, 1, 2, 2, 3, 4, 4, 5, 5, 6, 7]
  tgt = [2, 3, 3, 4, 4, 5, 6, 6, 7, 7, 8]
end
  
#adhesion 5,6
H₅₆ = @acset Graph begin
  V = 1
end

#bag 6
H₆ = @acset Graph begin
  V = 8
  E = 11
  src = [1, 1, 2, 2, 3, 4, 4, 5, 5, 6, 7]
  tgt = [2, 3, 3, 4, 4, 5, 6, 6, 7, 7, 8]
end

#adhesion 6,7
H₆₇ = @acset Graph begin
  V = 1
end

#bag 7
H₇ = @acset Graph begin
  V = 8
  E = 11
  src = [1, 1, 2, 2, 3, 4, 4, 5, 5, 6, 7]
  tgt = [2, 3, 3, 4, 4, 5, 6, 6, 7, 7, 8]
end
  
#adhesion 7,8
H₇₈ = @acset Graph begin
  V = 1
end
  
#bag 8
H₈ = @acset Graph begin
  V = 8
  E = 11
  src = [1, 1, 2, 2, 3, 4, 4, 5, 5, 6, 7]
  tgt = [2, 3, 3, 4, 4, 5, 6, 6, 7, 7, 8]
end
 
Gₛ = @acset Graph begin
  V = 8
  E = 7
  src = [1, 1, 1, 1, 1, 1, 1]
  tgt = [2, 3, 4, 5, 6, 7, 8]
end
    
#transformations
Γₛ⁰ = Dict(1 => H₁, 2 => H₂, 3 => H₃, 4 => H₄, 5 => H₅, 6 => H₆, 7 => H₇, 8 => H₈, 9 => H₁₂, 10 => H₂₃, 11 => H₃₄,
          12 => H₄₅, 13 => H₅₆, 14 => H₆₇, 15 => H₇₈)
Γₛ = FinDomFunctor(
  Γₛ⁰,
  Dict(
   1 => ACSetTransformation(Γₛ⁰[9], Γₛ⁰[1], V=[1]),
   2 => ACSetTransformation(Γₛ⁰[9], Γₛ⁰[2], V=[1]),
   3 => ACSetTransformation(Γₛ⁰[10], Γₛ⁰[2], V=[1]),
   4 => ACSetTransformation(Γₛ⁰[10], Γₛ⁰[3], V=[1]),
   5 => ACSetTransformation(Γₛ⁰[11], Γₛ⁰[3], V=[1]),
   6 => ACSetTransformation(Γₛ⁰[11], Γₛ⁰[4], V=[1]),
   7 => ACSetTransformation(Γₛ⁰[12], Γₛ⁰[4], V=[1]),
   8 => ACSetTransformation(Γₛ⁰[12], Γₛ⁰[5], V=[1]),
   9 => ACSetTransformation(Γₛ⁰[13], Γₛ⁰[5], V=[1]),
   10 => ACSetTransformation(Γₛ⁰[13], Γₛ⁰[6], V=[1]),
   11 => ACSetTransformation(Γₛ⁰[14], Γₛ⁰[6], V=[1]),
   12 => ACSetTransformation(Γₛ⁰[14], Γₛ⁰[7], V=[1]),
   13 => ACSetTransformation(Γₛ⁰[15], Γₛ⁰[7], V=[1]),
   14 => ACSetTransformation(Γₛ⁰[15], Γₛ⁰[8], V=[1])
 ),
 ∫(Gₛ)
)
    
bags8size8  = StrDecomp(Gₛ, Γₛ)

# benchmarks on size 4 bags
decide_sheaf_tree_shape(skeletalColoring(2), bags2size8)[1] == false
decide_sheaf_tree_shape(skeletalColoring(2), bags3size8)[1] == false
decide_sheaf_tree_shape(skeletalColoring(2), bags4size8)[1] == false
decide_sheaf_tree_shape(skeletalColoring(2), bags8size8)[1] == false
@benchmark decide_sheaf_tree_shape(skeletalColoring(2), bags2size8)[1]
@benchmark decide_sheaf_tree_shape(skeletalColoring(2), bags3size8)[1]
@benchmark decide_sheaf_tree_shape(skeletalColoring(2), bags4size8)[1]
@benchmark decide_sheaf_tree_shape(skeletalColoring(2), bags8size8)[1]
@benchmark colorability_test(2, bags2size8)
@benchmark colorability_test(2, bags3size8)
@benchmark colorability_test(2, bags4size8)
@benchmark colorability_test(2, bags8size8)

decide_sheaf_tree_shape(skeletalColoring(3), bags2size8)[1] == true
decide_sheaf_tree_shape(skeletalColoring(3), bags3size8)[1] == true
decide_sheaf_tree_shape(skeletalColoring(3), bags4size8)[1] == true
decide_sheaf_tree_shape(skeletalColoring(3), bags8size8)[1] == true
@benchmark decide_sheaf_tree_shape(skeletalColoring(3), bags2size8)[1]
@benchmark decide_sheaf_tree_shape(skeletalColoring(3), bags3size8)[1]
@benchmark decide_sheaf_tree_shape(skeletalColoring(3), bags4size8)[1]
@benchmark decide_sheaf_tree_shape(skeletalColoring(3), bags8size8)[1]
@benchmark colorability_test(3, bags2size8)
@benchmark colorability_test(3, bags3size8)
@benchmark colorability_test(3, bags4size8)
@benchmark colorability_test(3, bags8size8)

decide_sheaf_tree_shape(skeletalColoring(4), bags2size8)[1] == true
decide_sheaf_tree_shape(skeletalColoring(4), bags3size8)[1] == true
decide_sheaf_tree_shape(skeletalColoring(4), bags4size8)[1] == true
decide_sheaf_tree_shape(skeletalColoring(4), bags8size8)[1] == true
@benchmark decide_sheaf_tree_shape(skeletalColoring(4), bags2size8)[1]
@benchmark decide_sheaf_tree_shape(skeletalColoring(4), bags3size8)[1]
@benchmark decide_sheaf_tree_shape(skeletalColoring(4), bags4size8)[1]
@benchmark decide_sheaf_tree_shape(skeletalColoring(4), bags8size8)[1]
@benchmark colorability_test(4, bags2size8)
@benchmark colorability_test(4, bags3size8)
@benchmark colorability_test(4, bags4size8)
@benchmark colorability_test(4, bags8size8)