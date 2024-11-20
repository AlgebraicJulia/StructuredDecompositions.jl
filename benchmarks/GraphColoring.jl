# This file benchmarks the StructuredDecompositions method of graph coloring

# RESULTS = sheaf tree skeletal coloring is n linear for false returns, worse than n exponential for true returns
    # shouldn't be worse than exponential for true returns so I will take a look at it after considering k

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

# Benchmark 1 (small n per bag(4) and small bags k(2), 3 coloring)

#bag 1
H₁ = @acset Graph begin
  V = 4
  E = 4
  src = [1, 2, 3, 4]
  tgt = [2, 3, 4, 1]
end

#adhesion 1,2
H₁₂ = @acset Graph begin
  V = 2
end

#bag 2
H₂ = @acset Graph begin
  V = 4
  E = 5
  src = [1, 2, 3, 4, 1]
  tgt = [2, 3, 4, 1, 3]
end

Gₛ = @acset Graph begin
  V = 2
  E = 1
  src = [1]
  tgt = [2]
end

#transformations
Γₛ⁰ = Dict(1 => H₁, 2 => H₂, 3 => H₁₂)
Γₛ = FinDomFunctor(
  Γₛ⁰,
  Dict(
   1 => ACSetTransformation(Γₛ⁰[3], Γₛ⁰[1], V=[3, 4]),
   2 => ACSetTransformation(Γₛ⁰[3], Γₛ⁰[2], V=[2, 1]),
 ),
 ∫(Gₛ)
)

my_decomp1  = StrDecomp(Gₛ, Γₛ)

is_homomorphic(ob(colimit(my_decomp1)), K(2)) == false
@benchmark is_homomorphic(ob(colimit(my_decomp1)), K(2))

is_homomorphic(ob(colimit(my_decomp1)), K(3)) == true
@benchmark is_homomorphic(ob(colimit(my_decomp1)), K(3))

is_homomorphic(ob(colimit(my_decomp1)), K(4)) == true
@benchmark is_homomorphic(ob(colimit(my_decomp1)), K(4))

decide_sheaf_tree_shape(skeletalColoring(2), my_decomp1)[1] == false
@benchmark decide_sheaf_tree_shape(skeletalColoring(2), my_decomp1)[1]
@time decide_sheaf_tree_shape(skeletalColoring(2), my_decomp1)[1]

decide_sheaf_tree_shape(skeletalColoring(3), my_decomp1)[1] == true
@benchmark decide_sheaf_tree_shape(skeletalColoring(3), my_decomp1)[1]

decide_sheaf_tree_shape(skeletalColoring(4), my_decomp1)[1] == true
@benchmark decide_sheaf_tree_shape(skeletalColoring(4), my_decomp1)[1]

colorability_test(2, my_decomp1) == false
@benchmark colorability_test(2, my_decomp1)

colorability_test(3, my_decomp1) == true
@benchmark colorability_test(3, my_decomp1)

colorability_test(4, my_decomp1) == true
@benchmark colorability_test(4, my_decomp1)

# Benchmark 2 (medium n per bag(10) and small bags k(2), 3 coloring)

#bag 1
H₁ = @acset Graph begin
  V = 10
  E = 16
  src = [1, 1, 1, 2, 2, 3, 4, 4, 4, 6, 6, 7, 7, 7, 8, 9]
  tgt = [2, 5, 6, 3, 8, 4, 5, 9, 10, 7, 10, 8, 9, 10, 9, 10]
end

#adhesion 1,2
H₁₂ = @acset Graph begin
  V = 2
end

#bag 2
H₂ = @acset Graph begin
  V = 10
  E = 18
  src = [1, 1, 2, 2, 2, 2, 3, 4, 4, 5, 5, 6, 6, 7, 7, 7, 8, 9]
  tgt = [2, 9, 3, 4, 9, 10, 4, 5, 10, 6, 10, 7, 10, 8, 9, 10, 9, 10]
end

Gₛ = @acset Graph begin
  V = 2
  E = 1
  src = [1]
  tgt = [2]
end

#transformations
Γₛ⁰ = Dict(1 => H₁, 2 => H₂, 3 => H₁₂)
Γₛ = FinDomFunctor(
  Γₛ⁰,
  Dict(
    1 => ACSetTransformation(Γₛ⁰[3], Γₛ⁰[1], V=[4, 2]),
    2 => ACSetTransformation(Γₛ⁰[3], Γₛ⁰[2], V=[1, 3]),
  ),
  ∫(Gₛ)
)

my_decomp2  = StrDecomp(Gₛ, Γₛ)

is_homomorphic(ob(colimit(my_decomp2)), K(2)) == false
@benchmark is_homomorphic(ob(colimit(my_decomp2)), K(2))

is_homomorphic(ob(colimit(my_decomp2)), K(3)) == true
@benchmark is_homomorphic(ob(colimit(my_decomp2)), K(3))

is_homomorphic(ob(colimit(my_decomp2)), K(4)) == true
@benchmark is_homomorphic(ob(colimit(my_decomp2)), K(4))

decide_sheaf_tree_shape(skeletalColoring(2), my_decomp2)[1] == false
@benchmark decide_sheaf_tree_shape(skeletalColoring(2), my_decomp2)[1]

decide_sheaf_tree_shape(skeletalColoring(3), my_decomp2)[1] == true
@benchmark decide_sheaf_tree_shape(skeletalColoring(3), my_decomp2)[1]

decide_sheaf_tree_shape(skeletalColoring(4), my_decomp2)[1] == true
@benchmark decide_sheaf_tree_shape(skeletalColoring(4), my_decomp2)[1]

# Benchmark 3 (large n per bag(20) and small bags k(2), 3 coloring)

#bag 1
H₁ = @acset Graph begin
  V = 20
  E = 34
  src = [1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13, 14, 15, 15, 16, 16, 17, 17, 18]
  tgt = [2, 3, 4, 5, 4, 6, 5, 7, 8, 9, 7, 10, 8, 11, 9, 12, 13, 14, 11, 15, 12, 15, 13, 16, 14, 17, 18, 16, 19, 17, 19, 18, 19, 20]
end
  
#adhesion 1,2
H₁₂ = @acset Graph begin
  V = 2
end
  
#bag 2
H₂ = @acset Graph begin
  V = 20
  E = 34
  src = [1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13, 14, 15, 15, 16, 16, 17, 17, 18]
  tgt = [2, 3, 4, 5, 4, 6, 5, 7, 8, 9, 7, 10, 8, 11, 9, 12, 13, 14, 11, 15, 12, 15, 13, 16, 14, 17, 18, 16, 19, 17, 19, 18, 19, 20]
end
  
Gₛ = @acset Graph begin
  V = 2
  E = 1
  src = [1]
  tgt = [2]
end
  
#transformations
Γₛ⁰ = Dict(1 => H₁, 2 => H₂, 3 => H₁₂)
Γₛ = FinDomFunctor(
  Γₛ⁰,
  Dict(
    1 => ACSetTransformation(Γₛ⁰[3], Γₛ⁰[1], V=[1, 2]),
    2 => ACSetTransformation(Γₛ⁰[3], Γₛ⁰[2], V=[1, 2]),
  ),
  ∫(Gₛ)
)

my_decomp3  = StrDecomp(Gₛ, Γₛ)

is_homomorphic(ob(colimit(my_decomp3)), K(2)) == false
@benchmark is_homomorphic(ob(colimit(my_decomp3)), K(2))

is_homomorphic(ob(colimit(my_decomp3)), K(3)) == true
@benchmark is_homomorphic(ob(colimit(my_decomp3)), K(3))

is_homomorphic(ob(colimit(my_decomp3)), K(4)) == true
@benchmark is_homomorphic(ob(colimit(my_decomp3)), K(4))

decide_sheaf_tree_shape(skeletalColoring(2), my_decomp3)[1] == false
@benchmark decide_sheaf_tree_shape(skeletalColoring(2), my_decomp3)[1]

decide_sheaf_tree_shape(skeletalColoring(3), my_decomp3)[1] == true
@benchmark decide_sheaf_tree_shape(skeletalColoring(3), my_decomp3)[1]

decide_sheaf_tree_shape(skeletalColoring(4), my_decomp3)[1] == true
@benchmark decide_sheaf_tree_shape(skeletalColoring(4), my_decomp3)[1]

# Benchmark 4 (small n per bag(5) with medium bags k(4), 3 coloring)

#bag 1
H₁ = @acset Graph begin
  V = 5
  E = 7
  src = [1, 2, 2, 2, 3, 4, 5]
  tgt = [2, 3, 4, 5, 4, 5, 1]
end

#adhesion 1,2
H₁₂ = @acset Graph begin
  V = 1
end

#bag 2
H₂ = @acset Graph begin
  V = 5
  E = 7
  src = [1, 2, 2, 2, 3, 4, 5]
  tgt = [2, 3, 4, 5, 4, 5, 1]
end

#adhesion 3,4
H₂₃ = @acset Graph begin
  V = 1
end

#bag 3
H₃ = @acset Graph begin
  V = 5
  E = 7
  src = [1, 2, 2, 2, 3, 4, 5]
  tgt = [2, 3, 4, 5, 4, 5, 1]
end

#adhesion 3,4
H₃₄ = @acset Graph begin
  V = 1
end

#bag 4
H₄ = @acset Graph begin
  V = 5
  E = 7
  src = [1, 2, 2, 2, 3, 4, 5]
  tgt = [2, 3, 4, 5, 4, 5, 1]
end

Gₛ = @acset Graph begin
  V = 4
  E = 3
  src = [1, 2, 3]
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

my_decomp4  = StrDecomp(Gₛ, Γₛ)

decide_sheaf_tree_shape(skeletalColoring(2), my_decomp4)[1] == false
@benchmark decide_sheaf_tree_shape(skeletalColoring(2), my_decomp4)[1]

decide_sheaf_tree_shape(skeletalColoring(3), my_decomp4)[1] == true
@benchmark decide_sheaf_tree_shape(skeletalColoring(3), my_decomp4)[1]

decide_sheaf_tree_shape(skeletalColoring(4), my_decomp4)[1] == true
@benchmark decide_sheaf_tree_shape(skeletalColoring(4), my_decomp4)[1]

# Benchmark 5 (medium n per bag(10) with medium bags k(4), 3 coloring)

#bag 1
H₁ = @acset Graph begin
  V = 10
  E = 16
  src = [1, 1, 1, 2, 2, 3, 4, 4, 4, 6, 6, 7, 7, 7, 8, 9]
  tgt = [2, 5, 6, 3, 8, 4, 5, 9, 10, 7, 10, 8, 9, 10, 9, 10]
end

#adhesion 1,2
H₁₂ = @acset Graph begin
  V = 1
end

#bag 2
H₂ = @acset Graph begin
  V = 10
  E = 18
  src = [1, 1, 2, 2, 2, 2, 3, 4, 4, 5, 5, 6, 6, 7, 7, 7, 8, 9]
  tgt = [2, 9, 3, 4, 9, 10, 4, 5, 10, 6, 10, 7, 10, 8, 9, 10, 9, 10]
end

#adhesion 3,4
H₂₃ = @acset Graph begin
  V = 1
end

#bag 3
H₃ = @acset Graph begin
  V = 10
  E = 16
  src = [1, 1, 1, 2, 2, 3, 4, 4, 4, 6, 6, 7, 7, 7, 8, 9]
  tgt = [2, 5, 6, 3, 8, 4, 5, 9, 10, 7, 10, 8, 9, 10, 9, 10]
end

#adhesion 3,4
H₃₄ = @acset Graph begin
  V = 1
end

#bag 4
H₄ = @acset Graph begin
  V = 10
  E = 18
  src = [1, 1, 2, 2, 2, 2, 3, 4, 4, 5, 5, 6, 6, 7, 7, 7, 8, 9]
  tgt = [2, 9, 3, 4, 9, 10, 4, 5, 10, 6, 10, 7, 10, 8, 9, 10, 9, 10]
end

Gₛ = @acset Graph begin
  V = 4
  E = 3
  src = [1, 3, 4]
  tgt = [2, 2, 3]
end

#transformations
Γₛ⁰ = Dict(1 => H₁, 2 => H₂, 3 => H₃, 4 => H₄, 5 => H₁₂, 6 => H₂₃, 7 => H₃₄)
Γₛ = FinDomFunctor(
  Γₛ⁰,
  Dict(
    1 => ACSetTransformation(Γₛ⁰[5], Γₛ⁰[1], V=[4]),
    2 => ACSetTransformation(Γₛ⁰[5], Γₛ⁰[2], V=[1]),
    3 => ACSetTransformation(Γₛ⁰[6], Γₛ⁰[2], V=[1]),
    4 => ACSetTransformation(Γₛ⁰[6], Γₛ⁰[3], V=[1]),
    5 => ACSetTransformation(Γₛ⁰[7], Γₛ⁰[3], V=[1]),
    6 => ACSetTransformation(Γₛ⁰[7], Γₛ⁰[4], V=[4])
  ),
  ∫(Gₛ)
)

my_decomp5  = StrDecomp(Gₛ, Γₛ)

decide_sheaf_tree_shape(skeletalColoring(2), my_decomp5)[1] == false
@benchmark decide_sheaf_tree_shape(skeletalColoring(2), my_decomp5)[1]

decide_sheaf_tree_shape(skeletalColoring(3), my_decomp5)[1] == true
@benchmark decide_sheaf_tree_shape(skeletalColoring(3), my_decomp5)[1]

decide_sheaf_tree_shape(skeletalColoring(4), my_decomp5)[1] == true
@benchmark decide_sheaf_tree_shape(skeletalColoring(4), my_decomp5)[1]

# Benchmark 6 (large n per bag(20) with medium bags k(4), 3 coloring)

#bag 1
H₁ = @acset Graph begin
  V = 20
  E = 34
  src = [1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13, 14, 15, 15, 16, 16, 17, 17, 18]
  tgt = [2, 3, 4, 5, 4, 6, 5, 7, 8, 9, 7, 10, 8, 11, 9, 12, 13, 14, 11, 15, 12, 15, 13, 16, 14, 17, 18, 16, 19, 17, 19, 18, 19, 20]
end

#adhesion 1,2
H₁₂ = @acset Graph begin
  V = 1
end

#bag 2
H₂ = @acset Graph begin
  V = 20
  E = 34
  src = [1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13, 14, 15, 15, 16, 16, 17, 17, 18]
  tgt = [2, 3, 4, 5, 4, 6, 5, 7, 8, 9, 7, 10, 8, 11, 9, 12, 13, 14, 11, 15, 12, 15, 13, 16, 14, 17, 18, 16, 19, 17, 19, 18, 19, 20]
end

#adhesion 3,4
H₂₃ = @acset Graph begin
  V = 1
end

#bag 3
H₃ = @acset Graph begin
  V = 20
  E = 34
  src = [1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13, 14, 15, 15, 16, 16, 17, 17, 18]
  tgt = [2, 3, 4, 5, 4, 6, 5, 7, 8, 9, 7, 10, 8, 11, 9, 12, 13, 14, 11, 15, 12, 15, 13, 16, 14, 17, 18, 16, 19, 17, 19, 18, 19, 20]
end

#adhesion 3,4
H₃₄ = @acset Graph begin
  V = 1
end

#bag 4
H₄ = @acset Graph begin
  V = 20
  E = 34
  src = [1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13, 14, 15, 15, 16, 16, 17, 17, 18]
  tgt = [2, 3, 4, 5, 4, 6, 5, 7, 8, 9, 7, 10, 8, 11, 9, 12, 13, 14, 11, 15, 12, 15, 13, 16, 14, 17, 18, 16, 19, 17, 19, 18, 19, 20]
end

Gₛ = @acset Graph begin
  V = 4
  E = 3
  src = [1, 2, 3]
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

my_decomp6  = StrDecomp(Gₛ, Γₛ)

decide_sheaf_tree_shape(skeletalColoring(2), my_decomp6)[1] == false
@benchmark decide_sheaf_tree_shape(skeletalColoring(2), my_decomp6)[1]

decide_sheaf_tree_shape(skeletalColoring(3), my_decomp6)[1] == true
@benchmark decide_sheaf_tree_shape(skeletalColoring(3), my_decomp6)[1]

decide_sheaf_tree_shape(skeletalColoring(4), my_decomp6)[1] == true
@benchmark decide_sheaf_tree_shape(skeletalColoring(4), my_decomp6)[1]

# Benchmark 7 (small n per bag(5) with large bags k(10), 3 coloring)

#bag 1
H₁ = @acset Graph begin
  V = 5
  E = 7
  src = [1, 2, 2, 2, 3, 4, 5]
  tgt = [2, 3, 4, 5, 4, 5, 1]
end

#adhesion 1,2
H₁₂ = @acset Graph begin
  V = 1
end

#bag 2
H₂ = @acset Graph begin
  V = 5
  E = 7
  src = [1, 2, 2, 2, 3, 4, 5]
  tgt = [2, 3, 4, 5, 4, 5, 1]
end

#adhesion 3,4
H₂₃ = @acset Graph begin
  V = 1
end

#bag 3
H₃ = @acset Graph begin
  V = 5
  E = 7
  src = [1, 2, 2, 2, 3, 4, 5]
  tgt = [2, 3, 4, 5, 4, 5, 1]
end

#adhesion 3,4
H₃₄ = @acset Graph begin
  V = 1
end

#bag 4
H₄ = @acset Graph begin
  V = 5
  E = 7
  src = [1, 2, 2, 2, 3, 4, 5]
  tgt = [2, 3, 4, 5, 4, 5, 1]
end

#adhesion 4,5
H₄₅ = @acset Graph begin
  V = 1
end

#bag 5
H₅ = @acset Graph begin
  V = 5
  E = 7
  src = [1, 2, 2, 2, 3, 4, 5]
  tgt = [2, 3, 4, 5, 4, 5, 1]
end

#adhesion 5,6
H₅₆ = @acset Graph begin
  V = 1
end

#bag 6
H₆ = @acset Graph begin
  V = 5
  E = 7
  src = [1, 2, 2, 2, 3, 4, 5]
  tgt = [2, 3, 4, 5, 4, 5, 1]
end

#adhesion 6,7
H₆₇ = @acset Graph begin
  V = 1
end

#bag 7
H₇ = @acset Graph begin
  V = 5
  E = 7
  src = [1, 2, 2, 2, 3, 4, 5]
  tgt = [2, 3, 4, 5, 4, 5, 1]
end

#adhesion 7,8
H₇₈ = @acset Graph begin
  V = 1
end

#bag 8
H₈ = @acset Graph begin
  V = 5
  E = 7
  src = [1, 2, 2, 2, 3, 4, 5]
  tgt = [2, 3, 4, 5, 4, 5, 1]
end

#adhesion 8,9
H₈₉ = @acset Graph begin
  V = 1
end

#bag 9
H₉ = @acset Graph begin
  V = 5
  E = 7
  src = [1, 2, 2, 2, 3, 4, 5]
  tgt = [2, 3, 4, 5, 4, 5, 1]
end

#adhesion 9,10
H₉₁₀   = @acset Graph begin
  V = 1
end

#bag 10
H₁₀ = @acset Graph begin
  V = 5
  E = 7
  src = [1, 2, 2, 2, 3, 4, 5]
  tgt = [2, 3, 4, 5, 4, 5, 1]
end

Gₛ = @acset Graph begin
  V = 10
  E = 9
  src = [1, 2, 3, 4, 5, 6, 7, 8, 9]
  tgt = [2, 3, 4, 5, 6, 7, 8, 9, 10]
end

#transformations
Γₛ⁰ = Dict(1 => H₁, 2 => H₂, 3 => H₃, 4 => H₄, 5 => H₅, 6 => H₆, 7 => H₇, 8 => H₈, 9 => H₉, 10 => H₁₀, 11 => H₁₂, 12 => H₂₃, 13 => H₃₄, 14 => H₄₅,
          15 => H₅₆, 16 => H₆₇, 17 => H₇₈, 18 => H₈₉, 19 => H₉₁₀)
Γₛ = FinDomFunctor(
  Γₛ⁰,
  Dict(
    1 => ACSetTransformation(Γₛ⁰[11], Γₛ⁰[1], V=[1]),
    2 => ACSetTransformation(Γₛ⁰[11], Γₛ⁰[2], V=[1]),
    3 => ACSetTransformation(Γₛ⁰[12], Γₛ⁰[2], V=[1]),
    4 => ACSetTransformation(Γₛ⁰[12], Γₛ⁰[3], V=[1]),
    5 => ACSetTransformation(Γₛ⁰[13], Γₛ⁰[3], V=[1]),
    6 => ACSetTransformation(Γₛ⁰[13], Γₛ⁰[4], V=[1]),
    7 => ACSetTransformation(Γₛ⁰[14], Γₛ⁰[4], V=[1]),
    8 => ACSetTransformation(Γₛ⁰[14], Γₛ⁰[5], V=[1]),
    9 => ACSetTransformation(Γₛ⁰[15], Γₛ⁰[5], V=[1]),
    10 => ACSetTransformation(Γₛ⁰[15], Γₛ⁰[6], V=[1]),
    11 => ACSetTransformation(Γₛ⁰[16], Γₛ⁰[6], V=[1]),
    12 => ACSetTransformation(Γₛ⁰[16], Γₛ⁰[7], V=[1]),
    13 => ACSetTransformation(Γₛ⁰[17], Γₛ⁰[7], V=[1]),
    14 => ACSetTransformation(Γₛ⁰[17], Γₛ⁰[8], V=[1]),
    15 => ACSetTransformation(Γₛ⁰[18], Γₛ⁰[8], V=[1]),
    16 => ACSetTransformation(Γₛ⁰[18], Γₛ⁰[9], V=[1]),
    17 => ACSetTransformation(Γₛ⁰[19], Γₛ⁰[9], V=[1]),
    18 => ACSetTransformation(Γₛ⁰[19], Γₛ⁰[10], V=[1])
  ),
  ∫(Gₛ)
)

my_decomp7  = StrDecomp(Gₛ, Γₛ)

decide_sheaf_tree_shape(skeletalColoring(2), my_decomp7)[1] == false
@benchmark decide_sheaf_tree_shape(skeletalColoring(2), my_decomp7)[1]

decide_sheaf_tree_shape(skeletalColoring(3), my_decomp7)[1] == true
@benchmark decide_sheaf_tree_shape(skeletalColoring(3), my_decomp7)[1]

decide_sheaf_tree_shape(skeletalColoring(4), my_decomp7)[1] == true
@benchmark decide_sheaf_tree_shape(skeletalColoring(4), my_decomp7)[1]

# Benchmark 8 (medium n per bag(10) with large bags k(10), 3 coloring)

# Benchmark 9 (large n per bag(20) with large bags k(10), 3 coloring)