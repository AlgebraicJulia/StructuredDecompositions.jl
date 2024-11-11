# This file examines the relationship between run time, bag size, and number of bags for a fixed 40 node graph.
# Note since bags must overlap, the colorings will need to deal with more than 40 total nodes due to overlap.
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

# We will consider cases with 2(even), 4(even), 8(uneven), 16(uneven) bags
# We will first benchmark the coloring algorithm in graphs.jl

# el = GraphsPkg.Edge.([ (1,2), (1,3), (2,3), (2,4), (3,4), (4,5), (4,6), (5,6), (5,7), (6,7), (7,8), (7,9), (8,9), (8,10), (9,10), (10,11), 
# (10,12), (10,21), (10,22), (10,31), (10,32), (11,12), (11,13), (12,13), (13,14), (13,15), (14,15), (14,16), (15,16), (16,17), (16,18),
# (17,18), (17,19), (18,19), (19,20), (21,22), (21,23), (22,23), (23,24), (23,25), (24,25), (24,26), (25,26), (26,27), (26,28), (27,28), 
# (27,29), (28,29), (29,20), (31,32), (31,33), (32,33), (33,34), (33,35), (34,35), (34,36), (35,36), (36,37), (36,38), (37,38), (37,39),
# (38,39), (39,40)])

# graph = GraphsPkg.SimpleGraph(el)
# @benchmark GraphsPkg.degree_greedy_color(graph)

# 1 bag case (40 nodes, no adhesions)

#bag 1
bags1 = @acset Graph begin
  V = 40
  E = 63
  src = [1, 1, 2, 2, 3, 4, 4, 5, 5, 6, 7, 7, 8, 8, 9, 10, 10, 10, 10, 10, 10, 11, 11, 12, 13, 13, 14, 14, 15, 16, 16, 17, 17, 18, 19,
        21, 21, 22, 23, 23, 24, 24, 25, 26, 26, 27, 27, 28, 29, 31, 31, 32, 33, 33, 34, 34, 35, 36, 36, 37, 37, 38, 39]
  tgt = [2, 3, 3, 4, 4, 5, 6, 6, 7, 7, 8, 9, 9, 10, 10, 11, 12, 21, 22, 31, 32, 12, 13, 13, 14, 15, 15, 16, 16, 17, 18, 18, 19, 19, 20,
        22, 23, 23, 24, 25, 25, 26, 26, 27, 28, 28, 29, 29, 30, 32, 33, 33, 34, 35, 35, 36, 36, 37, 38, 38, 39, 39, 40]
end

bags2 = @acset Graph begin
  V = 1
  E = 0
  src = []
  tgt = []
end

bagsad = @acset Graph begin
  V = 1
  E = 0
  src = []
  tgt = []
end

Gₛ = @acset Graph begin
  V = 2
  E = 1
  src = [1]
  tgt = [1]
end

#transformations
Γₛ⁰ = Dict(1 => bags1, 2 => bags2, 3 => bagsad)
Γₛ = FinDomFunctor(
  Γₛ⁰,
  Dict(
   1 => ACSetTransformation(Γₛ⁰[3], Γₛ⁰[1], V=[1]),
   2 => ACSetTransformation(Γₛ⁰[3], Γₛ⁰[2], V=[1]),
 ),
 ∫(Gₛ)
)

bags1decomp  = StrDecomp(Gₛ, Γₛ)

@benchmark decide_sheaf_tree_shape(skeletalColoring(2), bags1decomp)[1]
@benchmark decide_sheaf_tree_shape(skeletalColoring(3), bags1decomp)[1]
@benchmark is_homomorphic(K(2), bag1)
@benchmark is_homomorphic(K(3), bag1)
@benchmark is_homomorphic(K(4), bag1)

# 2 bag case (21 nodes each, 1 node adhesion)

#bag 1
H₁ = @acset Graph begin
  V = 21
  E = 31
  src = [1, 1, 2, 2, 3, 4, 4, 5, 5, 6, 7, 7, 8, 8, 9, 10, 10, 11, 11, 12, 13, 13, 14, 14, 15, 16, 16, 17, 17, 18, 19]
  tgt = [2, 3, 3, 4, 4, 5, 6, 6, 7, 7, 8, 9, 9, 10, 10, 11, 12, 12, 13, 13, 14, 15, 15, 16, 16, 17, 18, 18, 19, 19, 20]
end

#adhesion 1,2
H₁₂ = @acset Graph begin
  V = 1
end

#bag 2
H₂ = @acset Graph begin
  V = 21
  E = 31
  src = [1, 1, 2, 2, 3, 4, 4, 5, 5, 6, 7, 7, 8, 8, 9, 10, 10, 11, 11, 12, 13, 13, 14, 14, 15, 16, 16, 17, 17, 18, 19]
  tgt = [2, 3, 3, 4, 4, 5, 6, 6, 7, 7, 8, 9, 9, 10, 10, 11, 12, 12, 13, 13, 14, 15, 15, 16, 16, 17, 18, 18, 19, 19, 20]
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

bags2decomp  = StrDecomp(Gₛ, Γₛ)

# 4 bag case (1 bag w/ 10 nodes, 3 bags w/ 11 nodes, 1 node adhesion)

#bag 1
H₁ = @acset Graph begin
  V = 10
  E = 15
  src = [1, 1, 2, 2, 3, 4, 4, 5, 5, 6, 7, 7, 8, 8, 9]
  tgt = [2, 3, 3, 4, 4, 5, 6, 6, 7, 7, 8, 9, 9, 10, 10]
end

#adhesion 1,2
H₁₂ = @acset Graph begin
  V = 1
end

#bag 2
H₂ = @acset Graph begin
  V = 11
  E = 16
  src = [1, 1, 2, 2, 3, 4, 4, 5, 5, 6, 7, 7, 8, 8, 9, 10]
  tgt = [2, 3, 3, 4, 4, 5, 6, 6, 7, 7, 8, 9, 9, 10, 10, 11]
end

#adhesion 3,4
H₂₃ = @acset Graph begin
  V = 1
end

#bag 3
H₃ = @acset Graph begin
  V = 11
  E = 16
  src = [1, 1, 2, 2, 3, 4, 4, 5, 5, 6, 7, 7, 8, 8, 9, 10]
  tgt = [2, 3, 3, 4, 4, 5, 6, 6, 7, 7, 8, 9, 9, 10, 10, 11]
end

#adhesion 3,4
H₃₄ = @acset Graph begin
  V = 1
end

#bag 4
H₄ = @acset Graph begin
  V = 11
  E = 16
  src = [1, 1, 2, 2, 3, 4, 4, 5, 5, 6, 7, 7, 8, 8, 9, 10]
  tgt = [2, 3, 3, 4, 4, 5, 6, 6, 7, 7, 8, 9, 9, 10, 10, 11]
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

bags4decomp  = StrDecomp(Gₛ, Γₛ)

# 8 bag case (bags between 4-7 nodes, 47 total nodes, 7 total adhesion nodes)

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
  V = 7
  E = 10
  src = [1, 1, 2, 2, 3, 4, 4, 5, 5, 6]
  tgt = [2, 3, 3, 4, 4, 5, 6, 6, 7, 7]
end

#adhesion 3,4
H₂₃ = @acset Graph begin
  V = 1
end

#bag 3
H₃ = @acset Graph begin
  V = 7
  E = 10
  src = [1, 1, 2, 2, 3, 4, 4, 5, 5, 6]
  tgt = [2, 3, 3, 4, 4, 5, 6, 6, 7, 7]
end

#adhesion 3,4
H₃₄ = @acset Graph begin
  V = 1
end

#bag 4
H₄ = @acset Graph begin
  V = 5
  E = 6
  src = [1, 1, 2, 2, 3, 4]
  tgt = [2, 3, 3, 4, 4, 5]
end

#adhesion 2,5
H₂₅ = @acset Graph begin
  V = 1
end

#bag 5
H₅ = @acset Graph begin
  V = 7
  E = 10
  src = [1, 1, 2, 2, 3, 4, 4, 5, 5, 6]
  tgt = [2, 3, 3, 4, 4, 5, 6, 6, 7, 7]
end

#adhesion 5,6
H₅₆ = @acset Graph begin
  V = 1
end

#bag 6
H₆ = @acset Graph begin
  V = 5
  E = 6
  src = [1, 1, 2, 2, 3, 4]
  tgt = [2, 3, 3, 4, 4, 5]
end

#adhesion 2,7
H₂₇ = @acset Graph begin
  V = 1
end

#bag 7
H₇ = @acset Graph begin
  V = 7
  E = 10
  src = [1, 1, 2, 2, 3, 4, 4, 5, 5, 6]
  tgt = [2, 3, 3, 4, 4, 5, 6, 6, 7, 7]
end

#adhesion 7,8
H₇₈ = @acset Graph begin
  V = 1
end

#bag 8
H₈ = @acset Graph begin
  V = 5
  E = 6
  src = [1, 1, 2, 2, 3, 4]
  tgt = [2, 3, 3, 4, 4, 5]
end

Gₛ = @acset Graph begin
  V = 8
  E = 7
  src = [1, 2, 2, 2, 3, 5, 7]
  tgt = [2, 3, 5, 7, 4, 6, 8]
end

#transformations
Γₛ⁰ = Dict(1 => H₁, 2 => H₂, 3 => H₃, 4 => H₄, 5 => H₅, 6 => H₆, 7 => H₇, 8 => H₈, 9 => H₁₂, 10 => H₂₃, 11 => H₃₄, 12 => H₂₅,
          13 => H₅₆, 14 => H₂₇, 15 => H₇₈)
Γₛ = FinDomFunctor(
  Γₛ⁰,
  Dict(
    1 => ACSetTransformation(Γₛ⁰[9], Γₛ⁰[1], V=[4]),
    2 => ACSetTransformation(Γₛ⁰[9], Γₛ⁰[2], V=[1]),
    3 => ACSetTransformation(Γₛ⁰[10], Γₛ⁰[2], V=[7]),
    4 => ACSetTransformation(Γₛ⁰[10], Γₛ⁰[3], V=[1]),
    5 => ACSetTransformation(Γₛ⁰[11], Γₛ⁰[3], V=[7]),
    6 => ACSetTransformation(Γₛ⁰[11], Γₛ⁰[4], V=[1]),
    7 => ACSetTransformation(Γₛ⁰[12], Γₛ⁰[2], V=[7]),
    8 => ACSetTransformation(Γₛ⁰[12], Γₛ⁰[5], V=[1]),
    9 => ACSetTransformation(Γₛ⁰[13], Γₛ⁰[5], V=[7]),
    10 => ACSetTransformation(Γₛ⁰[13], Γₛ⁰[6], V=[1]),
    11 => ACSetTransformation(Γₛ⁰[14], Γₛ⁰[2], V=[7]),
    12 => ACSetTransformation(Γₛ⁰[14], Γₛ⁰[7], V=[1]),
    13 => ACSetTransformation(Γₛ⁰[15], Γₛ⁰[7], V=[7]),
    14 => ACSetTransformation(Γₛ⁰[15], Γₛ⁰[8], V=[1]),
  ),
  ∫(Gₛ)
)

bags8decomp  = StrDecomp(Gₛ, Γₛ)

# 12 bag case (4-5 nodes per bag, 51 total nodes, 11 total adhesion nodes)

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
  V = 5
  E = 6
  src = [1, 1, 2, 2, 3, 4]
  tgt = [2, 3, 3, 4, 4, 5]
end

#adhesion 3,7
H₃₇ = @acset Graph begin
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

#adhesion 8,9
H₈₉ = @acset Graph begin
  V = 1
end

#bag 9
H₉ = @acset Graph begin
  V = 5
  E = 6
  src = [1, 1, 2, 2, 3, 4]
  tgt = [2, 3, 3, 4, 4, 5]
end

#adhesion 3,10
H₃₁₀   = @acset Graph begin
  V = 1
end

#bag 10
H₁₀ = @acset Graph begin
  V = 4
  E = 5
  src = [1, 1, 2, 2, 3]
  tgt = [2, 3, 3, 4, 4]
end

#adhesion 10,11
H₁₀₁₁   = @acset Graph begin
  V = 1
end

#bag 11
H₁₁ = @acset Graph begin
  V = 4
  E = 5
  src = [1, 1, 2, 2, 3]
  tgt = [2, 3, 3, 4, 4]
end

#adhesion 11,12
H₁₁₁₂   = @acset Graph begin
  V = 1
end

#bag 12
H₁₂₂ = @acset Graph begin
  V = 5
  E = 6
  src = [1, 1, 2, 2, 3, 4]
  tgt = [2, 3, 3, 4, 4, 5]
end

Gₛ = @acset Graph begin
  V = 12
  E = 11
  src = [1, 2, 3, 3, 3, 4, 5, 7, 8, 10, 11]
  tgt = [2, 3, 4, 7, 10, 5, 6, 8, 9, 11, 12]
end

#transformations
Γₛ⁰ = Dict(1 => H₁, 2 => H₂, 3 => H₃, 4 => H₄, 5 => H₅, 6 => H₆, 7 => H₇, 8 => H₈, 9 => H₉, 10 => H₁₀, 11 => H₁₁, 12 => H₁₂₂, 
          13 => H₁₂, 14 => H₂₃, 15 => H₃₄, 16 => H₄₅, 17 => H₅₆, 18 => H₃₇, 19 => H₇₈, 20 => H₈₉, 21 => H₃₁₀, 22 => H₁₀₁₁, 23 => H₁₁₁₂)
Γₛ = FinDomFunctor(
  Γₛ⁰,
  Dict(
    1 => ACSetTransformation(Γₛ⁰[13], Γₛ⁰[1], V=[4]),
    2 => ACSetTransformation(Γₛ⁰[13], Γₛ⁰[2], V=[1]),
    3 => ACSetTransformation(Γₛ⁰[14], Γₛ⁰[2], V=[4]),
    4 => ACSetTransformation(Γₛ⁰[14], Γₛ⁰[3], V=[1]),
    5 => ACSetTransformation(Γₛ⁰[15], Γₛ⁰[3], V=[4]),
    6 => ACSetTransformation(Γₛ⁰[15], Γₛ⁰[4], V=[1]),
    7 => ACSetTransformation(Γₛ⁰[16], Γₛ⁰[4], V=[4]),
    8 => ACSetTransformation(Γₛ⁰[16], Γₛ⁰[5], V=[1]),
    9 => ACSetTransformation(Γₛ⁰[17], Γₛ⁰[5], V=[4]),
    10 => ACSetTransformation(Γₛ⁰[17], Γₛ⁰[6], V=[1]),
    11 => ACSetTransformation(Γₛ⁰[18], Γₛ⁰[3], V=[4]),
    12 => ACSetTransformation(Γₛ⁰[18], Γₛ⁰[7], V=[1]),
    13 => ACSetTransformation(Γₛ⁰[19], Γₛ⁰[7], V=[4]),
    14 => ACSetTransformation(Γₛ⁰[19], Γₛ⁰[8], V=[1]),
    15 => ACSetTransformation(Γₛ⁰[20], Γₛ⁰[8], V=[4]),
    16 => ACSetTransformation(Γₛ⁰[20], Γₛ⁰[9], V=[1]),
    17 => ACSetTransformation(Γₛ⁰[21], Γₛ⁰[3], V=[4]),
    18 => ACSetTransformation(Γₛ⁰[21], Γₛ⁰[10], V=[1]),
    19 => ACSetTransformation(Γₛ⁰[22], Γₛ⁰[10], V=[4]),
    20 => ACSetTransformation(Γₛ⁰[22], Γₛ⁰[11], V=[1]),
    21 => ACSetTransformation(Γₛ⁰[23], Γₛ⁰[11], V=[4]),
    22 => ACSetTransformation(Γₛ⁰[23], Γₛ⁰[12], V=[1])
  ),
  ∫(Gₛ)
)

bags12decomp  = StrDecomp(Gₛ, Γₛ)

decide_sheaf_tree_shape(skeletalColoring(2), bags2decomp)[1] == false
decide_sheaf_tree_shape(skeletalColoring(2), bags4decomp)[1] == false
decide_sheaf_tree_shape(skeletalColoring(2), bags8decomp)[1] == false
decide_sheaf_tree_shape(skeletalColoring(2), bags12decomp)[1] == false
@benchmark decide_sheaf_tree_shape(skeletalColoring(2), bags2decomp)[1]
@benchmark decide_sheaf_tree_shape(skeletalColoring(2), bags4decomp)[1]
@benchmark decide_sheaf_tree_shape(skeletalColoring(2), bags8decomp)[1]
@benchmark decide_sheaf_tree_shape(skeletalColoring(2), bags12decomp)[1]

decide_sheaf_tree_shape(skeletalColoring(3), bags2decomp)[1] == true
decide_sheaf_tree_shape(skeletalColoring(3), bags4decomp)[1] == true
decide_sheaf_tree_shape(skeletalColoring(3), bags8decomp)[1] == true
decide_sheaf_tree_shape(skeletalColoring(3), bags12decomp)[1] == true
@benchmark decide_sheaf_tree_shape(skeletalColoring(3), bags2decomp)[1]
@benchmark decide_sheaf_tree_shape(skeletalColoring(3), bags4decomp)[1]
@benchmark decide_sheaf_tree_shape(skeletalColoring(3), bags8decomp)[1]
@benchmark decide_sheaf_tree_shape(skeletalColoring(3), bags12decomp)[1]

decide_sheaf_tree_shape(skeletalColoring(4), bags2decomp)[1] == true
decide_sheaf_tree_shape(skeletalColoring(4), bags4decomp)[1] == true
decide_sheaf_tree_shape(skeletalColoring(4), bags8decomp)[1] == true
decide_sheaf_tree_shape(skeletalColoring(4), bags12decomp)[1] == true
@benchmark decide_sheaf_tree_shape(skeletalColoring(4), bags2decomp)[1]
@benchmark decide_sheaf_tree_shape(skeletalColoring(4), bags4decomp)[1]
@benchmark decide_sheaf_tree_shape(skeletalColoring(4), bags8decomp)[1]
@benchmark decide_sheaf_tree_shape(skeletalColoring(4), bags12decomp)[1]