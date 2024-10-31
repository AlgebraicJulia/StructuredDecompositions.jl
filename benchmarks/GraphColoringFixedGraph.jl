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

el = GraphsPkg.Edge.([ (1,2), (1,3), (2,3), (2,4), (3,4), (4,5), (4,6), (5,6), (5,7), (6,7), (7,8), (7,9), (8,9), (8,10), (9,10), (10,11), 
(10,12), (10,21), (10,22), (10,31), (10,32), (11,12), (11,13), (12,13), (13,14), (13,15), (14,15), (14,16), (15,16), (16,17), (16,18),
(17,18), (17,19), (18,19), (19,20), (21,22), (21,23), (22,23), (23,24), (23,25), (24,25), (24,26), (25,26), (26,27), (26,28), (27,28), 
(27,29), (28,29), (29,20), (31,32), (31,33), (32,33), (33,34), (33,35), (34,35), (34,36), (35,36), (36,37), (36,38), (37,38), (37,39),
(38,39), (39,40)])

graph = GraphsPkg.SimpleGraph(el)
@benchmark GraphsPkg.degree_greedy_color(graph)