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

# colorability_test(n, the_test_case) = is_homomorphic(ob(colimit(the_test_case)), K(n)) == decide_sheaf_tree_shape(skeletalColoring(n), the_test_case)[1]

function adhesion_filter(tup::Tuple, d::StructuredDecomposition)
  # if d.decomp_type == Decomposition
  #   error("expecting ", CoDecomposition, " given ", Decomposition)
  # end
  # # d_csp is the cospan dx₁ -> de <- dx₂ corresp to some edge e = x₁x₂ in shape(d)
  # (csp, d_csp)      = tup  #unpack the tuple
  # # the pullback cone dx₁ <-l₁-- p --l₂ --> dx₂ with legs l₁ and l₂
  # p_cone            = pullback(d_csp)
  # p_legs            = legs(p_cone)
  # #for each leg lᵢ : p → xᵢ of the pullback cone, 
  # #compute its image ιᵢ : im lᵢ → dxᵢ
  # imgs              = map( f -> legs(image(f))[1], p_legs)
  # #now get the new desired cospan; 
  # #i.e.  im l₁ --ι₁--> dx₁ --l₁--> de <--l₂--dx₂ <--ι₂-- im l₂
  # new_d_csp         = map(t -> compose(t...), zip(imgs, d_csp))  
  # #get the domain of d 
  # d_dom             = dom(d.diagram)
  # #now make the new decomposition, call it δ
  # #start with the object map δ₀
  # function ob_replace(x)
  #   if x == dom(d_dom, csp[1])
  #     dom(new_d_csp[1])
  #   elseif x == dom(d_dom, csp[2])
  #     dom(new_d_csp[2])
  #   else 
  #     ob_map(d,x) 
  #   end
  # end
  # δ₀ = Dict( x => ob_replace(x) for x ∈ ob_generators(d_dom) )
  # #now do the same thing with the morphism map
  # function mor_replace(f) 
  #   if f == csp[1]
  #     return new_d_csp[1]
  #   elseif f == csp[2]
  #     return new_d_csp[2]
  #   else
  #     return hom_map(d,f)
  #   end 
  # end
  # δ₁ = Dict( f => mor_replace(f) for f ∈ hom_generators(d_dom) )
  # StrDecomp(d.decomp_shape, FinDomFunctor(δ₀, δ₁, d.domain), d.decomp_type)

  #attempt 3

  if d.decomp_type == Decomposition
    error("expecting ", CoDecomposition, " given ", Decomposition)
  end

  #unpack
  (csp, d_csp) = tup
  
  p_cone = pullback(d_csp)

  if isempty(p_cone)
    d_dom = FinSet(0)
  else
    p_legs = legs(p_cone)
    imgs = map(f->legs(image(f))[1], p_legs)
    new_d_csp = map(t->compose(t...), zip(imgs, d_csp))
  end

  d_dom = dom(d.diagram)

  function ob_replace(x)
    if x == dom(d_dom, csp[1])
      dom(new_d_csp[1])
    elseif x == dom(d_dom, csp[2])
      dom(new_d_csp[2])
    else
      ob_map(d,x)
    end
  end

  function mor_replace(f)
    if f == csp[1]
      return new_d_csp[1]
    elseif f == csp[2]
      return new_d_csp[2]
    else
      return hom_map(d,f)
    end
  end
  
  δ₀ = Dict( x => ob_replace(x) for x ∈ ob_generators(d_dom))
  δ₁ = Dict( f => mor_replace(f) for f ∈ hom_generators(d_dom))

  StrDecomp(d.decomp_shape, FinDomFunctor(δ₀, δ₁, d.domain), d.decomp_type)
end

#for some reason PartialFunctions is giving me an error here 
#and we have to explicitly Curry adhesion_filter.. 
# adhesion_filter(tup::Tuple) = d -> adhesion_filter(tup, d) 
#fixed adhesion_filter so no longer needed

"""Solve the decision problem encoded by a sheaf. 
The algorithm is as follows: 
  compute on each bag (optionally, if the decomposition of the solution space
                        is already known, then it can be passed as an argument),
  compute composites on edges, 
  project back down to bags
  answer (providing a witness)
    "no" if there is an empty bag; 
    "yes" otherwise.
"""
function decide_sheaf_tree_shape(f, d::StructuredDecomposition, solution_space_decomp::StructuredDecomposition = 𝐃(f, d, CoDecomposition))
  # witness = foldl(∘, map(adhesion_filter, adhesionSpans(solution_space_decomp, true)))(solution_space_decomp)

  # # witness = 
  # # ∘((adhesion_filter.(adhesionSpans(solution_space_decomp, true)))...)(solution_space_decomp)

  # # (foldr(&, map( !isempty, bags(witness))), witness)

  # (all(!isempty, bags(witness)), witness)

  witness = solution_space_decomp
  adhesion_spans = adhesionSpans(solution_space_decomp, true)
  # println(adhesion_spans)
  # count = 0
  for adhesion in adhesion_spans
    # count = count + 1
    witness = adhesion_filter(adhesion, witness)
    if any(isempty, bags(witness))
      return (false, witness)
    end
  end
  # println(count)
  return (true, witness)
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


# Benchmark 1-3 minimum and median homomorphic comparisons

# all of these are returning false
min_2_test_1_h = minimum(@benchmark is_homomorphic(ob(colimit(my_decomp1)), K(2)))
min_2_test_2_h = minimum(@benchmark is_homomorphic(ob(colimit(my_decomp2)), K(2)))
min_2_test_3_h = minimum(@benchmark is_homomorphic(ob(colimit(my_decomp3)), K(2)))
ratio(min_2_test_2_h, min_2_test_1_h)
ratio(min_2_test_3_h, min_2_test_2_h)
ratio(min_2_test_3_h, min_2_test_1_h)

med_2_test_1_h = median(@benchmark is_homomorphic(ob(colimit(my_decomp1)), K(2)))
med_2_test_2_h = median(@benchmark is_homomorphic(ob(colimit(my_decomp2)), K(2)))
med_2_test_3_h = median(@benchmark is_homomorphic(ob(colimit(my_decomp3)), K(2)))
ratio(med_2_test_2_h, med_2_test_1_h)
ratio(med_2_test_3_h, med_2_test_2_h)
ratio(med_2_test_3_h, med_2_test_1_h)

# all of these are returning true
min_3_test_1_h = minimum(@benchmark is_homomorphic(ob(colimit(my_decomp1)), K(3)))
min_3_test_2_h = minimum(@benchmark is_homomorphic(ob(colimit(my_decomp2)), K(3)))
min_3_test_3_h = minimum(@benchmark is_homomorphic(ob(colimit(my_decomp3)), K(3)))
ratio(min_3_test_2_h, min_3_test_1_h)
ratio(min_3_test_3_h, min_3_test_2_h)
ratio(min_3_test_3_h, min_3_test_1_h)

med_3_test_1_h = median(@benchmark is_homomorphic(ob(colimit(my_decomp1)), K(3)))
med_3_test_2_h = median(@benchmark is_homomorphic(ob(colimit(my_decomp2)), K(3)))
med_3_test_3_h = median(@benchmark is_homomorphic(ob(colimit(my_decomp3)), K(3)))
ratio(med_3_test_2_h, med_3_test_1_h)
ratio(med_3_test_3_h, med_3_test_2_h)
ratio(med_3_test_3_h, med_3_test_1_h)

# all of these are returning true
min_4_test_1_h = minimum(@benchmark is_homomorphic(ob(colimit(my_decomp1)), K(4)))
min_4_test_2_h = minimum(@benchmark is_homomorphic(ob(colimit(my_decomp2)), K(4)))
min_4_test_3_h = minimum(@benchmark is_homomorphic(ob(colimit(my_decomp3)), K(4)))
ratio(min_4_test_2_h, min_4_test_1_h)
ratio(min_4_test_3_h, min_4_test_2_h)
ratio(min_4_test_3_h, min_4_test_1_h)

med_4_test_1_h = median(@benchmark is_homomorphic(ob(colimit(my_decomp1)), K(4)))
med_4_test_2_h = median(@benchmark is_homomorphic(ob(colimit(my_decomp2)), K(4)))
med_4_test_3_h = median(@benchmark is_homomorphic(ob(colimit(my_decomp3)), K(4)))
ratio(med_4_test_2_h, med_4_test_1_h)
ratio(med_4_test_3_h, med_4_test_2_h)
ratio(med_4_test_3_h, med_4_test_1_h)


# Benchmark 1-3 minimum and median sheaf tree comparisons

# all of these are returning false
min_2_test_1_s = minimum(@benchmark decide_sheaf_tree_shape(skeletalColoring(2), my_decomp1)[1])
min_2_test_2_s = minimum(@benchmark decide_sheaf_tree_shape(skeletalColoring(2), my_decomp2)[1])
min_2_test_3_s = minimum(@benchmark decide_sheaf_tree_shape(skeletalColoring(2), my_decomp3)[1])
ratio(min_2_test_2_s, min_2_test_1_s)
ratio(min_2_test_3_s, min_2_test_2_s)
ratio(min_2_test_3_s, min_2_test_1_s)

med_2_test_1_s = median(@benchmark decide_sheaf_tree_shape(skeletalColoring(2), my_decomp1)[1])
med_2_test_2_s = median(@benchmark decide_sheaf_tree_shape(skeletalColoring(2), my_decomp2)[1])
med_2_test_3_s = median(@benchmark decide_sheaf_tree_shape(skeletalColoring(2), my_decomp3)[1])
ratio(med_2_test_2_s, med_2_test_1_s)
ratio(med_2_test_3_s, med_2_test_2_s)
ratio(med_2_test_3_s, med_2_test_1_s)

# all of these are returning true
min_3_test_1_s = minimum(@benchmark decide_sheaf_tree_shape(skeletalColoring(3), my_decomp1)[1])
min_3_test_2_s = minimum(@benchmark decide_sheaf_tree_shape(skeletalColoring(3), my_decomp2)[1])
min_3_test_3_s = minimum(@benchmark decide_sheaf_tree_shape(skeletalColoring(3), my_decomp3)[1])
ratio(min_3_test_2_s, min_3_test_1_s)
ratio(min_3_test_3_s, min_3_test_2_s)
ratio(min_3_test_3_s, min_3_test_1_s)

med_3_test_1_s = median(@benchmark decide_sheaf_tree_shape(skeletalColoring(3), my_decomp1)[1])
med_3_test_2_s = median(@benchmark decide_sheaf_tree_shape(skeletalColoring(3), my_decomp2)[1])
med_3_test_3_s = median(@benchmark decide_sheaf_tree_shape(skeletalColoring(3), my_decomp3)[1])
ratio(med_3_test_2_s, med_3_test_1_s)
ratio(med_3_test_3_s, med_3_test_2_s)
ratio(med_3_test_3_s, med_3_test_1_s)

# all of these are returning true
min_4_test_1_s = minimum(@benchmark decide_sheaf_tree_shape(skeletalColoring(4), my_decomp1)[1])
min_4_test_2_s = minimum(@benchmark decide_sheaf_tree_shape(skeletalColoring(4), my_decomp2)[1])
min_4_test_3_s = minimum(@benchmark decide_sheaf_tree_shape(skeletalColoring(4), my_decomp3)[1])
ratio(min_4_test_2_s, min_4_test_1_s)
ratio(min_4_test_3_s, min_4_test_2_s)
ratio(min_4_test_3_s, min_4_test_1_s)

med_4_test_1_s = median(@benchmark decide_sheaf_tree_shape(skeletalColoring(4), my_decomp1)[1])
med_4_test_2_s = median(@benchmark decide_sheaf_tree_shape(skeletalColoring(4), my_decomp2)[1])
med_4_test_3_s = median(@benchmark decide_sheaf_tree_shape(skeletalColoring(4), my_decomp3)[1])
ratio(med_4_test_2_s, med_4_test_1_s)
ratio(med_4_test_3_s, med_4_test_2_s)
ratio(med_4_test_3_s, med_4_test_1_s)

# Benchmark 4 (small n per bag(5) with medium bags k(4), 3 coloring)

# Benchmark 5 (medium n per bag(10) with medium bags k(4), 3 coloring)

# Benchmark 6 (large n per bag(20) with medium bags k(4), 3 coloring)

# Benchmark 7 (small n per bag(5) with large bags k(10), 3 coloring)

# Benchmark 8 (medium n per bag(10) with large bags k(10), 3 coloring)

# Benchmark 9 (large n per bag(20) with large bags k(10), 3 coloring)