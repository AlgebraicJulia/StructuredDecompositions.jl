module TestDecidingSheaves

using Test
using PartialFunctions
using MLStyle

using StructuredDecompositions.Decompositions
using StructuredDecompositions.DecidingSheaves
using StructuredDecompositions.FunctorUtils

using Catlab.Graphs
using Catlab.ACSetInterface
using Catlab.CategoricalAlgebra
using Catlab.Graphics

"""
An example: graph colorings
"""
#an H-coloring is a hom onto H
struct Coloring
  n     #the target graph
  func  #the function mapping opens to lists of homs from G to K_n
end

#construct an n-coloring
K(n) = complete_graph(Graph, n)
Coloring(n) = Coloring(n, g -> homomorphisms(g, K(n) ))
#make it callable
(c::Coloring)(X::Graph) = FinSet(c.func(X)) # given graph homos #f: G₁ → G₂ get morphism col(G₂) → col(G₁) by precomposition: take each λ₂ ∈ col(G₂) to hf ∈ col(G)
function (c::Coloring)(f::ACSetTransformation)  
  (G₁, G₂)   = (dom(f), codom(f)) 
  (cG₁, cG₂) = (c(G₁), c(G₂))
  FinFunction( λ₂ -> compose(f,λ₂), cG₂, cG₁ ) #note the contravariance
end

skeletalColoring(n) = skeleton ∘ Coloring(n)

colorability_test(n, the_test_case) = is_homomorphic(ob(colimit(the_test_case)), K(n)) == decide_sheaf_tree_shape(skeletalColoring(n), the_test_case)[1]

# TODO: Remove me?
# is_homomorphic(ob(colimit(my_decomp1)), K(2))

@testset "Test 1" begin

  ############################
  #     EXAMPLE INSTANCE 1 str decomp
  ############################
  # 7 node cycle

  #bag 1
  H₁ = @acset Graph begin
    V = 3
    E = 2
    src = [1, 2]
    tgt = [2, 3]
  end

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

  #transformations
  Γₛ⁰ = Dict(1 => H₁, 2 => H₂, 3 => H₁₂)
  Γₛ = FinDomFunctor(
    Γₛ⁰,
    Dict(
      1 => ACSetTransformation(Γₛ⁰[3], Γₛ⁰[1], V=[1, 3]),
      2 => ACSetTransformation(Γₛ⁰[3], Γₛ⁰[2], V=[4, 1]),
    ),
    ∫(Gₛ)
  )

  my_decomp1  = StrDecomp(Gₛ, Γₛ)

  @show "Decomp1"

  @test decide_sheaf_tree_shape(skeletalColoring(2), my_decomp1)[1] == false
  #evaluates if decomp1 2 coloring is possible

  @show "Decomp1 colorability_test"

  @test all(colorability_test(n, my_decomp1) for n ∈ range(1,3))
  #evaluate possible 1 thorugh 3 colorings

end

@testset "Test 2" begin

  ############################
  #     EXAMPLE INSTANCE 2 str decomp
  ############################
  # Triangle + one extension from a vertex

  #bag 1
  H₁ = @acset Graph begin
    V = 3
    E = 3
    src = [1, 2, 3]
    tgt = [2, 3, 1]
  end

  #adhesion 1,2
  H₁₂ = @acset Graph begin
    V = 2
  end

  #bag 2
  H₂ = @acset Graph begin
    V = 2
    E = 1
    src = [1]
    tgt = [2]
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
     1 => ACSetTransformation(Γₛ⁰[3], Γₛ⁰[1], V=[1, 3]),
     2 => ACSetTransformation(Γₛ⁰[3], Γₛ⁰[2], V=[2, 1]),
   ),
   ∫(Gₛ)
  )

  my_decomp2  = StrDecomp(Gₛ, Γₛ)

  @show "Decomp2"

  @test decide_sheaf_tree_shape(skeletalColoring(2), my_decomp2)[1] == false
  @test decide_sheaf_tree_shape(skeletalColoring(3), my_decomp2)[1] == true
  #evaluate if decomp2 2 and 3 colorings are possible

  @show "Decomp2 colorability_test"

  @test all(colorability_test(n, my_decomp2) for n ∈ range(1, 3))
  #evaluate possible 1 through 3 colorings

end

@testset "Test 3" begin

  ############################
  #     EXAMPLE INSTANCE 3 str decomp
  ############################
  # windmill shape
  # 4 triangles connected at single vertex

  #bag 1
  H₁ = @acset Graph begin
    V = 5
    E = 6
    src = [1, 2, 1, 4, 5, 4]
    tgt = [2, 3, 3, 3, 3, 5]
  end

  #adhesion 1,2
  H₁₂ = @acset Graph begin
    V = 1
  end

  #bag 2
  H₂ = @acset Graph begin
    V = 5
    E = 6
    src = [1, 2, 1, 4, 5, 4]
    tgt = [2, 3, 3, 3, 3, 5]
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
      1 => ACSetTransformation(Γₛ⁰[3], Γₛ⁰[1], V=[3]),
      2 => ACSetTransformation(Γₛ⁰[3], Γₛ⁰[2], V=[3]),
    ),
    ∫(Gₛ)
  )

  my_decomp3  = StrDecomp(Gₛ, Γₛ)
  
  @show "Decomp3"

  @test decide_sheaf_tree_shape(skeletalColoring(2), my_decomp3)[1] == false
  @test decide_sheaf_tree_shape(skeletalColoring(3), my_decomp3)[1] == true
  #evaluate if decomp3 2 and 3 colorings are possible

  @show "Decomp3 colorability_test"

  @test all(colorability_test(n, my_decomp3) for n ∈ range(1, 3))
  #evaluate possible 1 through 3 colorings

end

@testset "Test 4" begin

  ############################
  #     EXAMPLE INSTANCE 4 str decomp
  ############################
  # star inside of pentagon
  # connections at each star vertex to a pentagon vertex

  #bag 1
  H₁ = @acset Graph begin
    V = 10
    E = 10
    src = [1, 2, 3, 4, 5, 1, 2, 3, 4, 5]
    tgt = [2, 3, 4, 5, 1, 6, 7, 8, 9, 10]
  end

  #adhesion 1,2
  H₁₂ = @acset Graph begin
    V = 5
  end

  #bag 2
  H₂ = @acset Graph begin
    V = 5
    E = 5
    src = [1, 1, 2, 2, 3]
    tgt = [3, 4, 4, 5, 5]
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
      1 => ACSetTransformation(Γₛ⁰[3], Γₛ⁰[1], V=[6, 7, 8, 9, 10]),
      2 => ACSetTransformation(Γₛ⁰[3], Γₛ⁰[2], V=[1, 2, 3, 4, 5]),
    ),
    ∫(Gₛ)
  )

  my_decomp4  = StrDecomp(Gₛ, Γₛ)

  @show "Decomp4"

  @test decide_sheaf_tree_shape(skeletalColoring(2), my_decomp4)[1] == false
  @test decide_sheaf_tree_shape(skeletalColoring(3), my_decomp4)[1] == true
  #evaluate if decomp4 2 and 3 colorings are possible 

  @show "Decomp4 colorability_test"

  @test all(colorability_test(n, my_decomp4) for n ∈ range(1, 3))
  #evaluate possible 1 through 3 colorings

end

@testset "Test 5" begin

  ############################
  #     EXAMPLE INSTANCE 5 str decomp
  ############################
  # 2 by 2 square inside of 3 by 3 square
  # connections at each vertex of 2 by 2 to middle
  # vertices of 3 by 3 square

  #bag 1
  H₁ = @acset Graph begin
    V = 12
    E = 12
    src = [5, 6, 7, 8, 9, 10, 11, 12, 6, 8, 10, 12]
    tgt = [6, 7, 8, 9, 10, 11, 12, 5, 1, 2, 3, 4]
  end

  #adhesion 1,2
  H₁₂ = @acset Graph begin
    V = 4
  end

  #bag 2
  H₂ = @acset Graph begin
    V = 4
    E = 4
    src = [1, 2, 3, 4]
    tgt = [2, 3, 4, 1]
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
      1 => ACSetTransformation(Γₛ⁰[3], Γₛ⁰[1], V=[1, 2, 3, 4]),
      2 => ACSetTransformation(Γₛ⁰[3], Γₛ⁰[2], V=[1, 2, 3, 4]),
    ),
    ∫(Gₛ)
  )

  my_decomp5  = StrDecomp(Gₛ, Γₛ)

  @show "Decomp5"

  @test decide_sheaf_tree_shape(skeletalColoring(2), my_decomp5)[1] == false
  @test decide_sheaf_tree_shape(skeletalColoring(3), my_decomp5)[1] == true
  #evaluate if decomp5 2 and 3 colorings are possible

  @show "Decomp5 colorability_test"

  @test all(colorability_test(n, my_decomp5) for n ∈ range(1, 3))
  #evaluate possible 1 through 3 colorings

end

@testset "Test 6" begin

  ############################
  #     EXAMPLE INSTANCE 6 str decomp
  ############################
  # K4 inside of 3 by 3 square
  # connections at each vertex of K4 to middle
  # vertices of 3 by 3 square

  #bag 1
  H₁ = @acset Graph begin
    V = 12
    E = 12
    src = [5, 6, 7, 8, 9, 10, 11, 12, 6, 8, 10, 12]
    tgt = [6, 7, 8, 9, 10, 11, 12, 5, 1, 2, 3, 4]
  end

  #adhesion 1,2
  H₁₂ = @acset Graph begin
    V = 4
  end

  #bag 2
  H₂ = @acset Graph begin
    V = 4
    E = 6
    src = [1, 2, 3, 4, 1, 2]
    tgt = [2, 3, 4, 1, 3, 4]
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
      1 => ACSetTransformation(Γₛ⁰[3], Γₛ⁰[1], V=[1, 2, 3, 4]),
      2 => ACSetTransformation(Γₛ⁰[3], Γₛ⁰[2], V=[1, 2, 3, 4]),
    ),
    ∫(Gₛ)
  )

  my_decomp6  = StrDecomp(Gₛ, Γₛ)

  @show "Decomp6"

  @test decide_sheaf_tree_shape(skeletalColoring(2), my_decomp6)[1] == false
  @test decide_sheaf_tree_shape(skeletalColoring(3), my_decomp6)[1] == false
  #evaluate if decomp6 2 and 3 colorings are possible

  @show "Decomp6 colorability_test"

  @test all(colorability_test(n, my_decomp6) for n ∈ range(1, 3))
  #evaluate possible 1 through 3 colorings
end

@testset "Test 7" begin

  ############################
  #     EXAMPLE INSTANCE 7 str decomp
  ############################
  # slightly different variation of 6
  # connections are now to corners not center

  #bag 1
  H₁ = @acset Graph begin
    V = 12
    E = 12
    src = [5, 6, 7, 8, 9, 10, 11, 12, 5, 7, 9, 11]
    tgt = [6, 7, 8, 9, 10, 11, 12, 5, 1, 2, 3, 4]
  end

  #adhesion 1,2
  H₁₂ = @acset Graph begin
    V = 4
  end

  #bag 2
  H₂ = @acset Graph begin
    V = 4
    E = 6
    src = [1, 2, 3, 4, 1, 2]
    tgt = [2, 3, 4, 1, 3, 4]
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
      1 => ACSetTransformation(Γₛ⁰[3], Γₛ⁰[1], V=[1, 2, 3, 4]),
      2 => ACSetTransformation(Γₛ⁰[3], Γₛ⁰[2], V=[1, 2, 3, 4]),
    ),
    ∫(Gₛ)
  )

  my_decomp7  = StrDecomp(Gₛ, Γₛ)

  @show "Decomp7"

  @test decide_sheaf_tree_shape(skeletalColoring(2), my_decomp7)[1] == false
  @test decide_sheaf_tree_shape(skeletalColoring(3), my_decomp7)[1] == false
  #evaluate if decomp7 2 and 3 colorings are possible

  @show "Decomp7 colorability_test"

  @test all(colorability_test(n, my_decomp7) for n ∈ range(1, 3))
  #evaluate possible 1 through 3 colorings

end

@testset "Test 8" begin

  ############################
  #     EXAMPLE INSTANCE 8 str decomp
  ############################
  # large box with incomplete 3x3 lattice inside
  # incomplete lattice missing a corner

  #bag 1
  H₁ = @acset Graph begin
    V = 11
    E = 11
    src = [1, 2, 3, 4, 5, 1, 2, 2, 3, 4, 5]
    tgt = [2, 3, 4, 5, 1, 6, 7, 8, 9, 10, 11]
  end

  #adhesion 1,2
  H₁₂ = @acset Graph begin
    V = 6
  end

  #bag 2
  H₂ = @acset Graph begin
    V = 8
    E = 10
    src = [1, 1, 2, 3, 3, 4, 4, 5, 6, 7]
    tgt = [2, 3, 4, 4, 6, 5, 7, 8, 7, 8]
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
      1 => ACSetTransformation(Γₛ⁰[3], Γₛ⁰[1], V=[6, 7, 8, 9, 10, 11]),
      2 => ACSetTransformation(Γₛ⁰[3], Γₛ⁰[2], V=[1, 2, 5, 8, 7, 6]),
    ),
    ∫(Gₛ)
  )

  my_decomp8  = StrDecomp(Gₛ, Γₛ)

  @show "Decomp8"

  @test decide_sheaf_tree_shape(skeletalColoring(2), my_decomp8)[1] == false
  @test decide_sheaf_tree_shape(skeletalColoring(3), my_decomp8)[1] == true
  #evaluate if decomp8 2 and 3 colorings are possible

  @show "Decomp8 colorability_test"

  @test all(colorability_test(n, my_decomp8) for n ∈ range(1, 3))
  #evaluate possible 1 through 3 colorings

end

@testset "Test 9" begin

  ############################
  #     EXAMPLE INSTANCE 9 str decomp
  ############################
  # circle (defined by 3 vertices) with inscribed triangle 1
  # triangle 1 inscribed with triangle 2
  # triangle 2 inscribed with upsidedown triangle 3


  #bag 1
  H₁ = @acset Graph begin
    V = 12
    E = 15
    src = [1, 1, 2, 3, 3, 4, 5, 5, 6, 1, 2, 3, 4, 5, 6]
    tgt = [2, 3, 3, 4, 5, 5, 6, 1, 1, 7, 8, 9, 10, 11, 12]
  end

  #adhesion 1,2
  H₁₂ = @acset Graph begin
    V = 6
  end

  #bag 2
  H₂ = @acset Graph begin
    V = 6
    E = 9
    src = [1, 1, 2, 2, 2, 3, 4, 4, 5]
    tgt = [2, 6, 3, 4, 6, 4, 5, 6, 6]
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
      1 => ACSetTransformation(Γₛ⁰[3], Γₛ⁰[1], V=[7, 8, 9, 10, 11, 12]),
      2 => ACSetTransformation(Γₛ⁰[3], Γₛ⁰[2], V=[1, 2, 3, 4, 5, 6]),
    ),
    ∫(Gₛ)
  )

  my_decomp9  = StrDecomp(Gₛ, Γₛ)

  @show "Decomp9"

  @test decide_sheaf_tree_shape(skeletalColoring(2), my_decomp9)[1] == false
  @test decide_sheaf_tree_shape(skeletalColoring(3), my_decomp9)[1] == true
  #evaluate if decomp9 2 and 3 colorings are possible

  @show "Decomp9 colorability_test"

  @test all(colorability_test(n, my_decomp9) for n ∈ range(1, 3))
  #evaluate possible 1 through 3 colorings

end

end