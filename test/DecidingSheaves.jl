module TestDecidingSheaves


using Catlab
using MLStyle
using PartialFunctions
using StructuredDecompositions.Decompositions
using StructuredDecompositions.DecidingSheaves
using StructuredDecompositions.FunctorUtils
using StructuredDecompositions.JunctionTrees
using Test


import TreeWidthSolver


# fixing bug upstream
function Catlab.WiringDiagramAlgebras.make_homomorphism(row::AbstractVector{T}, X::StructACSet{S}, Y::StructACSet{S}) where {T, S}
  components = let i = 0
    NamedTuple{ob(S)}(T[row[i+=1] for _ in parts(X,c)] for c in ob(S))
  end
  ACSetTransformation(components, X, Y)
end


function K(n::Integer)
    complete_graph(Graph, n)
end


struct Coloring
    n::Int
end


function (coloring::Coloring)(graph::Graph)
    FinSet(homomorphisms(graph, K(coloring.n); alg=HomomorphismQuery()))
end


function (coloring::Coloring)(f::ACSetTransformation)
    FinFunction(λ -> compose(f, λ), coloring(codom(f)), coloring(dom(f)))
end


function skeletal_coloring(n::Integer)
    skeleton ∘ Coloring(n)
end


function test_colorability(n::Integer, decomp::StrDecomp)
    left = is_homomorphic(ob(colimit(decomp)), K(n))
    right = first(decide_sheaf_tree_shape(skeletal_coloring(n), decomp))
    isequal(left, right)
end

@testset "Test 1" begin
    ############################
    #     EXAMPLE INSTANCE 1 str decomp
    ############################
    # 5 node cycle

    # bag 1
    H1 = @acset Graph begin
        V = 3
        E = 2
        src = [1, 2]
        tgt = [2, 3]
    end

    # adhesion 1, 2
    H12 = @acset Graph begin
        V = 2
    end

    # bag 2
    H2 = @acset Graph begin
        V = 4
        E = 3
        src = [1, 2, 3]
        tgt = [2, 3, 4]
    end

    G = @acset Graph begin
        V = 2
        E = 1
        src = [1]
        tgt = [2]
    end

    # transformations
    Γ = FinDomFunctor(
        Dict(1 => H1, 2 => H2, 3 => H12),
        Dict(
            1 => ACSetTransformation(H12, H1, V=[1, 3]),
            2 => ACSetTransformation(H12, H2, V=[1, 4])),
        ∫(G))

    manual = StrDecomp(G, Γ)
    automatic = StrDecomp(ob(colimit(manual)))

    # evaluates if decomp1 2 coloring is possible
    @test decide_sheaf_tree_shape(skeletal_coloring(2), manual)[1] == false
    @test decide_sheaf_tree_shape(skeletal_coloring(2), automatic)[1] == false
    @test decide_sheaf_tree_shape(skeletal_coloring(3), manual)[1] == true
    @test decide_sheaf_tree_shape(skeletal_coloring(3), automatic)[1] == true

    # evaluate possible 1 through 3 colorings
    @test all(test_colorability(n, manual) for n ∈ range(1, 3))
    @test all(test_colorability(n, automatic) for n ∈ range(1, 3))
end

@testset "Test 1 Automatic as Manual" begin
    # bag 1
    H1 = @acset Graph begin
        V = 3
        E = 2
        src = [1, 1]
        tgt = [2, 3]
    end

    # adhesion 1, 2
    H12 = @acset Graph begin
        V = 2
    end

    # bag 2
    H2 = @acset Graph begin
        V = 3
        E = 1
        src = [1]
        tgt = [3]
    end

    H23 = @acset Graph begin
        V = 2
    end

    H3 = @acset Graph begin
        V = 3
        E = 2
        src = [3, 3]
        tgt = [1, 2]
    end

    G = @acset Graph begin
        V = 3
        E = 2
        src = [2, 2]
        tgt = [1, 3]
    end

    # transformations
    Γ = FinDomFunctor(
        Dict(1 => H1, 2 => H2, 3 => H3, 4 => H12, 5 => H23),
        Dict(
            1 => ACSetTransformation(H12, H1, V=[2, 3]),
            2 => ACSetTransformation(H12, H2, V=[2, 1]),
            3 => ACSetTransformation(H23, H2, V=[2, 3]),
            4 => ACSetTransformation(H23, H3, V=[2, 1])),
        ∫(G))

    #manual = StrDecomp(G, Γ) # codomains don't match up
    automatic1 = StrDecomp(ob(colimit(manual)))

    # evaluates if decomp1 2 coloring is possible
    #@test decide_sheaf_tree_shape(skeletal_coloring(2), manual)[1] == false
    @test decide_sheaf_tree_shape(skeletal_coloring(2), automatic1)[1] == false
    #@test decide_sheaf_tree_shape(skeletal_coloring(3), manual)[1] == true
    @test decide_sheaf_tree_shape(skeletal_coloring(3), automatic1)[1] == true

    # evaluate possible 1 thorugh 3 colorings
    #@test all(test_colorability(n, manual) for n ∈ range(1, 3))
    @test all(test_colorability(n, automatic1) for n ∈ range(1, 3))
end

@testset "Test 2" begin
    ############################
    #     EXAMPLE INSTANCE 2 str decomp
    ############################
    # Triangle + one extension from a vertex

    # bag 1
    H1 = @acset Graph begin
        V = 3
        E = 3
        src = [1, 2, 3]
        tgt = [2, 3, 1]
    end

    # adhesion 1, 2
    H12 = @acset Graph begin
        V = 1
    end

    # bag 2
    H2 = @acset Graph begin
        V = 2
        E = 1
        src = [1]
        tgt = [2]
    end

    G = @acset Graph begin
        V = 2
        E = 1
        src = [1]
        tgt = [2]
    end

    # transformations
    Γ = FinDomFunctor(
        Dict(1 => H1, 2 => H2, 3 => H12),
        Dict(
            1 => ACSetTransformation(H12, H1, V=[1]),
            2 => ACSetTransformation(H12, H2, V=[1])),
        ∫(G))

    manual = StrDecomp(G, Γ)
    automatic = StrDecomp(ob(colimit(manual)))

    # evaluate if decomp2 2 and 3 colorings are possible
    @test decide_sheaf_tree_shape(skeletal_coloring(2), manual)[1] == false
    @test decide_sheaf_tree_shape(skeletal_coloring(3), manual)[1] == true
    @test decide_sheaf_tree_shape(skeletal_coloring(2), automatic)[1] == false
    @test decide_sheaf_tree_shape(skeletal_coloring(3), automatic)[1] == true

    # evaluate possible 1 through 3 colorings
    @test all(test_colorability(n, manual) for n ∈ range(1, 3))
    @test all(test_colorability(n, automatic) for n ∈ range(1, 3))
end

@testset "Test 3" begin
    ############################
    #     EXAMPLE INSTANCE 3 str decomp
    ############################
    # windmill shape
    # 4 triangles connected at single vertex

    # bag 1
    H1 = @acset Graph begin
        V = 5
        E = 6
        src = [1, 2, 1, 4, 5, 4]
        tgt = [2, 3, 3, 3, 3, 5]
    end

    # adhesion 1, 2
    H12 = @acset Graph begin
        V = 1
    end

    # bag 2
    H2 = @acset Graph begin
        V = 5
        E = 6
        src = [1, 2, 1, 4, 5, 4]
        tgt = [2, 3, 3, 3, 3, 5]
    end

    G = @acset Graph begin
        V = 2
        E = 1
        src = [1]
        tgt = [2]
    end

    # transformations
    Γ = FinDomFunctor(
        Dict(1 => H1, 2 => H2, 3 => H12),
        Dict(
            1 => ACSetTransformation(H12, H1, V=[3]),
            2 => ACSetTransformation(H12, H2, V=[3])),
        ∫(G))

    manual = StrDecomp(G, Γ)
    automatic = StrDecomp(ob(colimit(manual)))

    # evaluate if decomp3 2 and 3 colorings are possible
    @test decide_sheaf_tree_shape(skeletal_coloring(2), manual)[1] == false
    @test decide_sheaf_tree_shape(skeletal_coloring(3), manual)[1] == true
    @test decide_sheaf_tree_shape(skeletal_coloring(2), automatic)[1] == false
    @test decide_sheaf_tree_shape(skeletal_coloring(3), automatic)[1] == true

    # evaluate possible 1 through 3 colorings
    @test all(test_colorability(n, manual) for n ∈ range(1, 3))
    @test all(test_colorability(n, automatic) for n ∈ range(1, 3))
end

@testset "Test 4" begin
    ############################
    #     EXAMPLE INSTANCE 4 str decomp
    ############################
    # star inside of pentagon
    # connections at each star vertex to a pentagon vertex

    # bag 1
    H1 = @acset Graph begin
        V = 10
        E = 10
        src = [1, 2, 3, 4, 5, 1, 2, 3, 4, 5]
        tgt = [2, 3, 4, 5, 1, 6, 7, 8, 9, 10]
    end

    # adhesion 1, 2
    H12 = @acset Graph begin
        V = 5
    end

    # bag 2
    H2 = @acset Graph begin
        V = 5
        E = 5
        src = [1, 1, 2, 2, 3]
        tgt = [3, 4, 4, 5, 5]
    end

    G = @acset Graph begin
        V = 2
        E = 1
        src = [1]
        tgt = [2]
    end

    # transformations
    Γ = FinDomFunctor(
        Dict(1 => H1, 2 => H2, 3 => H12),
        Dict(
            1 => ACSetTransformation(H12, H1, V=[6, 7, 8, 9, 10]),
            2 => ACSetTransformation(H12, H2, V=[1, 2, 3, 4, 5])),
        ∫(G))

    manual = StrDecomp(G, Γ)
    automatic = StrDecomp(ob(colimit(manual)))

    # evaluate if decomp4 2 and 3 colorings are possible 
    @test decide_sheaf_tree_shape(skeletal_coloring(2), manual)[1] == false
    @test decide_sheaf_tree_shape(skeletal_coloring(3), manual)[1] == true
    @test decide_sheaf_tree_shape(skeletal_coloring(2), automatic)[1] == false
    @test decide_sheaf_tree_shape(skeletal_coloring(3), automatic)[1] == true

    # evaluate possible 1 through 3 colorings
    @test all(test_colorability(n, manual) for n ∈ range(1, 3))
    @test all(test_colorability(n, automatic) for n ∈ range(1, 3))
end

@testset "Test 5" begin
    ############################
    #     EXAMPLE INSTANCE 5 str decomp
    ############################
    # 2 by 2 square inside of 3 by 3 square
    # connections at each vertex of 2 by 2 to middle
    # vertices of 3 by 3 square

    # bag 1
    H1 = @acset Graph begin
        V = 12
        E = 12
        src = [5, 6, 7, 8, 9, 10, 11, 12, 6, 8, 10, 12]
        tgt = [6, 7, 8, 9, 10, 11, 12, 5, 1, 2, 3, 4]
    end

    # adhesion 1, 2
    H12 = @acset Graph begin
        V = 4
    end

    # bag 2
    H2 = @acset Graph begin
        V = 4
        E = 4
        src = [1, 2, 3, 4]
        tgt = [2, 3, 4, 1]
    end

    G = @acset Graph begin
        V = 2
        E = 1
        src = [1]
        tgt = [2]
    end

    # transformations
    Γ = FinDomFunctor(
        Dict(1 => H1, 2 => H2, 3 => H12),
        Dict(
            1 => ACSetTransformation(H12, H1, V=[1, 2, 3, 4]),
            2 => ACSetTransformation(H12, H2, V=[1, 2, 3, 4])),
        ∫(G))

    manual = StrDecomp(G, Γ)
    automatic = StrDecomp(ob(colimit(manual)))

    # evaluate if decomp5 2 and 3 colorings are possible
    @test decide_sheaf_tree_shape(skeletal_coloring(2), manual)[1] == false
    @test decide_sheaf_tree_shape(skeletal_coloring(3), manual)[1] == true
    @test decide_sheaf_tree_shape(skeletal_coloring(2), automatic)[1] == false
    @test decide_sheaf_tree_shape(skeletal_coloring(3), automatic)[1] == true

    # evaluate possible 1 through 3 colorings
    @test all(test_colorability(n, manual) for n ∈ range(1, 3))
    @test all(test_colorability(n, automatic) for n ∈ range(1, 3))
end

@testset "Test 6" begin
    ############################
    #     EXAMPLE INSTANCE 6 str decomp
    ############################
    # K4 inside of 3 by 3 square
    # connections at each vertex of K4 to middle
    # vertices of 3 by 3 square

    # bag 1
    H1 = @acset Graph begin
        V = 12
        E = 12
        src = [5, 6, 7, 8, 9, 10, 11, 12, 6, 8, 10, 12]
        tgt = [6, 7, 8, 9, 10, 11, 12, 5, 1, 2, 3, 4]
    end

    # adhesion 1, 2
    H12 = @acset Graph begin
        V = 4
    end

    # bag 2
    H2 = @acset Graph begin
        V = 4
        E = 6
        src = [1, 2, 3, 4, 1, 2]
        tgt = [2, 3, 4, 1, 3, 4]
    end

    G = @acset Graph begin
        V = 2
        E = 1
        src = [1]
        tgt = [2]
    end

    # transformations
    Γ = FinDomFunctor(
        Dict(1 => H1, 2 => H2, 3 => H12),
        Dict(
            1 => ACSetTransformation(H12, H1, V=[1, 2, 3, 4]),
            2 => ACSetTransformation(H12, H2, V=[1, 2, 3, 4])),
        ∫(G))

    manual = StrDecomp(G, Γ)
    automatic = StrDecomp(ob(colimit(manual)))

    # evaluate if decomp6 2 and 3 colorings are possible
    @test decide_sheaf_tree_shape(skeletal_coloring(2), manual)[1] == false
    @test decide_sheaf_tree_shape(skeletal_coloring(3), manual)[1] == false
    @test decide_sheaf_tree_shape(skeletal_coloring(2), automatic)[1] == false
    @test decide_sheaf_tree_shape(skeletal_coloring(3), automatic)[1] == false

    # evaluate possible 1 through 3 colorings
    @test all(test_colorability(n, manual) for n ∈ range(1, 3))
    @test all(test_colorability(n, automatic) for n ∈ range(1, 3))
end

@testset "Test 7" begin
    ############################
    #     EXAMPLE INSTANCE 7 str decomp
    ############################
    # slightly different variation of 6
    # connections are now to corners not center

    # bag 1
    H1 = @acset Graph begin
        V = 12
        E = 12
        src = [5, 6, 7, 8, 9, 10, 11, 12, 5, 7, 9, 11]
        tgt = [6, 7, 8, 9, 10, 11, 12, 5, 1, 2, 3, 4]
    end

    # adhesion 1, 2
    H12 = @acset Graph begin
        V = 4
    end

    # bag 2
    H2 = @acset Graph begin
        V = 4
        E = 6
        src = [1, 2, 3, 4, 1, 2]
        tgt = [2, 3, 4, 1, 3, 4]
    end

    G = @acset Graph begin
        V = 2
        E = 1
        src = [1]
        tgt = [2]
    end

    # transformations
    Γ = FinDomFunctor(
        Dict(1 => H1, 2 => H2, 3 => H12),
        Dict(
            1 => ACSetTransformation(H12, H1, V=[1, 2, 3, 4]),
            2 => ACSetTransformation(H12, H2, V=[1, 2, 3, 4])),
        ∫(G))

    manual = StrDecomp(G, Γ)
    automatic = StrDecomp(ob(colimit(manual)))

    # evaluate if decomp7 2 and 3 colorings are possible
    @test decide_sheaf_tree_shape(skeletal_coloring(2), manual)[1] == false
    @test decide_sheaf_tree_shape(skeletal_coloring(3), manual)[1] == false
    @test decide_sheaf_tree_shape(skeletal_coloring(2), automatic)[1] == false
    @test decide_sheaf_tree_shape(skeletal_coloring(3), automatic)[1] == false

    # evaluate possible 1 through 3 colorings
    @test all(test_colorability(n, manual) for n ∈ range(1, 3))
    @test all(test_colorability(n, automatic) for n ∈ range(1, 3))
end

@testset "Test 8" begin
    ############################
    #     EXAMPLE INSTANCE 8 str decomp
    ############################
    # large box with incomplete 3x3 lattice inside
    # incomplete lattice missing a corner

    # bag 1
    H1 = @acset Graph begin
        V = 9
        E = 9
        src = [1, 2, 3, 4, 1, 2, 2, 3, 4]
        tgt = [2, 3, 4, 1, 5, 6, 7, 8, 9]
    end

    # adhesion 1, 2
    H12 = @acset Graph begin
        V = 5
    end

    # bag 2
    H2 = @acset Graph begin
        V = 8
        E = 10
        src = [1, 1, 2, 3, 3, 4, 4, 5, 6, 7]
        tgt = [2, 3, 4, 4, 6, 5, 7, 8, 7, 8]
    end

    G = @acset Graph begin
        V = 2
        E = 1
        src = [1]
        tgt = [2]
    end

    # transformations
    Γ = FinDomFunctor(
        Dict(1 => H1, 2 => H2, 3 => H12),
        Dict(
            1 => ACSetTransformation(H12, H1, V=[5, 6, 7, 8, 9]),
            2 => ACSetTransformation(H12, H2, V=[1, 2, 5, 8, 6])),
        ∫(G))

    manual = StrDecomp(G, Γ)
    automatic = StrDecomp(ob(colimit(manual)))

    # evaluate if decomp8 2 and 3 colorings are possible
    @test decide_sheaf_tree_shape(skeletal_coloring(2), manual)[1] == false
    @test decide_sheaf_tree_shape(skeletal_coloring(3), manual)[1] == true
    @test decide_sheaf_tree_shape(skeletal_coloring(2), automatic)[1] == false
    @test decide_sheaf_tree_shape(skeletal_coloring(3), automatic)[1] == true

    # evaluate possible 1 through 3 colorings
    @test all(test_colorability(n, manual) for n ∈ range(1, 3))
    @test all(test_colorability(n, automatic) for n ∈ range(1, 3))

    graph = ob(colimit(manual))

    #print(ob(colimit(automatic)))

    automatic = StrDecomp(graph; alg=BT())

    #print(automatic)
end

@testset "Test 9" begin
    ############################
    #     EXAMPLE INSTANCE 9 str decomp
    ############################
    # circle (defined by 3 vertices) with inscribed triangle 1
    # triangle 1 inscribed with triangle 2
    # triangle 2 inscribed with upsidedown triangle 3

    # bag 1
    H1 = @acset Graph begin
        V = 12
        E = 15
        src = [1, 1, 2, 3, 3, 4, 5, 5, 6, 1, 2, 3, 4, 5, 6]
        tgt = [2, 3, 3, 4, 5, 5, 6, 1, 1, 7, 8, 9, 10, 11, 12]
    end

    # adhesion 1, 2
    H12 = @acset Graph begin
        V = 6
    end

    # bag 2
    H2 = @acset Graph begin
        V = 6
        E = 9
        src = [1, 1, 2, 2, 2, 3, 4, 4, 5]
        tgt = [2, 6, 3, 4, 6, 4, 5, 6, 6]
    end

    G = @acset Graph begin
        V = 2
        E = 1
        src = [1]
        tgt = [2]
    end

    # transformations
    Γ = FinDomFunctor(
        Dict(1 => H1, 2 => H2, 3 => H12),
        Dict(
            1 => ACSetTransformation(H12, H1, V=[7, 8, 9, 10, 11, 12]),
            2 => ACSetTransformation(H12, H2, V=[1, 2, 3, 4, 5, 6])),
        ∫(G))

    manual = StrDecomp(G, Γ)
    automatic = StrDecomp(ob(colimit(manual)))

    # evaluate if decomp9 2 and 3 colorings are possible
    @test decide_sheaf_tree_shape(skeletal_coloring(2), manual)[1] == false
    @test decide_sheaf_tree_shape(skeletal_coloring(3), manual)[1] == true
    @test decide_sheaf_tree_shape(skeletal_coloring(2), automatic)[1] == false
    @test decide_sheaf_tree_shape(skeletal_coloring(3), automatic)[1] == true

    # evaluate possible 1 through 3 colorings
    @test all(test_colorability(n, manual) for n ∈ range(1, 3))
    @test all(test_colorability(n, automatic) for n ∈ range(1, 3))
end

end
