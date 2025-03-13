using Catlab
using StructuredDecompositions
using Test

# fix bug upstream
function Catlab.WiringDiagramAlgebras.make_homomorphism(row::AbstractVector{T}, X::StructACSet{S}, Y::StructACSet{S}) where {T, S}
    components = let i = 0
        NamedTuple{ob(S)}(T[row[i+=1] for _ in parts(X,c)] for c in ob(S))
    end

    ACSetTransformation(components, X, Y)
end

@testset "coloring" begin
    function K(n::Integer)
        complete_graph(SymmetricGraph, n)
    end

    struct Coloring
        n::Int
    end

    function (coloring::Coloring)(graph::SymmetricGraph)
        FinSet(homomorphisms(graph, K(coloring.n); alg=HomomorphismQuery()))
    end

    function (coloring::Coloring)(f::ACSetTransformation)
        FinFunction(λ -> compose(f, λ), coloring(codom(f)), coloring(dom(f)))
    end

    function skeletal_coloring(n::Integer)
        skeleton ∘ Coloring(n)
    end

    # Peterson graph
    graph = SymmetricGraph(10)

    add_edges!(graph,
        [1, 1, 1, 2, 2, 3, 3, 4, 4, 5, 6, 6, 7, 8, 9 ],
        [2, 5, 9, 3, 7, 10,4, 5, 8, 6, 7, 10,8, 9, 10])

    decomposition = StrDecomp(graph)
    @test !first(decide_sheaf_tree_shape(skeletal_coloring(1), decomposition))
    @test !first(decide_sheaf_tree_shape(skeletal_coloring(2), decomposition))
    @test  first(decide_sheaf_tree_shape(skeletal_coloring(3), decomposition)) 
    @test  first(decide_sheaf_tree_shape(skeletal_coloring(4), decomposition))

    # path graph
    graph = path_graph(SymmetricGraph, 10)

    decomposition = StrDecomp(graph)
    @test !first(decide_sheaf_tree_shape(skeletal_coloring(1), decomposition))
    @test  first(decide_sheaf_tree_shape(skeletal_coloring(2), decomposition))
    @test  first(decide_sheaf_tree_shape(skeletal_coloring(3), decomposition)) 
    @test  first(decide_sheaf_tree_shape(skeletal_coloring(4), decomposition))

    # cycle graph
    graph = cycle_graph(SymmetricGraph, 10)

    decomposition = StrDecomp(graph)
    @test !first(decide_sheaf_tree_shape(skeletal_coloring(1), decomposition))
    @test  first(decide_sheaf_tree_shape(skeletal_coloring(2), decomposition))
    @test  first(decide_sheaf_tree_shape(skeletal_coloring(3), decomposition)) 
    @test  first(decide_sheaf_tree_shape(skeletal_coloring(4), decomposition))

    # star graph
    graph = star_graph(SymmetricGraph, 10)

    decomposition = StrDecomp(graph)
    @test !first(decide_sheaf_tree_shape(skeletal_coloring(1), decomposition))
    @test  first(decide_sheaf_tree_shape(skeletal_coloring(2), decomposition))
    @test  first(decide_sheaf_tree_shape(skeletal_coloring(3), decomposition)) 
    @test  first(decide_sheaf_tree_shape(skeletal_coloring(4), decomposition))

    # wheel graph
    graph = wheel_graph(SymmetricGraph, 10)

    decomposition = StrDecomp(graph)
    @test !first(decide_sheaf_tree_shape(skeletal_coloring(1), decomposition))
    @test !first(decide_sheaf_tree_shape(skeletal_coloring(2), decomposition))
    @test !first(decide_sheaf_tree_shape(skeletal_coloring(3), decomposition)) 
    @test  first(decide_sheaf_tree_shape(skeletal_coloring(4), decomposition))
end

