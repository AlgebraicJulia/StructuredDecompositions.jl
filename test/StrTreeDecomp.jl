using Test

@testset "Tree Generation" begin
    @test prufer_gen(3) = [(1,), (2,), (3,)]
end

@testset "Graph Generation" begin
    @test true
end

