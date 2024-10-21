using Test

# @testset "Decompositions" begin
#   include("Decompositions.jl")
# end

@testset "DecidingSheaves" begin
  include("DecidingSheaves.jl")
end

@testset "FunctorUtils" begin
  include("FunctorUtils.jl")
end

@testset "JunctionTrees" begin
  include("JunctionTrees.jl")
end

@testset "NestedUWDs" begin
  include("NestedUWDs.jl")
end
