module TestFunctorUtils

using Test
using PartialFunctions
using MLStyle

using ..FunctorUtils

using Catlab
using Catlab.CategoricalAlgebra

isempty(FinSet(0))

unique_initial = FinFunction(Int64[], FinSet(0), FinSet(4))
unique_initial

skeleton(unique_initial)

evens(i, j) = j ≥ 2i ? FinFunction(n -> 2n, FinSet(i), FinSet(j)) : error("try again")

f = evens(30, 60)
g = evens(10, 30)

#idempotence
@test skeleton(skeleton(f)) == skeleton(f)

t = (g,f)
compose(t...)

#=
quadruplesList(i) = map(x -> 4x, collect(FinSet(i)))
quadruplesList(10)

@test [f_restricted_to_g(x) for x ∈ collect(dom(f_restricted_to_g))] == quadruplesList(10)

@test skeleton(f_restricted_to_g) == FinFunction(quadruplesList(10), FinSet(10), FinSet(60))
=#
end
