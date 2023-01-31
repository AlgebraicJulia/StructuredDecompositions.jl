module TestFunctorUtils

using Test
using PartialFunctions
using MLStyle

using StructuredDecompositions.FunctorUtils

using Catlab
using Catlab.CategoricalAlgebra

isempty(FinSet(0))

unique_initial = FinFunction(Int64[], FinSet(0), FinSet(4))
unique_initial

skeleton(unique_initial)

evens(i, j) = j ≥ 2i ? FinFunction(n -> 2n, FinSet(i), FinSet(j)) : error(i, " is greater than half of ", j)

f = evens(30, 60)
g = evens(10, 30)

#idempotence
@test skeleton(skeleton(f)) == skeleton(f)

end
