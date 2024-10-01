# [Functor Utils] (@id FunctorUtils)

Functor Utils only includes 4 functions and builds off of Decompositions.

We first define the forgetful functors vs.

```julia
function vs(X::Graph) FinSet(length(vertices(X))) end
function vs(f::ACSetTransformation) components(f)[1] end
```

We also define the functor skeleton taking set to the skeleton of the set.

```julia
function skeleton(s::FinSet) FinSet(length(s)) end
function skeleton(f::FinFunction)
  (dd, cc) = (dom(f), codom(f))
  ℓ = isempty(dd) ? Int[] : [findfirst(item -> item == f(x), collect(cc)) for x ∈ collect(dd)]
  FinFunction(ℓ, skeleton(dd), skeleton(cc))
end
```