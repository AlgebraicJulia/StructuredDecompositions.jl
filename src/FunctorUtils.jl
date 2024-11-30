module FunctorUtils

export vs, skeleton

using ..Decompositions

using Catlab
using Catlab.CategoricalAlgebra
using Catlab.Graphs

#forgetful functor vs: Gr → Set taking G to VG
function vs(X::Graph) FinSet(length(vertices(X))) end
function vs(f::ACSetTransformation) components(f)[1] end

#define the functor skeleton : Set → Skel(Set)
function skeleton(s::FinSet) FinSet(length(s)) end
function skeleton(f::FinFunction)
  (dd, cc) = (dom(f), codom(f))
  #(skel_dom, skel_cod) = (skeleton(dd), skeleton(cc))
  ℓ = isempty(dd) ? Int[] : [findfirst(item -> item == f(x), collect(cc)) for x ∈ collect(dd)]
  FinFunction(ℓ, skeleton(dd), skeleton(cc))
end

end
