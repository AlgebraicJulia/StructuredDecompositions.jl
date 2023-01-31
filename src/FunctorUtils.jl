module FunctorUtils

export vs, skeleton, restriction, ↓

using ..Decompositions
using ..DecidingSheaves

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



#=
"""Given a function f: a → b, compute its image in b.
"""
function image(f::FinFunction)::FinSet
  FinSet( map( x -> f(x) , collect(dom(f))) )
end
=#
#=
"""restrict a function f: a → b to the image of a function s: a' → a (s need not be monic...).
"""
function (restriction_to)(f::FinFunction, s::FinFunction)::FinFunction
  (dom(f)== codom(s)) ? compose(legs(image(s))[1], f) : error("Domain Error: can only restrict a function to a subset of its domain")
end


(↓)(f,g) = restriction(f,g)

=#

end