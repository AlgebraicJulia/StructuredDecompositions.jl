module FunctorUtils

export vs, skeleton

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
  ℓ = [findfirst(item -> item == f(x), collect(cc)) for x ∈ collect(dd)]
  FinFunction(ℓ, skeleton(dd), skeleton(cc))
  #=
  #make function from dictionary
  zip_iso(xs, ys) = begin
    pairs = zip(collect(xs), collect(ys))
    my_dict = Dict(pairs)
    x -> my_dict[x]
  end 
  left_iso  = FinFunction(zip_iso(skel_dom, dd), skel_dom, dd)
  right_iso = FinFunction(zip_iso(cc, skel_cod), cc, skel_cod)
  left_iso ⋅ f ⋅ right_iso
  =#
end

end