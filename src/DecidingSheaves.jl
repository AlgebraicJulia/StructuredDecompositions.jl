module DecidingSheaves

export Presheaf, Sheaf, decide, decide_tree_shape

using ..Decompositions

#####################
#   COMPUTATION
#####################
"""(Pre)-Sheaves
"""
abstract type Presheaf end

abstract type Sheaf <: Presheaf end

"""
Filtering algorithm. 
Note: we are assuming that we only know how to work with FinSet(Int) !

INPUT: a Finset^{op}-valued structured decomposition d : FG → Span Finset^{op} 
          (which is expected to be in co-decomposition form; 
            i.e. as a diagram d : FG → Cospan Finset )
       and an indexed span ( (ℓ, r), ( d(ℓ), d(r) ) ) in d 
          (i.e a pair consisting of span (ℓ, r) in ∫G and its image under d)

OUTPUT: a structured decomposition obtained by replacing the span de in d 
        by the span obtained by projecting the pullback of de (i.e. taking images)
"""
function adhesion_filter(tup, d::StructuredDecomposition)
  # d_csp is the cospan dx₁ -> de <- dx₂ corresp to some edge e = x₁x₂ in shape(d)
  (csp, d_csp)      = tup  #unpack the tuple
  # the pullback cone dx₁ <-l₁-- p --l₂ --> dx₂ with legs l₁ and l₂
  p_cone            = pullback(d_csp)
  p, p_legs         = ob(p_cone), legs(p_cone)
  #for each leg lᵢ : p → xᵢ of the pullback cone, 
  #compute its image ιᵢ : im lᵢ → dxᵢ
  imgs              = map( f -> legs(image(f))[1], p_legs)
  #now get the new desired cospan; i.e.  im l₁ --ι₁--> dx₁ --l₁--> de <--l₂--dx₂ <--ι₂-- im l₂
  #note, you need to take the skeleton for good measure
  new_d_csp         = map(skeleton ∘ compose, zip(imgs, p_legs))
  #get the domain of d 
  d_dom             = dom(d.diagram)
  #now make the new decomposition, call it δ
  #start with the object map δ₀
  ob_replace(i, x)  = x == dom(d, csp[i]) ? dom(new_d_csp[i]) : ob_map(d,x) 
  ob_replace(x)     = (ob_replace $ 1) ∘ (ob_replace $ 2)
  δ₀                = Dict( x => ob_replace(x) for x ∈ ob_generators(d_dom))
  #now do the same thing with the morphism map
  mor_replace(i, f) = f == csp[i] ? new_d_csp[i] : hom_map(d,f) 
  mor_replace(f)    = (mor_replace $ 1) ∘ (mor_replace $ 2)
  δ₁                = Dict( f => mor_replace(f) for f ∈ hom_generators(d_dom))
  StrDecom(d.decomp_shape, d.domain, FinDomFunctor(δ₀, δ₁, d.domain), d.decomp_type)
end

#given a decomposition d and a list ℓ of indexed_spans of the decomposition,
#run adhesion_filter for each indexed_span in ℓ
function adhesion_filter(ℓ, d::StructuredDecomposition) foldr( ∘, reverse([adhesion_filter $ ll for ll ∈ ℓ])) end

"""Solve the decision problem encoded by a sheaf. 
The algorithm is as follows: 
  compute on each bag,
  compute composites on edges, 
  project back down to bags
  answer 
    "no" if there is an empty bag; 
    "yes" otherwise.
"""
function decide_tree_shape(f::Sheaf, d::StructuredDecomposition)::Bool
  d_filtered = adhesion_filter(adhesionSpans(d, true), d)
  foldr(&, map( !isempty, bags(d_filtered)))
end

end