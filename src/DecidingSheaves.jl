module DecidingSheaves

export Presheaf, Sheaf, decide_sheaf_tree_shape

using ..Decompositions
using ..FunctorUtils

using Catlab
using Catlab.CategoricalAlgebra

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
function adhesion_filter(tup::Tuple, d::StructuredDecomposition)
  if d.decomp_type == Decomposition
    error("expecting ", CoDecomposition, " given ", Decomposition)
  end
  # d_csp is the cospan dx₁ -> de <- dx₂ corresp to some edge e = x₁x₂ in shape(d)
  (csp, d_csp)      = tup  #unpack the tuple
  # the pullback cone dx₁ <-l₁-- p --l₂ --> dx₂ with legs l₁ and l₂
  p_cone            = pullback(d_csp)
  p_legs            = legs(p_cone)
  #for each leg lᵢ : p → xᵢ of the pullback cone, 
  #compute its image ιᵢ : im lᵢ → dxᵢ
  imgs              = map( f -> legs(image(f))[1], p_legs)
  #now get the new desired cospan; 
  #i.e.  im l₁ --ι₁--> dx₁ --l₁--> de <--l₂--dx₂ <--ι₂-- im l₂
  new_d_csp         = map(t -> compose(t...), zip(imgs, d_csp))  
  #get the domain of d 
  d_dom             = dom(d.diagram)
  #now make the new decomposition, call it δ
  #start with the object map δ₀
  function ob_replace(x)
    if x == dom(d_dom, csp[1])
      dom(new_d_csp[1])
    elseif x == dom(d_dom, csp[2])
      dom(new_d_csp[2])
    else 
      ob_map(d,x) 
    end
  end
  δ₀ = Dict( x => ob_replace(x) for x ∈ ob_generators(d_dom) )
  #now do the same thing with the morphism map
  function mor_replace(f) 
    if f == csp[1]
      return new_d_csp[1]
    elseif f == csp[2]
      return new_d_csp[2]
    else
      return hom_map(d,f)
    end 
  end
  δ₁ = Dict( f => mor_replace(f) for f ∈ hom_generators(d_dom) )
  StrDecomp(d.decomp_shape, FinDomFunctor(δ₀, δ₁, d.domain), d.decomp_type)
end

#for some reason PartialFunctions is giving me an error here 
#and we have to explicitly Curry adhesion_filter.. 
adhesion_filter(tup::Tuple) = d -> adhesion_filter(tup, d) 

"""Solve the decision problem encoded by a sheaf. 
The algorithm is as follows: 
  compute on each bag (optionally, if the decomposition of the solution space
                        is already known, then it can be passed as an argument),
  compute composites on edges, 
  project back down to bags
  answer (providing a witness)
    "no" if there is an empty bag; 
    "yes" otherwise.
"""
function decide_sheaf_tree_shape(f, d::StructuredDecomposition, solution_space_decomp::StructuredDecomposition = 𝐃(f, d, CoDecomposition))
  witness = foldl(∘, map(adhesion_filter, adhesionSpans(solution_space_decomp, true)))(solution_space_decomp)
  (foldr(&, map( !isempty, bags(witness))), witness)
end


end