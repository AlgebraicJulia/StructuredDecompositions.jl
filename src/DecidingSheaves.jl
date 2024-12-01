module DecidingSheaves

export Presheaf, Sheaf, decide_sheaf_tree_shape

using ..Decompositions

using Catlab
using Catlab.CategoricalAlgebra

# TODO type this
# for each leg lᵢ : p → xᵢ of the pullback cone, 
#   # compute its image ιᵢ : im lᵢ → dxᵢ
function images(d_cospan::Any)
    pb = pullback(d_cospan);
    if isempty(pb)
        d_dom = FinSet(0)
    else
        imgs = map(f -> first(legs(image(f))), legs(pb))
        new_dual_cospan = map(t -> compose(t...), zip(imgs, d_cospan))
    end
end

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
  (csp, d_cospan) = tup
  
  # the pullback cone dx₁ <-l₁-- p --l₂ --> dx₂ with legs l₁ and l₂
  new_d_cospan = images(d_cospan)

  # get the domain of d
  d_dom = dom(d.diagram)

  # now make the new decomposition, call it δ
  # us ob_replace and mor_replace
  function ob_replace(x)
    if x == dom(d_dom, csp[1])
      dom(new_d_cospan[1])
    elseif x == dom(d_dom, csp[2])
      dom(new_d_cospan[2])
    else
      ob_map(d,x)
    end
  end

  function mor_replace(f)
    if f == csp[1]
      return new_d_cospan[1]
    elseif f == csp[2]
      return new_d_cospan[2]
    else
      hom_map(d,f)
    end
  end

  # start with the object map δ₀
  δ₀ = Dict( x => ob_replace(x) for x ∈ ob_generators(d_dom))
  # now do the same thing with the morphism map
  δ₁ = Dict( f => mor_replace(f) for f ∈ hom_generators(d_dom))

  StrDecomp(d.decomp_shape, FinDomFunctor(δ₀, δ₁, d.domain), d.decomp_type)
end

"""Solve the decision problem encoded by a sheaf. 
The algorithm is as follows: 
  compute on each bag 
  (optionally, if the decomposition of the solution space
  is already known, then it can be passed as an argument),
  compute composites on edges, 
  project back down to bags
  answer (providing a witness)
    "no" if there is an empty bag; 
    "yes" otherwise.
"""
function decide_sheaf_tree_shape(f, d::StructuredDecomposition, solution_space_decomp::StructuredDecomposition = 𝐃(f, d, CoDecomposition))
  witness = solution_space_decomp
  adhesion_spans = adhesionSpans(solution_space_decomp, true)
  for adhesion in adhesion_spans
    witness = adhesion_filter(adhesion, witness)
    if any(isempty, bags(witness))
      return (false, witness)
    end
  end
  return (true, witness)
end

end
