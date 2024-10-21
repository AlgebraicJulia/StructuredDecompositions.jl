module DecidingSheaves

export Presheaf, Sheaf, decide_sheaf_tree_shape

using ..Decompositions
using ..FunctorUtils

using Catlab
using Catlab.CategoricalAlgebra

struct DecompError <: Exception end

Base.showerror(io, e::DecompError) = print(io, "Expecting Codecomposition but received decomposition")

# TODO FinSet(Int) assumed
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
    d.decomp_type == Decomposition && throw(DecompError)
    # d_cospan is the cospan dx₁ -> de <- dx₂ corresp to some edge e = x₁x₂ in shape(d)
    (cospan, d_cospan) = tup  #unpack the tuple
    # the pullback cone dx₁ <-l₁-- p --l₂ --> dx₂ with legs l₁ and l₂
    p_legs             = (legs∘pullback)(d_cospan)
    # for each leg lᵢ : p → xᵢ of the pullback cone, 
    # compute its image ιᵢ : im lᵢ → dxᵢ
    imgs               = (first∘legs∘image).(p_legs)
    # now get the new desired cospan; 
    # i.e.  im l₁ --ι₁--> dx₁ --l₁--> de <--l₂-- dx₂ <--ι₂-- im l₂
    new_d_cospan       = map(t -> compose(t...), zip(imgs, d_cospan))  
    # get the domain of d 
    d_dom              = dom(d.diagram)
  
    # TODO is there ever a time when "out" has length > 1 
    function ob_replace(x)
        out = dom.(new_d_cospan[x .== dom.(Ref(d_dom), cospan)])
        !isempty(out) ? first(out) : ob_map(d, x)
    end

    function mor_replace(f)
        out = new_d_cospan[f .== cospan]
        !isempty(out) ? first(out) : hom_map(d, f)
    end

    # now make the new decomposition, call it δ
    # start with the object map δ₀
    δ₀ = Dict( x => ob_replace(x) for x ∈ ob_generators(d_dom) )
    # now do the same thing with the morphism map
    δ₁ = Dict( f => mor_replace(f) for f ∈ hom_generators(d_dom) )

    StrDecomp(d.decomp_shape, FinDomFunctor(δ₀, δ₁, d.domain), d.decomp_type)
end

# for some reason PartialFunctions is giving me an error here 
# and we have to explicitly Curry adhesion_filter.. 
adhesion_filter(tup::Tuple) = d -> adhesion_filter(tup, d) 

export adhesion_filter

"""    decide_sheaf_tree_shape(f, 
                d::StructuredDecomposition, 
                solution_space_decomp::StructuredDecomposition)

Solve the decision problem encoded by a sheaf. 
The algorithm is as follows: 
  1. compute on each bag. Optionally, if the decomposition of the solution space
     is already known, then it can be passed as an argument.
  2. compute composites on edges 
  3. project back down to bags
  4. answer (providing a witness)
    "no" if there is an empty bag; 
    "yes" otherwise.

"""
function decide_sheaf_tree_shape(f,
        d::StructuredDecomposition,
        solution_space_decomp::StructuredDecomposition = 𝐃(f, d, CoDecomposition))

    # witness = foldl(∘, map(adhesion_filter, adhesionSpans(solution_space_decomp, true)))(solution_space_decomp)

    @info "Structured Space Decomposition: $solution_space_decomp"
    @info "Adhesion Spans: $(adhesionSpans(solution_space_decomp, true))"
    witness = 
    ∘((adhesion_filter.(adhesionSpans(solution_space_decomp, true)))...)(solution_space_decomp)
  
    (all(!isempty, bags(witness)), witness)
end


end
