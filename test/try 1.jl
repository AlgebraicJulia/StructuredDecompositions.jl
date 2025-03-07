using StructuredDecompositions.Decompositions
using Catlab
using MLStyle
using PartialFunctions
using StructuredDecompositions.Decompositions
using StructuredDecompositions.DecidingSheaves
using StructuredDecompositions.FunctorUtils
using StructuredDecompositions.JunctionTrees
using Test


import TreeWidthSolver

function Catlab.WiringDiagramAlgebras.make_homomorphism(row::AbstractVector{T}, X::StructACSet{S}, Y::StructACSet{S}) where {T, S}
  components = let i = 0
    NamedTuple{ob(S)}(T[row[i+=1] for _ in parts(X,c)] for c in ob(S))
  end
  ACSetTransformation(components, X, Y)
end


function K(n::Integer)
    complete_graph(Graph, n)
end


struct Coloring
    n::Int
end


function (coloring::Coloring)(graph::Graph)
    FinSet(homomorphisms(graph, K(coloring.n); alg=HomomorphismQuery()))
end


function (coloring::Coloring)(f::ACSetTransformation)
    FinFunction(λ -> compose(f, λ), coloring(codom(f)), coloring(dom(f)))
end


function skeletal_coloring(n::Integer)
    skeleton ∘ Coloring(n)
end


function test_colorability(n::Integer, decomp::StrDecomp)
    left = is_homomorphic(ob(colimit(decomp)), K(n))
    right = first(decide_sheaf_tree_shape(skeletal_coloring(n), decomp))
    isequal(left, right)
end

# bag 1
    H1 = @acset Graph begin
        V = 3
        E = 2
        src = [1, 1]
        tgt = [2, 3]
    end

    # adhesion 1, 2
    H12 = @acset Graph begin
        V = 2
    end

    # bag 2
    H2 = @acset Graph begin
        V = 3
        E = 1
        src = [1]
        tgt = [3]
    end

    H23 = @acset Graph begin
        V = 2
    end

    H3 = @acset Graph begin
        V = 3
        E = 2
        src = [3, 3]
        tgt = [1, 2]
    end

    G = @acset Graph begin
        V = 3
        E = 2
        src = [1, 2]
        tgt = [2, 3]
    end

    # transformations
    Γ = FinDomFunctor(
        Dict(1 => H1, 2 => H2, 3 => H3, 4 => H12, 5 => H23),
        Dict(
            1 => ACSetTransformation(H12, H1, V=[2, 3]),
            2 => ACSetTransformation(H12, H2, V=[2, 1]),
            3 => ACSetTransformation(H23, H2, V=[2, 3]),
            4 => ACSetTransformation(H23, H3, V=[2, 1])),
        ∫(G))

    manual = StrDecomp(G, Γ)
    automatic = StrDecomp(ob(colimit(manual)))

decide_sheaf_tree_shape(skeletal_coloring(2), manual)[1] == false

decide_sheaf_tree_shape(skeletal_coloring(2), automatic)[1] == false

decide_sheaf_tree_shape(skeletal_coloring(3), manual)[1] == true

decide_sheaf_tree_shape(skeletal_coloring(3), automatic)[1] == true

all(test_colorability(n, manual) for n ∈ range(1, 3))

function adhesion_filter(tup::Tuple, d::StructuredDecomposition)
  if d.decomp_type == Decomposition
    error("expecting ", CoDecomposition, " given ", Decomposition)
  end
  # d_csp is the cospan dx₁ -> de <- dx₂ corresp to some edge e = x₁x₂ in shape(d)
  (csp, d_csp)      = tup  #unpack the tuple
  # the pullback cone dx₁ <-l₁-- p --l₂ --> dx₂ with legs l₁ and l₂
  p_cone            = pullback(d_csp)
  p_legs            = legs(p_cone)
  # for each leg lᵢ : p → xᵢ of the pullback cone, 
  # compute its image ιᵢ : im lᵢ → dxᵢ
  imgs              = force.(map( f -> legs(image(f))[1], p_legs))
  # now get the new desired cospan; 
  # i.e.  im l₁ --ι₁--> dx₁ --l₁--> de <--l₂--dx₂ <--ι₂-- im l₂
  new_d_csp         = force.(map(t -> compose(t...), zip(imgs, d_csp)))  
  # get the domain of d 
  d_dom             = dom(d.diagram)
  d_codom           = codom(d.diagram)
  # now make the new decomposition, call it δ
  # start with the object map δ₀
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
  # now do the same thing with the morphism map
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

  witness = 𝐃(skeletal_coloring(2), manual, CoDecomposition)

  adhesion_spans = adhesionSpans(𝐃(skeletal_coloring(2), manual, CoDecomposition), true)

  for adhesion in adhesion_spans
    witness = adhesion_filter(adhesion, witness)
    if any(isempty, bags(witness))
      return (false, witness)
    end
  end
  return (true, witness)

tup = ([1, 3], [FinFunction([4, 1], 2, 4), FinFunction([3, 4, 1, 2], 4, 4)])

d = witness

(csp, d_csp) = tup

p_cone = pullback(d_csp)

p_legs = legs(p_cone)

imgs = force.(map(f->legs(image(f))[1], p_legs))

new_d_csp = force.(map(t->compose(t...), zip(imgs, d_csp)))

d_dom = dom(d.diagram)

  function ob_replace(x)
    if x == dom(d_dom, csp[1])
      return dom(new_d_csp[1])
    elseif x == dom(d_dom, csp[2])
      return dom(new_d_csp[2])
    elseif x == codom(d_dom, csp[1])
      return codom(new_d_csp[1])
    else
      ob_map(d,x)
    end
  end

  function mor_replace(f)
    if f == csp[1]
      return new_d_csp[1]
    elseif f == csp[2]
      return new_d_csp[2]
    else
      return hom_map(d,f)
    end
  end

  # start with the object map δ₀
  δ₀ = Dict( x => ob_replace(x) for x ∈ ob_generators(d_dom))

  # now do the same thing with the morphism map
  δ₁ = Dict( f => mor_replace(f) for f ∈ hom_generators(d_dom))

StrDecomp(d.decomp_shape, FinDomFunctor(δ₀, δ₁, d.domain), d.decomp_type)

any(isempty, bags(witness))


