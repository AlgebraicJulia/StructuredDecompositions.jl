# [Deciding Sheaves] (@id DecidingSheaves)

There are two functions that are used here. 

The first one called adhesion_filter will take an input Finset^{op}-valued structured decomposition and return a new structured decompostion replacing the span de in d by the span obtained by projecting the pullback of de. 

```julia
function adhesion_filter(tup::Tuple, d::StructuredDecomposition)
```

The second function is called decide_sheaf_tree_shape and solves a decision problem that is encoded by a sheaf. This function works by computing first on each of the bags, then computes composites on edges and projects back down to bags.

```julia
function decide_sheaf_tree_shape(f, d::StructuredDecomposition, solution_space_decomp::StructuredDecomposition = 𝐃(f, d, CoDecomposition))
```

## Graph Coloring Structure and Examples

One of the many decision problems that can be solved using this system of functions is figuring out graph colorings. To put it simply, a graph coloring is an assignment of labels(colors) to each vertex of a graph so that no two adjacent vertices have the same label(color).

We first define the structure of a coloring.

```julia
#an H-coloring is a hom onto H
struct Coloring
  n     #the target graph
  func  #the function mapping opens to lists of homs from G to K_n
end
```

We then define how we are going to create and test colorings.

```julia
#construct an n-coloring
K(n) = complete_graph(Graph, n)
Coloring(n) = Coloring(n, g -> homomorphisms(g, K(n) ))
#make it callable
(c::Coloring)(X::Graph) = FinSet(c.func(X)) # given graph homos #f: G₁ → G₂ get morphism col(G₂) → col(G₁) by precomposition: take each λ₂ ∈ col(G₂) to hf ∈ col(G)
function (c::Coloring)(f::ACSetTransformation)  
  (G₁, G₂)   = (dom(f), codom(f)) 
  (cG₁, cG₂) = (c(G₁), c(G₂))
  FinFunction( λ₂ -> compose(f,λ₂), cG₂, cG₁ ) #note the contravariance
end

skeletalColoring(n) = skeleton ∘ Coloring(n)

colorability_test(n, the_test_case) = is_homomorphic(ob(colimit(the_test_case)), K(n)) == decide_sheaf_tree_shape(skeletalColoring(n), the_test_case)[1]
```

We can consider the example of a ring with seven nodes as our graph. 
We first seperate the nodes into bags with adhesions and what our adhesions look like.

```julia
#bag 1
H₁ = @acset Graph begin
    V = 3
    E = 2
    src = [1, 2]
    tgt = [2, 3]
end

#adhesion 1,2
H₁₂ = @acset Graph begin
    V = 2
end

#bag 2
H₂ = @acset Graph begin
    V = 4
    E = 3
    src = [1, 2, 3]
    tgt = [2, 3, 4]
end

Gₛ = @acset Graph begin
    V = 2
    E = 1
    src = [1]
    tgt = [2]
end
```

We can then construct our structured decomposition through transformations of these bags and adhesions.

```julia
#transformations
Γₛ⁰ = Dict(1 => H₁, 2 => H₂, 3 => H₁₂)
Γₛ = FinDomFunctor(
    Γₛ⁰,
    Dict(
      1 => ACSetTransformation(Γₛ⁰[3], Γₛ⁰[1], V=[1, 3]),
      2 => ACSetTransformation(Γₛ⁰[3], Γₛ⁰[2], V=[4, 1]),
    ),
    ∫(Gₛ)
)

my_decomp1  = StrDecomp(Gₛ, Γₛ)
```