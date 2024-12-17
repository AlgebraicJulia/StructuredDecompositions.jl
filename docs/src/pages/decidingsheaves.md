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

We can consider the example of a graph with 5 vertices forming an 'X'. 
We first seperate the nodes into bags with adhesions and what our adhesions look like.
Notice that the number of nodes in the graph is equal to the sum of the nodes in each bag minus the adhesions between each bag.

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
    V = 1
end

#bag 2
H₂ = @acset Graph begin
    V = 3
    E = 2
    src = [1, 2]
    tgt = [2, 3]
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
      1 => ACSetTransformation(Γₛ⁰[3], Γₛ⁰[1], V=[2]),
      2 => ACSetTransformation(Γₛ⁰[3], Γₛ⁰[2], V=[2]),
    ),
    ∫(Gₛ)
)

my_decomp1  = StrDecomp(Gₛ, Γₛ)
```

The complete graph is pictured below with the adhesion vertex colored yellow.

<img src='docs\src\assets\example graph.png' width='512' alt='Complete example graph'>

We can recognize the bags that are formed. The bags are pictured as red and blue circles below, overlapping at the yellow adhesion vertex.

<img src='docs\src\assets\example colored bags.png' width='512' alt='Example graph with highlighted bags'>

We can thus split the graph by these bags to form the decomposition shown below.

<img src='docs\src\assets\example bags and adhesion.png' width='512' alt='Example graph split by bags and adhesions'>

and the resulting tree decomposition below.

<img src='docs\src\assets\example tree.png' width='512' alt='Example tree decomposition'>

The algorithm will work by coloring the individual bags then checking if the results will connect at the adhesion point. If no match is found at the adhesion point then there will be a no result. 

## Understanding the Complexity

For a c-coloring with tree-width w and total vertices |VT|, we know that to determine a brute force solution we are at complexity O(c<sup>|VT|</sup>). With this dynamic solution we expect complexity of O(c<sup>2w</sup>)|VT|.

We will also consider some other factors. The brute force solution will be better in some edge cases. In a yes instance, if brute force is looking for only 1 solution, a yes instance could get lucky, especially if there are many possible c-colorings. In a yes instance, if brute force is looking for all solutions, the dynamic solution will always be slower in the 2-bag case which can be easily confirmed in the expected complexity calculation. In a no instance, again, the dynamic solution will be slower in the 2-bag case due to the same reason.

## Benchmarking Complexity

As stated above, due to the structure of the algorithm, we expect a complexity of O(c<sup>2w</sup>)|VT| for the dynamic solution. Thus we want to see an exponential decrease in run time as a fixed graph is decomposed into more bags. We also want to show a linear increase in runtime as the number of bags increases linearly and we want to show an exponential increase in runtime as the number of vertices in the graph increases linearly. 

First consider the fixed graph case with the following fixed graph.

<img src='docs\src\assets\fixed graph.png' width='512' alt='Example fixed graph'>

We considered the runtime when the graph above is decomposed into smaller and smaller bags as shown below.

<img src='docs\src\assets\fixed graph bags.png' width='512' alt='Example fixed graph decompositions'>

The runtimes for computing if the fixed graph is 3-colorable with varying bag sizes is shown below.

<img src='docs\src\assets\Fixed Graph 3-Coloring.png' width='256' alt='Runtime graph for fixed graph 3-coloring'>

The results show an exponential decrease in runtime as the number of bags is increased on a fixed graph as expected from the complexity calculation above.

Now consider the following benchmarks from increasing the number of bags of a fixed size at a linear rate. When using the following bag of size 4 with some examples below

<img src='docs\src\assets\increasing bag.png' width='128' alt='An example of 1 possible bag used in the benchmarking'>

<img src='docs\src\assets\bags increasing example.png' width='128' alt='Example graphs used in the increasing bag case'>

The following are the results for fixed bags of size 4 and 8 for 2, 3, and 4 coloring tests.

First the results for the fixed bag of size 4 and for a 2-coloring test. The graph shows a linear increase in runtime.

<img src='docs\src\assets\Bags Increasing Linearly of size 4 2 color.png' width='256' alt='Bags increasing size 4 2 coloring results'>

The following results are for the fixed bag of size 4 and for a 3-coloring test. The graph shows a linear increase in runtime.

<img src='docs\src\assets\Bags Increasing Linearly of size 4 3 color.png' width='256' alt='Bags increasing size 4 3 coloring results'>

The following results are for the fixed bag of size 4 and for a 4-coloring test. The graph shows a linear increase in runtime.

<img src='docs\src\assets\Bags Increasing Linearly of size 4 4 color.png' width='256' alt='Bags increasing size 4 4 coloring results'>

The following results are for the fixed bag of size 8 and for a 2-coloring test. The graph shows a linear increase in runtime.

<img src='docs\src\assets\Bags Increasing Linearly of size 8 2 color.png' width='256' alt='Bags increasing size 8 2 coloring results'>

The following results are for the fixed bag of size 8 and for a 3-coloring test. The graph shows a linear increase in runtime.

<img src='docs\src\assets\Bags Increasing Linearly of size 8 3 color.png' width='256' alt='Bags increasing size 8 3 coloring results'>

The following results are for the fixed bag of size 8 and for a 4-coloring test. The graph shows a linear increase in runtime.

<img src='docs\src\assets\Bags Increasing Linearly of size 8 4 color.png' width='256' alt='Bags increasing size 8 4 coloring results'>

These results all confirm the expected complexity from the calculation above.

Now finally consider the following benchmarks for increasing the total vertices in a graph while keeping the number of bags at 2 and the number of adhesions at 1. Below are some examples to illustrate the idea.

<img src='docs\src\assets\nodes increasing.png' width='256' alt='Example of graphs with increasing nodes'>

The following are the results for 2, 3, and 4 coloring tests with similar graphs.

<img src='docs\src\assets\Increase Vertices Linearly 2 color.png' width='256' alt='Nodes increasing results for 2 coloring'>

<img src='docs\src\assets\Increase Vertices Linearly 3 color.png' width='256' alt='Nodes increasing results for 3 coloring'>

<img src='docs\src\assets\Increase Vertices Linearly 4 color.png' width='256' alt='Nodes increasing results for 4 coloring'>