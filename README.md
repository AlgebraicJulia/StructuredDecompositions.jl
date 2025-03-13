# StructuredDecompositions.jl

Many graph algorithms can be sped up be working with [tree decompositions][1] of graphs. These are generalized by [structured decompositions][2], which can be formed for arbitary data structures. StructuredDecompositions.jl is a package for manupulating and working with structured decompositions. It is part of the AlgebraicJulia ecosystem.

# Sheaves

Since the definition of a structured decompositions is functorial, one can easily lift computational problems (defined as functors mapping inputs to solution spaces) to functors between from decompositions of the inputs to decompositions of solution spaces. This package allows one to leverage insights to solve decision problems that are encoded as [sheaves][3] efficiently (i.e. in [fixed-parameter-tractable][4] time parameterized by the width of the decompositions).

  [1]: https://en.wikipedia.org/wiki/Tree_decomposition
  [2]: https://arxiv.org/abs/2207.06091
  [3]: https://en.wikipedia.org/wiki/Sheaf_(mathematics)
  [4]: https://en.wikipedia.org/wiki/Parameterized_complexity
