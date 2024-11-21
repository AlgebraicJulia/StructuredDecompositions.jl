# StructuredDecompositions.jl

HELLO

Structural graph theorists and algorithmicists alike know that it's usually a smart idea to decompose graphs into smaller and simpler parts before trying answer difficult computational questions. [Tree decompositions][1] are one of the best-known ways of systematically chopping graphs up and, as such they have been key tools for establishing deep results in many areas of discrete mathematics including graph minor theory and algorithmic meta-theorems.

### But what happens if we want to compute on other kinds of mathematical structures?

This project consists of an implementation of the category-theoretic notion of [structured decompositions][2]. These provide a formalism for decomposing decomposing arbitrary mathematical objects (not just graphs) and morphisms between them into smaller constituent parts. Since the definition of a structured decompositions is functorial, one can easily lift computational problems (defined as functors mapping inputs to solution spaces) to functors between from decompositions of the inputs to decompositions of solution spaces. This package allows one to leverage insights to solve decision problems that are encoded as [sheaves][3] efficiently (i.e. in [fixed-parameter-tractable][4] time parameterized by the width of the decompositions).


[1]: https://en.wikipedia.org/wiki/Tree_decomposition
[2]: https://arxiv.org/abs/2207.06091
[3]: https://en.wikipedia.org/wiki/Sheaf_(mathematics)
[4]: https://en.wikipedia.org/wiki/Parameterized_complexity
