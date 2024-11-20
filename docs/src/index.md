# StructuredDecompositions.jl

```@meta
CurrentModule = StructuredDecompositions
```

Structural graph theorists and algorithmicists alike know that it's usually a smart idea to decompose graphs into smaller and simpler parts before trying answer difficult computational questions. [Tree decompositions](https://en.wikipedia.org/wiki/Tree_decomposition) are one of the best-known ways of systematically chopping graphs up and, as such they have been key tools for establishing deep results in many areas of discrete mathematics including graph minor theory and algorithmic meta-theorems. 

## Generality
  
This project consists of an implementation of the category-theoretic notion of [structured decompositions](https://arxiv.org/abs/2207.06091). These provide a formalism for decomposing arbitrary mathematical objects (not just graphs) and morphisms between them into smaller constituent parts. Since the definition of a structured decompositions is functorial, one can easily lift computational problems (defined as functors mapping inputs to solution spaces) to functors between from decompositions of the inputs to decompositions of solution spaces.

## What is in this package?

This package allows one to leverage insights to solve decision problems that are encoded as [sheaves](https://en.wikipedia.org/wiki/Sheaf_(mathematics)) efficiently (i.e. in [fixed-parameter-tractable](https://en.wikipedia.org/wiki/Parameterized_complexity) time parameterized by the width of the decompositions). Currently, this packages includes many general tools that can be used to decompose arbitrary mathematical objects. One of the many applications of this package, exampled in the Deciding Sheaves module, is graph colorings.

```@contents
Pages = [
    "pages/decompostions.md",
    "pages/decidingsheaves.md",
    "pages/junction_trees.md",
    "pages/nested_uwds.md",
    "pages/functor_utils.md",
    ]
Depth = 2
```
