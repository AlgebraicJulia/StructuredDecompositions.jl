# StructuredDecompositions.jl

Many graph algorithms can be sped up be working with [tree decompositions][1] of graphs. These are generalized by [structured decompositions][2], which can be formed for arbitary data structures. StructuredDecompositions.jl is a package for manupulating and working with structured decompositions. It is part of the AlgebraicJulia ecosystem.

# Sheaves

Since the definition of a structured decompositions is functorial, one can easily lift computational problems (defined as functors mapping inputs to solution spaces) to functors between from decompositions of the inputs to decompositions of solution spaces. This package allows one to leverage insights to solve decision problems that are encoded as [sheaves][3] efficiently (i.e. in [fixed-parameter-tractable][4] time parameterized by the width of the decompositions).

# Benchmarks

You can make a file `benchmarks.md` with the following tables by navigating to the project directory and typing the following in to your command line.
```
julia makebenchmarks.jl benchmarks.md
```

## 2 Coloring

| library | bags | time | gctime | memory | allocs |
| :------ | :--- | :----| :----- | :----- | :----- |
| StructuredDecompositions | 1 | 3.208 ms | 0.000 ns | 3.07 MiB | 18444 |
| StructuredDecompositions | 2 | 3.006 ms | 0.000 ns | 2.93 MiB | 16935 |
| StructuredDecompositions | 4 | 4.211 ms | 0.000 ns | 3.95 MiB | 25902 |
| StructuredDecompositions | 8 | 5.770 ms | 0.000 ns | 5.17 MiB | 39942 |
| StructuredDecompositions | 12 | 6.530 ms | 0.000 ns | 5.57 MiB | 50131 |
| HomSearch |  | 1.348 ms | 0.000 ns | 1.37 MiB | 6528 |

## 3 Coloring

| library | bags | time | gctime | memory | allocs |
| :------ | :--- | :----| :----- | :----- | :----- |
| StructuredDecompositions | 1 | 98.464 s | 5.066 s | 85.64 GiB | 472577846 |
| StructuredDecompositions | 2 | 1.141 s | 42.218 ms | 972.17 MiB | 12140720 |
| StructuredDecompositions | 4 | 182.616 ms | 5.627 ms | 168.18 MiB | 952391 |
| StructuredDecompositions | 8 | 91.642 ms | 2.794 ms | 83.41 MiB | 500101 |
| StructuredDecompositions | 12 | 52.070 ms | 0.000 ns | 46.25 MiB | 330108 |
| HomSearch |  | 21.064 ms | 0.000 ns | 21.09 MiB | 99689 |

## 4 Coloring

| library | bags | time | gctime | memory | allocs |
| :------ | :--- | :--- | :----- | :----- | :----- |
| StructuredDecompositions | 4 | 156.684 s | 8.318 s | 119.82 GiB | 6085756624 |
| StructuredDecompositions | 8 | 2.009 s | 83.526 ms | 1.61 GiB | 14517740 |
| StructuredDecompositions | 12 | 354.173 ms | 14.021 ms | 292.50 MiB | 2184943 |
| HomSearch |  | 41.627 ms | 877.875 μs | 40.61 MiB | 191057 |


  [1]: https://en.wikipedia.org/wiki/Tree_decomposition
  [2]: https://arxiv.org/abs/2207.06091
  [3]: https://en.wikipedia.org/wiki/Sheaf_(mathematics)
  [4]: https://en.wikipedia.org/wiki/Parameterized_complexity
