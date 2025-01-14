# StructuredDecompositions.jl

Many graph algorithms can be sped up be working with [tree decompositions][1] of graphs. These are generalized by [structured decompositions][2], which can be formed for arbitary data structures. StructuredDecompositions.jl is a package for manupulating and working with structured decompositions. It is part of the AlgebraicJulia ecosystem.

# Sheaves

Since the definition of a structured decompositions is functorial, one can easily lift computational problems (defined as functors mapping inputs to solution spaces) to functors between from decompositions of the inputs to decompositions of solution spaces. This package allows one to leverage insights to solve decision problems that are encoded as [sheaves][3] efficiently (i.e. in [fixed-parameter-tractable][4] time parameterized by the width of the decompositions).

# Benchmarks

You can make a file `benchmarks.md` with the following tables by navigating to the `benchmark` directory and running the following command.
```
julia --project make.jl
```

## 2 Coloring

| library | bags | time | gctime | memory | allocs |
| :------ | :--- | :----| :----- | :----- | :----- |
| StructuredDecompositions | 1 | 4.060 ms | 0.000 ns | 3.07 MiB | 18444 |
| StructuredDecompositions | 2 | 3.832 ms | 0.000 ns | 2.93 MiB | 16935 |
| StructuredDecompositions | 4 | 5.286 ms | 0.000 ns | 3.95 MiB | 25902 |
| StructuredDecompositions | 8 | 7.053 ms | 0.000 ns | 5.17 MiB | 39942 |
| StructuredDecompositions | 12 | 7.865 ms | 0.000 ns | 5.57 MiB | 50131 |
| HomSearch |  | 1.732 ms | 0.000 ns | 1.37 MiB | 6528 |
| SimpleGraphAlgorithms |  | 7.716 ns | 0.000 ns | 0 bytes | 0 |

## 3 Coloring

| library | bags | time | gctime | memory | allocs |
| :------ | :--- | :----| :----- | :----- | :----- |
| StructuredDecompositions | 1 | 123.718 s | 5.662 s | 85.64 GiB | 472577846 |
| StructuredDecompositions | 2 | 1.428 s | 45.898 ms | 972.17 MiB | 12140720 |
| StructuredDecompositions | 4 | 230.669 ms | 6.021 ms | 168.18 MiB | 952391 |
| StructuredDecompositions | 8 | 114.980 ms | 2.573 ms | 83.41 MiB | 500101 |
| StructuredDecompositions | 12 | 64.231 ms | 0.000 ns | 46.25 MiB | 330108 |
| HomSearch |  | 27.330 ms | 0.000 ns | 21.09 MiB | 99689 |
| SimpleGraphAlgorithms |  | 7.674 ns | 0.000 ns | 0 bytes | 0 |

## 4 Coloring

| library | bags | time | gctime | memory | allocs |
| :------ | :--- | :--- | :----- | :----- | :----- |
| StructuredDecompositions | 4 | 199.385 s | 9.183 s | 119.82 GiB | 6085756624 |
| StructuredDecompositions | 8 | 2.441 s | 95.750 ms | 1.61 GiB | 14517740 |
| StructuredDecompositions | 12 | 429.789 ms | 14.376 ms | 292.50 MiB | 2184943 |
| HomSearch |  | 52.656 ms | 0.000 ns | 40.61 MiB | 191057 |
| SimpleGraphAlgorithms |  | 7.675 ns | 0.000 ns | 0 bytes | 0 |


  [1]: https://en.wikipedia.org/wiki/Tree_decomposition
  [2]: https://arxiv.org/abs/2207.06091
  [3]: https://en.wikipedia.org/wiki/Sheaf_(mathematics)
  [4]: https://en.wikipedia.org/wiki/Parameterized_complexity
