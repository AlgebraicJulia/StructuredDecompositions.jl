# StructuredDecompositions.jl

Many graph algorithms can be sped up be working with [tree decompositions][1] of graphs. These are generalized by [structured decompositions][2], which can be formed for arbitary data structures. StructuredDecompositions.jl is a package for manupulating and working with structured decompositions. It is part of the AlgebraicJulia ecosystem.

# Sheaves

Since the definition of a structured decompositions is functorial, one can easily lift computational problems (defined as functors mapping inputs to solution spaces) to functors between from decompositions of the inputs to decompositions of solution spaces. This package allows one to leverage insights to solve decision problems that are encoded as [sheaves][3] efficiently (i.e. in [fixed-parameter-tractable][4] time parameterized by the width of the decompositions).

# Benchmarks

## Junction Tree Construction

| library | name | vertices | edges | time | gctime | memory | allocs |
| :------ | :--- | :--------| :-----| :----| :----- | :----- | :----- |
| StructuredDecompositions | mycielskian2 | 2 | 1 | 1.542 μs | 0.000 ns | 3.69 KiB | 94 |
| QDLDL | mycielskian2 | 2 | 1 | 590.199 ns | 0.000 ns | 1.81 KiB | 49 |
| StructuredDecompositions | mycielskian4 | 11 | 23 | 2.398 μs | 0.000 ns | 6.70 KiB | 97 |
| QDLDL | mycielskian4 | 11 | 23 | 995.900 ns | 0.000 ns | 4.70 KiB | 49 |
| StructuredDecompositions | dwt_59 | 59 | 104 | 7.646 μs | 0.000 ns | 25.78 KiB | 97 |
| QDLDL | dwt_59 | 59 | 104 | 3.453 μs | 0.000 ns | 19.86 KiB | 49 |
| StructuredDecompositions | can_292 | 292 | 1124 | 46.708 μs | 0.000 ns | 146.36 KiB | 132 |
| QDLDL | can_292 | 292 | 1124 | 27.750 μs | 0.000 ns | 146.08 KiB | 69 |
| StructuredDecompositions | lshp3466 | 3466 | 10215 | 637.667 μs | 0.000 ns | 1.49 MiB | 142 |
| QDLDL | lshp3466 | 3466 | 10215 | 789.083 μs | 0.000 ns | 2.32 MiB | 70 |
| StructuredDecompositions | wing | 62032 | 121544 | 18.225 ms | 0.000 ns | 28.59 MiB | 143 |
| QDLDL | wing | 62032 | 121544 | 97.209 ms | 0.000 ns | 177.01 MiB | 71 |
| StructuredDecompositions | 144 | 144649 | 1074393 | 67.180 ms | 0.000 ns | 98.68 MiB | 143 |
| QDLDL | 144 | 144649 | 1074393 | 1.082 s | 757.083 μs | 1.47 GiB | 71 |
| StructuredDecompositions | 333SP | 3712815 | 11108633 | 1.191 s | 13.445 ms | 1.64 GiB | 143 |
| QDLDL | 333SP | 3712815 | 11108633 | 2.326 s | 5.432 ms | 3.89 GiB | 71 |

## Vertex Coloring

| library | colors | bags | time | gctime | memory | allocs |
| :------ | :----- | :--- | :----| :----- | :----- | :----- |
| StructuredDecompositions | 2 | 1 | 4.220 ms | 0.000 ns | 3.07 MiB | 18444 |
| StructuredDecompositions | 2 | 2 | 3.997 ms | 0.000 ns | 2.93 MiB | 16935 |
| StructuredDecompositions | 2 | 4 | 5.500 ms | 0.000 ns | 3.95 MiB | 25902 |
| StructuredDecompositions | 2 | 8 | 7.358 ms | 0.000 ns | 5.17 MiB | 39942 |
| StructuredDecompositions | 2 | 12 | 8.137 ms | 0.000 ns | 5.57 MiB | 50131 |
| Catlab | 2 |     | 1.822 ms | 0.000 ns | 1.37 MiB | 6528 |
| SimpleGraphAlgorithms | 2 |     | 7.708 ns | 0.000 ns | 0 bytes | 0 |
| StructuredDecompositions | 3 | 1 | 129.478 s | 5.500 s | 85.64 GiB | 472577846 |
| StructuredDecompositions | 3 | 2 | 1.482 s | 50.695 ms | 972.17 MiB | 12140720 |
| StructuredDecompositions | 3 | 4 | 241.190 ms | 0.000 ns | 168.18 MiB | 952391 |
| StructuredDecompositions | 3 | 8 | 117.753 ms | 0.000 ns | 83.41 MiB | 500101 |
| StructuredDecompositions | 3 | 12 | 65.992 ms | 0.000 ns | 46.25 MiB | 330108 |
| Catlab | 3 |     | 28.544 ms | 0.000 ns | 21.09 MiB | 99689 |
| SimpleGraphAlgorithms | 3 |     | 7.708 ns | 0.000 ns | 0 bytes | 0 |

  [1]: https://en.wikipedia.org/wiki/Tree_decomposition
  [2]: https://arxiv.org/abs/2207.06091
  [3]: https://en.wikipedia.org/wiki/Sheaf_(mathematics)
  [4]: https://en.wikipedia.org/wiki/Parameterized_complexity
