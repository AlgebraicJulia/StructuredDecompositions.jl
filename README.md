# StructuredDecompositions.jl

Many graph algorithms can be sped up be working with [tree decompositions][1] of graphs. These are generalized by [structured decompositions][2], which can be formed for arbitary data structures. StructuredDecompositions.jl is a package for manupulating and working with structured decompositions. It is part of the AlgebraicJulia ecosystem.

# Sheaves

Since the definition of a structured decompositions is functorial, one can easily lift computational problems (defined as functors mapping inputs to solution spaces) to functors between from decompositions of the inputs to decompositions of solution spaces. This package allows one to leverage insights to solve decision problems that are encoded as [sheaves][3] efficiently (i.e. in [fixed-parameter-tractable][4] time parameterized by the width of the decompositions).

# Benchmarks

## Junction Tree Construction

| library | name | supernode partition | vertices | edges | time | gctime | memory | allocs |
| :------ | :--- | :------------------ | :------- | :---- | :--- | :----- | :----- | :----- |
| StructuredDecompositions | mycielskian2 | maximal | 2 | 1 | 1.508 μs | 0.000 ns | 3.69 KiB | 94 |
| StructuredDecompositions | mycielskian2 | fundamental | 2 | 1 | 1.508 μs | 0.000 ns | 3.69 KiB | 94 |
| StructuredDecompositions | mycielskian2 | nodal | 2 | 1 | 1.254 μs | 0.000 ns | 3.34 KiB | 83 |
| QDLDL | mycielskian2 |      | 2 | 1 | 604.539 ns | 0.000 ns | 1.81 KiB | 49 |
| StructuredDecompositions | mycielskian4 | maximal | 11 | 23 | 2.356 μs | 0.000 ns | 6.70 KiB | 97 |
| StructuredDecompositions | mycielskian4 | fundamental | 11 | 23 | 2.532 μs | 0.000 ns | 6.97 KiB | 97 |
| StructuredDecompositions | mycielskian4 | nodal | 11 | 23 | 2.282 μs | 0.000 ns | 6.33 KiB | 83 |
| QDLDL | mycielskian4 |      | 11 | 23 | 1.017 μs | 0.000 ns | 4.70 KiB | 49 |
| StructuredDecompositions | dwt_59 | maximal | 59 | 104 | 7.364 μs | 0.000 ns | 25.78 KiB | 97 |
| StructuredDecompositions | dwt_59 | fundamental | 59 | 104 | 7.531 μs | 0.000 ns | 25.78 KiB | 97 |
| StructuredDecompositions | dwt_59 | nodal | 59 | 104 | 7.239 μs | 0.000 ns | 23.89 KiB | 83 |
| QDLDL | dwt_59 |      | 59 | 104 | 3.479 μs | 0.000 ns | 19.86 KiB | 49 |
| StructuredDecompositions | can_292 | maximal | 292 | 1124 | 46.917 μs | 0.000 ns | 146.36 KiB | 132 |
| StructuredDecompositions | can_292 | fundamental | 292 | 1124 | 47.625 μs | 0.000 ns | 148.39 KiB | 132 |
| StructuredDecompositions | can_292 | nodal | 292 | 1124 | 50.250 μs | 0.000 ns | 146.75 KiB | 121 |
| QDLDL | can_292 |      | 292 | 1124 | 28.083 μs | 0.000 ns | 146.08 KiB | 69 |
| StructuredDecompositions | lshp3466 | maximal | 3466 | 10215 | 642.666 μs | 0.000 ns | 1.49 MiB | 142 |
| StructuredDecompositions | lshp3466 | fundamental | 3466 | 10215 | 637.292 μs | 0.000 ns | 1.49 MiB | 142 |
| StructuredDecompositions | lshp3466 | nodal | 3466 | 10215 | 895.542 μs | 0.000 ns | 2.38 MiB | 122 |
| QDLDL | lshp3466 |      | 3466 | 10215 | 787.583 μs | 0.000 ns | 2.32 MiB | 70 |
| StructuredDecompositions | wing | maximal | 62032 | 121544 | 18.386 ms | 0.000 ns | 28.59 MiB | 143 |
| StructuredDecompositions | wing | fundamental | 62032 | 121544 | 18.762 ms | 0.000 ns | 29.89 MiB | 143 |
| StructuredDecompositions | wing | nodal | 62032 | 121544 | 54.225 ms | 0.000 ns | 179.76 MiB | 123 |
| QDLDL | wing |      | 62032 | 121544 | 96.856 ms | 0.000 ns | 177.01 MiB | 71 |
| StructuredDecompositions | 144 | maximal | 144649 | 1074393 | 67.015 ms | 0.000 ns | 98.68 MiB | 143 |
| StructuredDecompositions | 144 | fundamental | 144649 | 1074393 | 67.013 ms | 0.000 ns | 98.85 MiB | 143 |
| StructuredDecompositions | 144 | nodal | 144649 | 1074393 | 474.307 ms | 1.025 ms | 1.45 GiB | 123 |
| QDLDL | 144 |      | 144649 | 1074393 | 1.123 s | 15.623 ms | 1.47 GiB | 71 |
| StructuredDecompositions | 333SP | maximal | 3712815 | 11108633 | 1.209 s | 33.611 ms | 1.64 GiB | 143 |
| StructuredDecompositions | 333SP | fundamental | 3712815 | 11108633 | 1.200 s | 14.837 ms | 1.64 GiB | 143 |
| StructuredDecompositions | 333SP | nodal | 3712815 | 11108633 | 2.047 s | 11.187 ms | 3.95 GiB | 123 |
| QDLDL | 333SP |      | 3712815 | 11108633 | 2.566 s | 5.187 ms | 3.89 GiB | 71 |

## Vertex Coloring

| library | colors | bags | time | gctime | memory | allocs |
| :------ | :----- | :--- | :----| :----- | :----- | :----- |
| StructuredDecompositions | 2 | 1 | 4.230 ms | 0.000 ns | 3.07 MiB | 18444 |
| StructuredDecompositions | 2 | 2 | 3.986 ms | 0.000 ns | 2.93 MiB | 16935 |
| StructuredDecompositions | 2 | 4 | 5.520 ms | 0.000 ns | 3.95 MiB | 25902 |
| StructuredDecompositions | 2 | 8 | 7.420 ms | 0.000 ns | 5.17 MiB | 39942 |
| StructuredDecompositions | 2 | 12 | 8.208 ms | 0.000 ns | 5.57 MiB | 50131 |
| Catlab | 2 |     | 1.808 ms | 0.000 ns | 1.37 MiB | 6528 |
| SimpleGraphAlgorithms | 2 |     | 7.716 ns | 0.000 ns | 0 bytes | 0 |
| StructuredDecompositions | 3 | 1 | 133.087 s | 6.036 s | 85.64 GiB | 472577846 |
| StructuredDecompositions | 3 | 2 | 1.502 s | 54.585 ms | 972.17 MiB | 12140720 |
| StructuredDecompositions | 3 | 4 | 242.934 ms | 0.000 ns | 168.18 MiB | 952391 |
| StructuredDecompositions | 3 | 8 | 119.023 ms | 0.000 ns | 83.41 MiB | 500101 |
| StructuredDecompositions | 3 | 12 | 66.521 ms | 0.000 ns | 46.25 MiB | 330108 |
| Catlab | 3 |     | 28.519 ms | 0.000 ns | 21.09 MiB | 99689 |
| SimpleGraphAlgorithms | 3 |     | 7.717 ns | 0.000 ns | 0 bytes | 0 |

  [1]: https://en.wikipedia.org/wiki/Tree_decomposition
  [2]: https://arxiv.org/abs/2207.06091
  [3]: https://en.wikipedia.org/wiki/Sheaf_(mathematics)
  [4]: https://en.wikipedia.org/wiki/Parameterized_complexity
