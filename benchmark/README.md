# Benchmarks

This file was automatically generated on 2025-01-28. To regenerate it, navigate to the ``benchmark`` directory and run the following command.
```
julia --project make.jl
```

## Chordal Completion

| library | name | vertices | edges | time | memory |
| :------ | :----| :------- | :---- | :--- | :----- |
| StructuredDecompositions | mycielskian4 | 11 | 23 | 3.229 μs | 7.45 KiB |
| QDLDL | mycielskian4 | 11 | 23 | 1.062 μs | 4.70 KiB |
| StructuredDecompositions | dwt_59 | 59 | 104 | 9.625 μs | 27.84 KiB |
| QDLDL | dwt_59 | 59 | 104 | 3.500 μs | 19.86 KiB |
| StructuredDecompositions | can_292 | 292 | 1124 | 53.042 μs | 162.38 KiB |
| QDLDL | can_292 | 292 | 1124 | 28.333 μs | 146.08 KiB |
| StructuredDecompositions | lshp3466 | 3466 | 10215 | 721.125 μs | 2.07 MiB |
| QDLDL | lshp3466 | 3466 | 10215 | 792.708 μs | 2.32 MiB |
| StructuredDecompositions | wing | 62032 | 121544 | 29.302 ms | 117.51 MiB |
| QDLDL | wing | 62032 | 121544 | 98.131 ms | 177.01 MiB |
| StructuredDecompositions | 144 | 144649 | 1074393 | 149.070 ms | 889.06 MiB |
| QDLDL | 144 | 144649 | 1074393 | 1.095 s | 1.47 GiB |
| StructuredDecompositions | 333SP | 3712815 | 11108633 | 1.441 s | 3.06 GiB |
| QDLDL | 333SP | 3712815 | 11108633 | 2.325 s | 3.89 GiB |

## Vertex Coloring

| library | bags | colors | time | memory |
| :------ | :--- | :----- | :----| :----- |
| StructuredDecompositions | 1 | 2 | 46.419 ms | 35.33 MiB |
| StructuredDecompositions | 2 | 2 | 29.778 ms | 21.42 MiB |
| StructuredDecompositions | 4 | 2 | 27.078 ms | 17.87 MiB |
| StructuredDecompositions | 8 | 2 | 28.743 ms | 17.41 MiB |
| StructuredDecompositions | 12 | 2 | 28.494 ms | 16.03 MiB |
| Catlab |     | 2 | 1.824 ms | 1.37 MiB |
| SimpleGraphAlgorithms |     | 2 | 919.375 μs | 65.93 KiB |
| GenericTensorNetworks |     | 2 | 756.833 μs | 439.58 KiB |
| StructuredDecompositions | 1 | 3 | 6.209 s | 3.57 GiB |
| StructuredDecompositions | 2 | 3 | 388.327 ms | 217.43 MiB |
| StructuredDecompositions | 4 | 3 | 43.685 ms | 26.93 MiB |
| StructuredDecompositions | 8 | 3 | 39.913 ms | 23.98 MiB |
| StructuredDecompositions | 12 | 3 | 38.812 ms | 21.99 MiB |
| Catlab |     | 3 | 28.411 ms | 21.09 MiB |
| SimpleGraphAlgorithms |     | 3 | 2.390 ms | 628.40 KiB |
| GenericTensorNetworks |     | 3 | 769.083 μs | 637.53 KiB |
