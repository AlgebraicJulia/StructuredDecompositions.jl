# Benchmarks

This file was automatically generated on 2025-01-30. To regenerate it, navigate to the ``benchmark`` directory and run the following command.
```
julia --project make.jl
```

## Chordal Completion

| library | name | vertices | edges | time | memory |
| :------ | :----| :------- | :---- | :--- | :----- |
| StructuredDecompositions | mycielskian4 | 11 | 23 | 2.630 μs | 6.84 KiB |
| QDLDL | mycielskian4 | 11 | 23 | 1.113 μs | 4.70 KiB |
| StructuredDecompositions | dwt_59 | 59 | 104 | 7.625 μs | 25.58 KiB |
| QDLDL | dwt_59 | 59 | 104 | 3.609 μs | 19.86 KiB |
| StructuredDecompositions | can_292 | 292 | 1124 | 42.584 μs | 144.56 KiB |
| QDLDL | can_292 | 292 | 1124 | 26.709 μs | 146.08 KiB |
| StructuredDecompositions | lshp3466 | 3466 | 10215 | 573.084 μs | 1.95 MiB |
| QDLDL | lshp3466 | 3466 | 10215 | 792.958 μs | 2.32 MiB |
| StructuredDecompositions | wing | 62032 | 121544 | 23.451 ms | 114.05 MiB |
| QDLDL | wing | 62032 | 121544 | 98.276 ms | 177.01 MiB |
| StructuredDecompositions | 144 | 144649 | 1074393 | 131.527 ms | 874.30 MiB |
| QDLDL | 144 | 144649 | 1074393 | 1.113 s | 1.47 GiB |
| StructuredDecompositions | 333SP | 3712815 | 11108633 | 1.163 s | 2.91 GiB |
| QDLDL | 333SP | 3712815 | 11108633 | 2.424 s | 3.89 GiB |

## Vertex Coloring

| library | bags | colors | time | memory |
| :------ | :--- | :----- | :----| :----- |
| StructuredDecompositions | 1 | 2 | 47.227 ms | 35.33 MiB |
| StructuredDecompositions | 2 | 2 | 29.149 ms | 21.42 MiB |
| StructuredDecompositions | 4 | 2 | 26.474 ms | 17.87 MiB |
| StructuredDecompositions | 8 | 2 | 28.630 ms | 17.41 MiB |
| StructuredDecompositions | 12 | 2 | 29.182 ms | 16.03 MiB |
| Catlab |     | 2 | 1.840 ms | 1.37 MiB |
| SimpleGraphAlgorithms |     | 2 | 861.667 μs | 65.93 KiB |
| GenericTensorNetworks |     | 2 | 773.250 μs | 449.08 KiB |
| StructuredDecompositions | 1 | 3 | 6.361 s | 3.57 GiB |
| StructuredDecompositions | 2 | 3 | 394.797 ms | 217.43 MiB |
| StructuredDecompositions | 4 | 3 | 42.849 ms | 26.93 MiB |
| StructuredDecompositions | 8 | 3 | 39.746 ms | 23.98 MiB |
| StructuredDecompositions | 12 | 3 | 38.268 ms | 21.99 MiB |
| Catlab |     | 3 | 28.564 ms | 21.09 MiB |
| SimpleGraphAlgorithms |     | 3 | 2.562 ms | 628.40 KiB |
| GenericTensorNetworks |     | 3 | 798.834 μs | 650.52 KiB |
