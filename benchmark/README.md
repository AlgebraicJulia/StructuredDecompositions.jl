# Benchmarks

This file was automatically generated on 2025-01-29. To regenerate it, navigate to the ``benchmark`` directory and run the following command.
```
julia --project make.jl
```

## Chordal Completion

| library | name | vertices | edges | time | memory |
| :------ | :----| :------- | :---- | :--- | :----- |
| StructuredDecompositions | mycielskian4 | 11 | 23 | 2.828 μs | 6.88 KiB |
| QDLDL | mycielskian4 | 11 | 23 | 1.096 μs | 4.70 KiB |
| StructuredDecompositions | dwt_59 | 59 | 104 | 7.792 μs | 25.45 KiB |
| QDLDL | dwt_59 | 59 | 104 | 3.505 μs | 19.86 KiB |
| StructuredDecompositions | can_292 | 292 | 1124 | 42.959 μs | 142.41 KiB |
| QDLDL | can_292 | 292 | 1124 | 28.458 μs | 146.08 KiB |
| StructuredDecompositions | lshp3466 | 3466 | 10215 | 582.166 μs | 1.93 MiB |
| QDLDL | lshp3466 | 3466 | 10215 | 793.417 μs | 2.32 MiB |
| StructuredDecompositions | wing | 62032 | 121544 | 23.325 ms | 113.82 MiB |
| QDLDL | wing | 62032 | 121544 | 99.133 ms | 177.01 MiB |
| StructuredDecompositions | 144 | 144649 | 1074393 | 138.158 ms | 872.25 MiB |
| QDLDL | 144 | 144649 | 1074393 | 1.084 s | 1.47 GiB |
| StructuredDecompositions | 333SP | 3712815 | 11108633 | 1.287 s | 2.89 GiB |
| QDLDL | 333SP | 3712815 | 11108633 | 2.582 s | 3.89 GiB |

## Vertex Coloring

| library | bags | colors | time | memory |
| :------ | :--- | :----- | :----| :----- |
| StructuredDecompositions | 1 | 2 | 45.973 ms | 35.33 MiB |
| StructuredDecompositions | 2 | 2 | 29.229 ms | 21.42 MiB |
| StructuredDecompositions | 4 | 2 | 26.571 ms | 17.87 MiB |
| StructuredDecompositions | 8 | 2 | 28.301 ms | 17.41 MiB |
| StructuredDecompositions | 12 | 2 | 27.890 ms | 16.03 MiB |
| Catlab |     | 2 | 1.790 ms | 1.37 MiB |
| SimpleGraphAlgorithms |     | 2 | 819.125 μs | 65.93 KiB |
| GenericTensorNetworks |     | 2 | 772.500 μs | 441.95 KiB |
| StructuredDecompositions | 1 | 3 | 6.441 s | 3.57 GiB |
| StructuredDecompositions | 2 | 3 | 395.514 ms | 217.43 MiB |
| StructuredDecompositions | 4 | 3 | 43.277 ms | 26.93 MiB |
| StructuredDecompositions | 8 | 3 | 39.771 ms | 23.98 MiB |
| StructuredDecompositions | 12 | 3 | 38.253 ms | 21.99 MiB |
| Catlab |     | 3 | 28.554 ms | 21.09 MiB |
| SimpleGraphAlgorithms |     | 3 | 2.543 ms | 628.40 KiB |
| GenericTensorNetworks |     | 3 | 785.208 μs | 635.42 KiB |
