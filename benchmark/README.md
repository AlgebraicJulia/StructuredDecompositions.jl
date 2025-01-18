# Benchmarks

This file was automatically generated on 2025-01-17. To regenerate it, navigate to the ``benchmark`` directory and run the following command.
```
julia --project make.jl
```

## Junction Tree Construction

| library | supernode | name | edges | time | memory |
| :------ | :-------- | :----| :---- | :--- | :----- |
| StructuredDecompositions | maximal | mycielskian4 | 23 | 2.491 μs | 6.70 KiB |
| StructuredDecompositions | fundamental | mycielskian4 | 23 | 2.620 μs | 6.97 KiB |
| StructuredDecompositions | nodal | mycielskian4 | 23 | 2.444 μs | 6.33 KiB |
| QDLDL |      | mycielskian4 | 23 | 1.033 μs | 4.70 KiB |
| StructuredDecompositions | maximal | can_292 | 1124 | 46.708 μs | 146.36 KiB |
| StructuredDecompositions | fundamental | can_292 | 1124 | 47.291 μs | 148.39 KiB |
| StructuredDecompositions | nodal | can_292 | 1124 | 49.542 μs | 146.75 KiB |
| QDLDL |      | can_292 | 1124 | 28.000 μs | 146.08 KiB |
| StructuredDecompositions | maximal | wing | 121544 | 18.023 ms | 28.59 MiB |
| StructuredDecompositions | fundamental | wing | 121544 | 18.042 ms | 29.89 MiB |
| StructuredDecompositions | nodal | wing | 121544 | 55.980 ms | 179.76 MiB |
| QDLDL |      | wing | 121544 | 96.785 ms | 177.01 MiB |
| StructuredDecompositions | maximal | 333SP | 11108633 | 1.164 s | 1.64 GiB |
| StructuredDecompositions | fundamental | 333SP | 11108633 | 1.199 s | 1.64 GiB |
| StructuredDecompositions | nodal | 333SP | 11108633 | 1.945 s | 3.95 GiB |
| QDLDL |      | 333SP | 11108633 | 2.307 s | 3.89 GiB |

## Vertex Coloring

| library | bags | colors | time | memory |
| :------ | :--- | :----- | :----| :----- |
| StructuredDecompositions | 1 | 2 | 52.604 ms | 35.33 MiB |
| StructuredDecompositions | 2 | 2 | 34.256 ms | 21.42 MiB |
| StructuredDecompositions | 4 | 2 | 30.434 ms | 17.87 MiB |
| StructuredDecompositions | 8 | 2 | 32.644 ms | 17.41 MiB |
| StructuredDecompositions | 12 | 2 | 32.090 ms | 16.03 MiB |
| Catlab |     | 2 | 1.934 ms | 1.37 MiB |
| SimpleGraphAlgorithms |     | 2 | 34.944 ms | 1.69 MiB |
| StructuredDecompositions | 1 | 3 | 7.584 s | 3.57 GiB |
| StructuredDecompositions | 2 | 3 | 417.496 ms | 217.43 MiB |
| StructuredDecompositions | 4 | 3 | 49.003 ms | 26.93 MiB |
| StructuredDecompositions | 8 | 3 | 45.527 ms | 23.98 MiB |
| StructuredDecompositions | 12 | 3 | 43.611 ms | 21.99 MiB |
| Catlab |     | 3 | 30.607 ms | 21.09 MiB |
| SimpleGraphAlgorithms |     | 3 | 35.462 ms | 1.69 MiB |
