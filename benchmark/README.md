# Benchmarks

This file was automatically generated on 2025-01-19. To regenerate it, navigate to the ``benchmark`` directory and run the following command.
```
julia --project make.jl
```

## Junction Tree Construction

| library | supernode | name | edges | time | memory |
| :------ | :-------- | :----| :---- | :--- | :----- |
| StructuredDecompositions | maximal | mycielskian4 | 23 | 2.417 μs | 6.70 KiB |
| StructuredDecompositions | fundamental | mycielskian4 | 23 | 2.537 μs | 6.97 KiB |
| StructuredDecompositions | nodal | mycielskian4 | 23 | 2.333 μs | 6.33 KiB |
| QDLDL |      | mycielskian4 | 23 | 1.079 μs | 4.70 KiB |
| StructuredDecompositions | maximal | can_292 | 1124 | 47.625 μs | 146.36 KiB |
| StructuredDecompositions | fundamental | can_292 | 1124 | 48.667 μs | 148.39 KiB |
| StructuredDecompositions | nodal | can_292 | 1124 | 50.584 μs | 146.75 KiB |
| QDLDL |      | can_292 | 1124 | 28.083 μs | 146.08 KiB |
| StructuredDecompositions | maximal | wing | 121544 | 18.623 ms | 28.59 MiB |
| StructuredDecompositions | fundamental | wing | 121544 | 17.968 ms | 29.89 MiB |
| StructuredDecompositions | nodal | wing | 121544 | 55.779 ms | 179.76 MiB |
| QDLDL |      | wing | 121544 | 98.545 ms | 177.01 MiB |
| StructuredDecompositions | maximal | 333SP | 11108633 | 1.154 s | 1.64 GiB |
| StructuredDecompositions | fundamental | 333SP | 11108633 | 1.212 s | 1.64 GiB |
| StructuredDecompositions | nodal | 333SP | 11108633 | 1.865 s | 3.95 GiB |
| QDLDL |      | 333SP | 11108633 | 2.293 s | 3.89 GiB |

## Vertex Coloring

| library | bags | colors | time | memory |
| :------ | :--- | :----- | :----| :----- |
| StructuredDecompositions | 1 | 2 | 48.854 ms | 35.33 MiB |
| StructuredDecompositions | 2 | 2 | 30.523 ms | 21.42 MiB |
| StructuredDecompositions | 4 | 2 | 26.950 ms | 17.87 MiB |
| StructuredDecompositions | 8 | 2 | 28.435 ms | 17.41 MiB |
| StructuredDecompositions | 12 | 2 | 27.363 ms | 16.03 MiB |
| Catlab |     | 2 | 1.598 ms | 1.37 MiB |
| SimpleGraphAlgorithms |     | 2 | 1.067 ms | 65.93 KiB |
| GenericTensorNetworks |     | 2 | 759.791 μs | 446.41 KiB |
| StructuredDecompositions | 1 | 3 | 6.639 s | 3.57 GiB |
| StructuredDecompositions | 2 | 3 | 387.453 ms | 217.43 MiB |
| StructuredDecompositions | 4 | 3 | 43.337 ms | 26.93 MiB |
| StructuredDecompositions | 8 | 3 | 39.029 ms | 23.98 MiB |
| StructuredDecompositions | 12 | 3 | 37.165 ms | 21.99 MiB |
| Catlab |     | 3 | 25.115 ms | 21.09 MiB |
| SimpleGraphAlgorithms |     | 3 | 2.346 ms | 628.30 KiB |
| GenericTensorNetworks |     | 3 | 728.667 μs | 637.53 KiB |
