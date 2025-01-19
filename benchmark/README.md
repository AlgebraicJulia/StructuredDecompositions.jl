# Benchmarks

This file was automatically generated on 2025-01-19. To regenerate it, navigate to the ``benchmark`` directory and run the following command.
```
julia --project make.jl
```

## Junction Tree Construction

| library | supernode | name | edges | time | memory |
| :------ | :-------- | :----| :---- | :--- | :----- |
| StructuredDecompositions | maximal | mycielskian4 | 23 | 2.347 μs | 6.70 KiB |
| StructuredDecompositions | fundamental | mycielskian4 | 23 | 2.505 μs | 6.97 KiB |
| StructuredDecompositions | nodal | mycielskian4 | 23 | 2.301 μs | 6.33 KiB |
| QDLDL |      | mycielskian4 | 23 | 1.062 μs | 4.70 KiB |
| StructuredDecompositions | maximal | can_292 | 1124 | 46.792 μs | 146.36 KiB |
| StructuredDecompositions | fundamental | can_292 | 1124 | 47.000 μs | 148.39 KiB |
| StructuredDecompositions | nodal | can_292 | 1124 | 50.083 μs | 146.75 KiB |
| QDLDL |      | can_292 | 1124 | 28.541 μs | 146.08 KiB |
| StructuredDecompositions | maximal | wing | 121544 | 18.299 ms | 28.59 MiB |
| StructuredDecompositions | fundamental | wing | 121544 | 18.617 ms | 29.89 MiB |
| StructuredDecompositions | nodal | wing | 121544 | 54.528 ms | 179.76 MiB |
| QDLDL |      | wing | 121544 | 97.459 ms | 177.01 MiB |
| StructuredDecompositions | maximal | 333SP | 11108633 | 1.227 s | 1.64 GiB |
| StructuredDecompositions | fundamental | 333SP | 11108633 | 1.190 s | 1.64 GiB |
| StructuredDecompositions | nodal | 333SP | 11108633 | 1.991 s | 3.95 GiB |
| QDLDL |      | 333SP | 11108633 | 2.704 s | 3.89 GiB |

## Vertex Coloring

| library | bags | colors | time | memory |
| :------ | :--- | :----- | :----| :----- |
| StructuredDecompositions | 1 | 2 | 48.644 ms | 35.33 MiB |
| StructuredDecompositions | 2 | 2 | 40.179 ms | 21.42 MiB |
| StructuredDecompositions | 4 | 2 | 27.130 ms | 17.87 MiB |
| StructuredDecompositions | 8 | 2 | 28.679 ms | 17.41 MiB |
| StructuredDecompositions | 12 | 2 | 27.599 ms | 16.03 MiB |
| Catlab |     | 2 | 1.584 ms | 1.37 MiB |
| SimpleGraphAlgorithms |     | 2 | 953.708 μs | 65.93 KiB |
| GenericTensorNetworks |     | 2 | 756.084 μs | 441.59 KiB |
| StructuredDecompositions | 1 | 3 | 6.244 s | 3.57 GiB |
| StructuredDecompositions | 2 | 3 | 394.211 ms | 217.43 MiB |
| StructuredDecompositions | 4 | 3 | 43.305 ms | 26.93 MiB |
| StructuredDecompositions | 8 | 3 | 39.150 ms | 23.98 MiB |
| StructuredDecompositions | 12 | 3 | 37.484 ms | 21.99 MiB |
| Catlab |     | 3 | 24.795 ms | 21.09 MiB |
| SimpleGraphAlgorithms |     | 3 | 2.203 ms | 628.30 KiB |
| GenericTensorNetworks |     | 3 | 733.042 μs | 638.47 KiB |
