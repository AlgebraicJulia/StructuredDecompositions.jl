# Benchmarks

To regenerate this file, navigate to the ``benchmark`` directory and run the following command.
```
julia --project make.jl
```

## Junction Tree Construction

| library | supernode | name | edges | time | memory |
| :------ | :-------- | :----| :---- | :--- | :----- |
| StructuredDecompositions | maximal | mycielskian4 | 23 | 2.426 μs | 6.70 KiB |
| StructuredDecompositions | fundamental | mycielskian4 | 23 | 2.565 μs | 6.97 KiB |
| StructuredDecompositions | nodal | mycielskian4 | 23 | 2.333 μs | 6.33 KiB |
| QDLDL |      | mycielskian4 | 23 | 1.046 μs | 4.70 KiB |
| StructuredDecompositions | maximal | can_292 | 1124 | 46.708 μs | 146.36 KiB |
| StructuredDecompositions | fundamental | can_292 | 1124 | 47.792 μs | 148.39 KiB |
| StructuredDecompositions | nodal | can_292 | 1124 | 50.042 μs | 146.75 KiB |
| QDLDL |      | can_292 | 1124 | 28.250 μs | 146.08 KiB |
| StructuredDecompositions | maximal | wing | 121544 | 17.985 ms | 28.59 MiB |
| StructuredDecompositions | fundamental | wing | 121544 | 17.644 ms | 29.89 MiB |
| StructuredDecompositions | nodal | wing | 121544 | 54.461 ms | 179.76 MiB |
| QDLDL |      | wing | 121544 | 97.953 ms | 177.01 MiB |
| StructuredDecompositions | maximal | 333SP | 11108633 | 1.206 s | 1.64 GiB |
| StructuredDecompositions | fundamental | 333SP | 11108633 | 1.197 s | 1.64 GiB |
| StructuredDecompositions | nodal | 333SP | 11108633 | 2.519 s | 3.95 GiB |
| QDLDL |      | 333SP | 11108633 | 2.572 s | 3.89 GiB |

## Vertex Coloring

| library | bags | colors | time | memory |
| :------ | :--- | :----- | :----| :----- |
| StructuredDecompositions | 1 | 2 | 53.562 ms | 35.33 MiB |
| StructuredDecompositions | 2 | 2 | 35.535 ms | 21.42 MiB |
| StructuredDecompositions | 4 | 2 | 31.310 ms | 17.87 MiB |
| StructuredDecompositions | 8 | 2 | 33.498 ms | 17.41 MiB |
| StructuredDecompositions | 12 | 2 | 32.965 ms | 16.03 MiB |
| Catlab |     | 2 | 1.919 ms | 1.37 MiB |
| SimpleGraphAlgorithms |     | 2 | 7.716 ns | 0 bytes |
| StructuredDecompositions | 1 | 3 | 7.147 s | 3.57 GiB |
| StructuredDecompositions | 2 | 3 | 419.376 ms | 217.43 MiB |
| StructuredDecompositions | 4 | 3 | 50.728 ms | 26.93 MiB |
| StructuredDecompositions | 8 | 3 | 46.970 ms | 23.98 MiB |
| StructuredDecompositions | 12 | 3 | 45.236 ms | 21.99 MiB |
| Catlab |     | 3 | 30.398 ms | 21.09 MiB |
| SimpleGraphAlgorithms |     | 3 | 7.716 ns | 0 bytes |
