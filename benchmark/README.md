# Benchmarks

To regenerate this file, navigate to the ``benchmark`` directory and run the following command.
```
julia --project make.jl
```

## Junction Tree Construction

| library | name | supernode | edges | time | memory |
| :------ | :--- | :---------| :---- | :--- | :----- |
| StructuredDecompositions | mycielskian2 | maximal | 1 | 1.621 μs | 3.69 KiB |
| StructuredDecompositions | mycielskian2 | fundamental | 1 | 1.621 μs | 3.69 KiB |
| StructuredDecompositions | mycielskian2 | nodal | 1 | 1.371 μs | 3.34 KiB |
| QDLDL | mycielskian2 |      | 1 | 633.485 ns | 1.81 KiB |
| StructuredDecompositions | mycielskian4 | maximal | 23 | 2.431 μs | 6.70 KiB |
| StructuredDecompositions | mycielskian4 | fundamental | 23 | 2.699 μs | 6.97 KiB |
| StructuredDecompositions | mycielskian4 | nodal | 23 | 2.403 μs | 6.33 KiB |
| QDLDL | mycielskian4 |      | 23 | 1.038 μs | 4.70 KiB |
| StructuredDecompositions | dwt_59 | maximal | 104 | 7.511 μs | 25.78 KiB |
| StructuredDecompositions | dwt_59 | fundamental | 104 | 7.604 μs | 25.78 KiB |
| StructuredDecompositions | dwt_59 | nodal | 104 | 7.323 μs | 23.89 KiB |
| QDLDL | dwt_59 |      | 104 | 3.479 μs | 19.86 KiB |
| StructuredDecompositions | can_292 | maximal | 1124 | 47.166 μs | 146.36 KiB |
| StructuredDecompositions | can_292 | fundamental | 1124 | 47.667 μs | 148.39 KiB |
| StructuredDecompositions | can_292 | nodal | 1124 | 49.792 μs | 146.75 KiB |
| QDLDL | can_292 |      | 1124 | 28.500 μs | 146.08 KiB |
| StructuredDecompositions | lshp3466 | maximal | 10215 | 642.250 μs | 1.49 MiB |
| StructuredDecompositions | lshp3466 | fundamental | 10215 | 634.792 μs | 1.49 MiB |
| StructuredDecompositions | lshp3466 | nodal | 10215 | 895.125 μs | 2.38 MiB |
| QDLDL | lshp3466 |      | 10215 | 787.208 μs | 2.32 MiB |
| StructuredDecompositions | wing | maximal | 121544 | 18.058 ms | 28.59 MiB |
| StructuredDecompositions | wing | fundamental | 121544 | 18.351 ms | 29.89 MiB |
| StructuredDecompositions | wing | nodal | 121544 | 54.642 ms | 179.76 MiB |
| QDLDL | wing |      | 121544 | 97.489 ms | 177.01 MiB |
| StructuredDecompositions | 144 | maximal | 1074393 | 66.672 ms | 98.68 MiB |
| StructuredDecompositions | 144 | fundamental | 1074393 | 66.247 ms | 98.85 MiB |
| StructuredDecompositions | 144 | nodal | 1074393 | 478.761 ms | 1.45 GiB |
| QDLDL | 144 |      | 1074393 | 1.170 s | 1.47 GiB |
| StructuredDecompositions | 333SP | maximal | 11108633 | 1.198 s | 1.64 GiB |
| StructuredDecompositions | 333SP | fundamental | 11108633 | 1.205 s | 1.64 GiB |
| StructuredDecompositions | 333SP | nodal | 11108633 | 2.146 s | 3.95 GiB |
| QDLDL | 333SP |      | 11108633 | 2.458 s | 3.89 GiB |

## Vertex Coloring

| library | colors | bags | time | memory |
| :------ | :----- | :--- | :----| :----- |
| StructuredDecompositions | 2 | 1 | 4.179 ms | 3.07 MiB |
| StructuredDecompositions | 2 | 2 | 3.936 ms | 2.93 MiB |
| StructuredDecompositions | 2 | 4 | 5.452 ms | 3.95 MiB |
| StructuredDecompositions | 2 | 8 | 7.264 ms | 5.17 MiB |
| StructuredDecompositions | 2 | 12 | 8.072 ms | 5.57 MiB |
| Catlab | 2 |     | 1.789 ms | 1.37 MiB |
| SimpleGraphAlgorithms | 2 |     | 7.716 ns | 0 bytes |
| StructuredDecompositions | 3 | 1 | 129.200 s | 85.64 GiB |
| StructuredDecompositions | 3 | 2 | 1.510 s | 972.17 MiB |
| StructuredDecompositions | 3 | 4 | 241.629 ms | 168.18 MiB |
| StructuredDecompositions | 3 | 8 | 120.200 ms | 83.41 MiB |
| StructuredDecompositions | 3 | 12 | 65.602 ms | 46.25 MiB |
| Catlab | 3 |     | 28.123 ms | 21.09 MiB |
| SimpleGraphAlgorithms | 3 |     | 7.716 ns | 0 bytes |
