# Benchmarks

To regenerate this file, navigate to the ``benchmark`` directory and run the following command.
```
julia --project make.jl
```

## Junction Tree Construction

| library | supernode | name | edges | time | memory |
| :------ | :-------- | :----| :---- | :--- | :----- |
| StructuredDecompositions | maximal | mycielskian2 | 1 | 1.608 μs | 3.69 KiB |
| StructuredDecompositions | fundamental | mycielskian2 | 1 | 1.587 μs | 3.69 KiB |
| StructuredDecompositions | nodal | mycielskian2 | 1 | 1.379 μs | 3.34 KiB |
| QDLDL |      | mycielskian2 | 1 | 622.503 ns | 1.81 KiB |
| StructuredDecompositions | maximal | mycielskian4 | 23 | 2.449 μs | 6.70 KiB |
| StructuredDecompositions | fundamental | mycielskian4 | 23 | 2.625 μs | 6.97 KiB |
| StructuredDecompositions | nodal | mycielskian4 | 23 | 2.412 μs | 6.33 KiB |
| QDLDL |      | mycielskian4 | 23 | 1.029 μs | 4.70 KiB |
| StructuredDecompositions | maximal | dwt_59 | 104 | 7.490 μs | 25.78 KiB |
| StructuredDecompositions | fundamental | dwt_59 | 104 | 7.688 μs | 25.78 KiB |
| StructuredDecompositions | nodal | dwt_59 | 104 | 7.375 μs | 23.89 KiB |
| QDLDL |      | dwt_59 | 104 | 3.432 μs | 19.86 KiB |
| StructuredDecompositions | maximal | can_292 | 1124 | 47.500 μs | 146.36 KiB |
| StructuredDecompositions | fundamental | can_292 | 1124 | 47.709 μs | 148.39 KiB |
| StructuredDecompositions | nodal | can_292 | 1124 | 49.916 μs | 146.75 KiB |
| QDLDL |      | can_292 | 1124 | 28.958 μs | 146.08 KiB |
| StructuredDecompositions | maximal | lshp3466 | 10215 | 640.333 μs | 1.49 MiB |
| StructuredDecompositions | fundamental | lshp3466 | 10215 | 635.625 μs | 1.49 MiB |
| StructuredDecompositions | nodal | lshp3466 | 10215 | 882.875 μs | 2.38 MiB |
| QDLDL |      | lshp3466 | 10215 | 786.000 μs | 2.32 MiB |
| StructuredDecompositions | maximal | wing | 121544 | 18.089 ms | 28.59 MiB |
| StructuredDecompositions | fundamental | wing | 121544 | 19.034 ms | 29.89 MiB |
| StructuredDecompositions | nodal | wing | 121544 | 54.193 ms | 179.76 MiB |
| QDLDL |      | wing | 121544 | 97.048 ms | 177.01 MiB |
| StructuredDecompositions | maximal | 144 | 1074393 | 67.279 ms | 98.68 MiB |
| StructuredDecompositions | fundamental | 144 | 1074393 | 66.697 ms | 98.85 MiB |
| StructuredDecompositions | nodal | 144 | 1074393 | 458.174 ms | 1.45 GiB |
| QDLDL |      | 144 | 1074393 | 1.097 s | 1.47 GiB |
| StructuredDecompositions | maximal | 333SP | 11108633 | 1.214 s | 1.64 GiB |
| StructuredDecompositions | fundamental | 333SP | 11108633 | 1.240 s | 1.64 GiB |
| StructuredDecompositions | nodal | 333SP | 11108633 | 2.442 s | 3.95 GiB |
| QDLDL |      | 333SP | 11108633 | 2.728 s | 3.89 GiB |

## Vertex Coloring

| library | bags | colors | time | memory |
| :------ | :--- | :----- | :----| :----- |
| StructuredDecompositions | 1 | 2 | 4.099 ms | 3.07 MiB |
| StructuredDecompositions | 2 | 2 | 3.890 ms | 2.93 MiB |
| StructuredDecompositions | 4 | 2 | 5.371 ms | 3.95 MiB |
| StructuredDecompositions | 8 | 2 | 7.249 ms | 5.17 MiB |
| StructuredDecompositions | 12 | 2 | 8.049 ms | 5.57 MiB |
| Catlab |     | 2 | 1.758 ms | 1.37 MiB |
| SimpleGraphAlgorithms |     | 2 | 7.716 ns | 0 bytes |
| StructuredDecompositions | 1 | 3 | 128.516 s | 85.64 GiB |
| StructuredDecompositions | 2 | 3 | 1.474 s | 972.17 MiB |
| StructuredDecompositions | 4 | 3 | 234.185 ms | 168.18 MiB |
| StructuredDecompositions | 8 | 3 | 117.473 ms | 83.41 MiB |
| StructuredDecompositions | 12 | 3 | 64.865 ms | 46.25 MiB |
| Catlab |     | 3 | 27.626 ms | 21.09 MiB |
| SimpleGraphAlgorithms |     | 3 | 7.758 ns | 0 bytes |
