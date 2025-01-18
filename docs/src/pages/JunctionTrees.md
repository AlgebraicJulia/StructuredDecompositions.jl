# [Junction Trees] (@id JunctionTrees)

JunctionTrees.jl is a Julia package for constructing [tree decompositions](https://en.wikipedia.org/wiki/Tree_decomposition) of simple graphs. You can use it as follows.

```julia
julia> using StructuredDecompositions

julia> graph = [
           0 1 1 0 0 0 0 0
           1 0 1 0 0 1 0 0
           1 1 0 1 1 0 0 0
           0 0 1 0 1 0 0 0
           0 0 1 1 0 0 1 1
           0 1 0 0 0 0 1 0
           0 0 0 0 1 1 0 1
           0 0 0 0 1 0 1 0
       ];

julia> label, tree = junctiontree(graph);

julia> tree
6-element JunctionTree:
[6, 7, 8]
├─ [1, 6, 7]
├─ [4, 6, 8]
│  └─ [3, 4, 6]
│     └─ [2, 3, 6]
└─ [5, 7, 8]
```

A junction tree is vector of bags, so you can retrieve the bag at node 3 by typing `tree[3]`.
```julia
julia> bag = tree[3]
3-element Bag:
 3
 4
 6
```

Notice that the bag is sorted. Its elements correspond to the vertices `label[bag]`.
```julia
julia> vertices = label[bag]
3-element Vector{Int64}:
 7
 6
 5
```

The width of a junction tree is computed by the function `treewidth`.
```julia
julia> treewidth(tree)
2
```
