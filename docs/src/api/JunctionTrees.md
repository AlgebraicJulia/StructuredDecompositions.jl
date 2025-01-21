# JunctionTrees

```@meta
CurrentModule = StructuredDecompositions.JunctionTrees
```

## Trees
```@docs
Tree
AbstractTree
eliminationtree
eliminationtree!
rootindex
parentindex
firstchildindex
nextsiblingindex
rootindices
childindices
```

## Junction Trees

```@docs
Bag
JunctionTree
junctiontree!
junctiontree
treewidth!
treewidth
separator
residual
relative
```

## Chordal Graphs
```@docs
filledgraph
ischordal
isfilled
isperfect
```

## Elimination Algorithms

```@docs
EliminationAlgorithm
PermutationOrAlgorithm
DEFAULT_ELIMINATION_ALGORITHM
MCS
RCM
AMD
SymAMD
MMD
NodeND
Spectral
BT
permutation   
```

## Supernodes
```@docs
SupernodeType
DEFAULT_SUPERNODE_TYPE
Nodal
Maximal
Fundamental
```
