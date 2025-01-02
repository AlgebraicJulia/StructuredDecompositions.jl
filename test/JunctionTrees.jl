using StructuredDecompositions.JunctionTrees


using AbstractTrees
using Catlab.BasicGraphs
using Catlab.RelationalPrograms
using Catlab.UndirectedWiringDiagrams
using Test


# Chordal Graphs and Semidefinite Optimization
# Vandenberghe and Andersen
graph = SymmetricGraph(17)

add_edges!(graph,
    [1, 1, 1, 1,  2, 2, 5, 5,  6, 6,  7, 7, 7,  10, 10, 10, 10, 12, 12, 12, 12, 15],
    [3, 4, 5, 15, 3, 4, 9, 16, 9, 16, 8, 9, 15, 11, 13, 14, 17, 13, 14, 16, 17, 17])

order = Permutation(graph, CuthillMcKeeJL_RCM())
@test length(order) == 17

order = Permutation(graph, SymRCMJL_RCM())
@test length(order) == 17

order = Permutation(graph, AMDJL_AMD())
@test length(order) == 17

order = Permutation(graph, AMDJL_SYMAMD())
@test length(order) == 17

order = Permutation(graph, MetisJL_ND())
@test length(order) == 17

order = Permutation(graph, TreeWidthSolverJL_BT())
@test length(order) == 17

order = Permutation(graph, MCS())
@test length(order) == 17

order = Permutation(1:17)
@test length(order) == 17

# Figure 4.3
permutation, tree = jtree(graph, order, Node())
@test treewidth(tree) == 4

@test map(i -> parentindex(tree, i), 1:17)  == [
    2,
    4,
    4,
    5,
    16,
    8,
    8,
    9,
    10,
    14,
    14,
    13,
    14,
    15,
    16,
    17,
    nothing,
]

@test map(i -> collect(childindices(tree, i)), 1:17) == [
    [],
    [1],
    [],
    [2, 3],
    [4],
    [],
    [],
    [6, 7],
    [8],
    [9],
    [],
    [],
    [12],
    [10, 11, 13],
    [14],
    [5, 15],
    [16],
]

@test map(i -> view(permutation, residual(tree, i)), 1:17) == [
    [10], # j
    [11], # k
    [12], # l
    [13], # m
    [14], # n
    [1],  # a
    [2],  # b
    [3],  # c
    [4],  # d
    [5],  # e
    [6],  # f
    [7],  # g
    [8],  # h
    [9],  # i
    [15], # o
    [16], # p
    [17], # q
]

@test map(i -> view(permutation, seperator(tree, i)), 1:17) == [
    [11, 13, 14, 17], # k m n q
    [13, 14, 17],     # m n q
    [13, 14, 16, 17], # m n p q
    [14, 16, 17],     # n p q
    [16, 17],         # p q
    [3, 4, 5, 15],    # c d e o
    [3, 4],           # c d
    [4, 5, 15],       # d e o
    [5, 15],          # e o
    [9, 15, 16],      # i o p
    [9, 16],          # i p
    [8, 9, 15],       # h i o
    [9, 15],          # i o
    [15, 16],         # o p
    [16, 17],         # p q
    [17],             # q
    [],               #
] 


# Figure 4.7 (left)
permutation, tree = jtree(graph, order, Maximal())
@test treewidth(tree) == 4

@test map(i -> parentindex(tree, i), 1:8)  == [
    8,
    3,
    6,
    6,
    6,
    7,
    8,
    nothing
]

@test map(i -> collect(childindices(tree, i)), 1:8)  == [
    [],
    [],
    [2],
    [],
    [],
    [3, 4, 5],
    [6],
    [1, 7],
]

@test map(i -> view(permutation, residual(tree, i)), 1:8) == [
    [10, 11],             # j k
    [2],                  # b
    [1, 3, 4],            # a c d
    [6],                  # f
    [7, 8],               # g h
    [5, 9],               # e i
    [15],                 # o
    [12, 13, 14, 16, 17], # l m n p q
]

@test map(i -> view(permutation, seperator(tree, i)), 1:8) == [
    [13, 14, 17], # m n q
    [3, 4],       # c d
    [5, 15],      # e o
    [9, 16],      # i p
    [9, 15],      # i o
    [15, 16],     # o p
    [16, 17],     # p q
    [],           #
]

# Figure 4.9
permutation, tree = jtree(graph, order, Fundamental())
@test treewidth(tree) == 4

@test map(i -> parentindex(tree, i), 1:12)  == [
    3,
    3,
    12,
    6,
    6,
    7,
    10,
    10,
    10,
    11,
    12,
    nothing,
]

@test map(i -> collect(childindices(tree, i)), 1:12)  == [
    [],
    [],
    [1, 2],
    [],
    [],
    [4, 5],
    [6],
    [],
    [],
    [7, 8, 9],
    [10],
    [3, 11],
]

@test map(i -> view(permutation, residual(tree, i)), 1:12) == [
    [10, 11], # j k
    [12],     # l
    [13, 14], # m n
    [1],      # a
    [2],      # b
    [3, 4],   # c d
    [5],      # e
    [6],      # f
    [7, 8],   # g h
    [9],      # i
    [15],     # o
    [16, 17], # p q
]

@test map(i -> view(permutation, seperator(tree, i)), 1:12) == [
    [13, 14, 17],     # m n q
    [13, 14, 16, 17], # m n p q
    [16, 17],         # p q
    [3, 4, 5, 15],    # c d e o
    [3, 4],           # c d
    [5, 15],          # e o
    [9, 15, 16],      # i o p
    [9, 16],          # i p
    [9, 15],          # i o
    [15, 16],         # o p
    [16, 17],         # p q
    [],               #
]
