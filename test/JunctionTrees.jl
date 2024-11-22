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

order = JunctionTrees.Order(graph, CuthillMcKeeJL_RCM())
@test length(order) == 17

order = JunctionTrees.Order(graph, AMDJL_AMD())
@test length(order) == 17

order = JunctionTrees.Order(graph, MetisJL_ND())
@test length(order) == 17

order = JunctionTrees.Order(graph, TreeWidthSolverJL_BT())
@test length(order) == 17

order = JunctionTrees.Order(graph, MCS())
@test length(order) == 17

order = JunctionTrees.Order(1:17)
@test length(order) == 17

# Figure 4.3
jtree = JunctionTree(graph, order, Node())
@test width(jtree) == 4
@test treeheight(jtree) == 7
@test treesize(jtree) == 17

@test map(i -> parentindex(jtree, i), 1:17)  == [
    2,
    9,
    9,
    6,
    6,
    7,
    8,
    9,
    10,
    16,
    14,
    13,
    14,
    15,
    16,
    17,
    nothing,
]

@test map(i -> childindices(jtree, i), 1:17) == [
    [],
    [1],
    [],
    [],
    [],
    [4, 5],
    [6],
    [7],
    [2, 3, 8],
    [9],
    [],
    [],
    [12],
    [11, 13],
    [14],
    [10, 15],
    [16],
]

@test map(i -> residual(jtree, i), 1:17) == [
    [7],  # g
    [8],  # h
    [6],  # f
    [2],  # b
    [1],  # a
    [3],  # c
    [4],  # d
    [5],  # e
    [9],  # i
    [15], # o
    [12], # l
    [10], # j
    [11], # k
    [13], # m
    [14], # n
    [16], # p
    [17], # q
]

@test map(i -> seperator(jtree, i), 1:17) == [
    [8, 9, 15],       # h i o
    [9, 15],          # i o
    [9, 16],          # i p
    [3, 4],           # c d
    [3, 4, 5, 15],    # c d e o
    [4, 5, 15],       # d e o
    [5, 15],          # e o
    [9, 15, 16],      # i o p
    [15, 16],         # o p
    [16, 17],         # p q
    [13, 14, 16, 17], # m n p q
    [11, 13, 14, 17], # k m n q
    [13, 14, 17],     # m n q
    [14, 16, 17],     # n p q
    [16, 17],         # p q
    [17],             # q
    [],               #
] 


# Figure 4.7 (left)
jtree = JunctionTree(graph, order, Maximal())
@test width(jtree) == 4
@test treeheight(jtree) == 4
@test treesize(jtree) == 8

@test map(i -> parentindex(jtree, i), 1:8)  == [
    5,
    5,
    4,
    5,
    6,
    8,
    8,
    nothing
]

@test map(i -> childindices(jtree, i), 1:8)  == [
    [],
    [],
    [],
    [3],
    [1, 2, 4],
    [5],
    [],
    [6, 7],
]

@test map(i -> residual(jtree, i), 1:8) == [
    [7, 8],               # g h
    [6],                  # f
    [2],                  # b
    [1, 3, 4],            # a c d
    [5, 9],               # e i
    [15],                 # o
    [10, 11],             # j k
    [12, 13, 14, 16, 17], # l m n p q
]

@test map(i -> seperator(jtree, i), 1:8) == [
    [9, 15],      # i o
    [9, 16],      # i p
    [3, 4],       # c d
    [5, 15],      # e o
    [15, 16],     # o p
    [16, 17],     # p q
    [13, 14, 17], # m n q
    [],           #
]

# Figure 4.9
jtree = JunctionTree(graph, order, Fundamental())
@test width(jtree) == 4
@test treeheight(jtree) == 5
@test treesize(jtree) == 12

@test map(i -> parentindex(jtree, i), 1:12)  == [
    7,
    7,
    5,
    5,
    6,
    7,
    8,
    12,
    11,
    11,
    12,
    nothing,
]

@test map(i -> childindices(jtree, i), 1:12)  == [
    [],
    [],
    [],
    [],
    [3, 4],
    [5],
    [1, 2, 6],
    [7],
    [],
    [],
    [9, 10],
    [8, 11],
]

@test map(i -> residual(jtree, i), 1:12) == [
    [7, 8],   # g h
    [6],      # f
    [2],      # b
    [1],      # a
    [3, 4],   # c d
    [5],      # e
    [9],      # i
    [15],     # o
    [12],     # l
    [10, 11], # j k
    [13, 14], # m n
    [16, 17], # p q
]

@test map(i -> seperator(jtree, i), 1:12) == [
    [9, 15],          # i o
    [9, 16],          # i p
    [3, 4],           # c d
    [3, 4, 5, 15],    # c d e o
    [5, 15],          # e o
    [9, 15, 16],      # i o p
    [15, 16],         # o p
    [16, 17],         # p q
    [13, 14, 16, 17], # m n p q
    [13, 14, 17],     # m n q
    [16, 17],         # p q
    [],               #
]
