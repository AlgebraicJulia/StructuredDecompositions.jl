using AbstractTrees
using Catlab.BasicGraphs
using Catlab.RelationalPrograms
using Catlab.UndirectedWiringDiagrams
using StructuredDecompositions.JunctionTrees
using Test


# Vandenberghe and Andersen
# Chordal Graphs and Semidefinite Optimization
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

order = JunctionTrees.Order(1:17)
@test length(order) == 17


# Figure 4.2
etree = EliminationTree(graph, order)
@test parentindex.([etree], 1:17) == [3, 3, 4, 5, 9, 9, 8, 9, 15, 11, 13, 13, 14, 16, 16, 17, 17]


# Figure 4.3
stree = SupernodeTree(etree, Node())


#=
@test supernode.([stree], order(stree, 1:17))



parent = JunctionTrees.makeetree(graph, order)

# Figure 4.2
@test parent == [3, 3, 4, 5, 9, 9, 8, 9, 15, 11, 13, 13, 14, 16, 16, 17, 0]

etree = JunctionTrees.Tree(17, parent)
indegrees, outdegrees = JunctionTrees.getdegrees(graph, order, etree)

@test indegrees == [0, 0, 2, 3, 3, 0, 0, 1, 4, 0, 1, 0, 3, 4, 7, 7, 7]
@test outdegrees == [4, 2, 3, 2, 3, 2, 3, 2, 2, 4, 3, 4, 3, 2, 2, 1, 0]

# Figure 4.3 
snd, sbt, q, a = JunctionTrees.makestree(etree, outdegrees, Node()) 

@test snd == [
    [1],
    [2],
    [3],
    [4],
    [5],
    [6],
    [7],
    [8],
    [9],
    [10],
    [11],
    [12],
    [13],
    [14],
    [15],
    [16],
    [17]]

@test sbt == 1:17
@test q == [3, 3, 4, 5, 9, 9, 8, 9, 15, 11, 13, 13, 14, 16, 16, 17, 0]
@test a == [3, 3, 4, 5, 9, 9, 8, 9, 15, 11, 13, 13, 14, 16, 16, 17, 0]

# Figure 4.7 (left)
snd, sbt, q, a = JunctionTrees.makestree(etree, outdegrees, MaximalSupernode()) 

@test snd == [
    [1, 3, 4],
    [2],
    [5, 9],
    [6],
    [7, 8],
    [10, 11],
    [12, 13, 14, 16, 17],
    [15] ]

@test sbt == [1, 2, 1, 1, 3, 4, 5, 5, 3, 6, 6, 7, 7, 7, 8, 7, 7]
@test q == [3, 1, 8, 3, 3, 7, 0, 7]
@test a == [5, 3, 15, 9, 9, 13, 0, 16]

# Figure 4.9
snd, sbt, q, a = JunctionTrees.makestree(etree, outdegrees, FundamentalSupernode()) 

@test snd == [
    [1],
    [2],
    [3, 4],
    [5],
    [6],
    [7, 8],
    [9],
    [10, 11],
    [12],
    [13, 14],
    [15],
    [16, 17] ]

@test sbt == [1, 2, 3, 3, 4, 5, 6, 6, 7, 8, 8, 9, 10, 10, 11, 12, 12]
@test q == [3, 3, 4, 7, 7, 7, 11, 10, 10, 12, 12, 0]
@test a == [3, 3, 5, 9, 9, 9, 15, 13, 13, 16, 16, 0]

# Figure 4.3
jtree = JunctionTree(graph, order, Node())

@test getresidual.([jtree], getsubtree.([jtree], 1:17)) == [
    [1],
    [2],
    [3],
    [4],
    [5],
    [6],
    [7],
    [8],
    [9],
    [10],
    [11],
    [12],
    [13],
    [14],
    [15],
    [16],
    [17] ]

@test getseperator.([jtree], getsubtree.([jtree], 1:17)) == [
    [3, 4, 5, 15],
    [3, 4],
    [4, 5, 15],
    [5, 15],
    [9, 15, 16],
    [9, 16],
    [8, 9, 15],
    [9, 15],
    [15, 16],
    [11, 13, 14, 17],
    [13, 14, 17],
    [13, 14, 16, 17],
    [14, 16, 17],
    [16, 17],
    [16, 17],
    [17],
    [] ]

@test getlevel.([jtree], getsubtree.([jtree], 1:17)) == [
    7,
    7,
    6,
    5,
    4,
    4,
    5,
    4,
    3,
    5,
    4,
    4,
    3,
    2,
    2,
    1,
    0 ]

@test isdescendant(jtree, getsubtree(jtree, 5), getsubtree(jtree, 15))
@test !isdescendant(jtree, getsubtree(jtree, 15), getsubtree(jtree, 5))
@test !isdescendant(jtree, getsubtree(jtree, 10), getsubtree(jtree, 15))
@test !isdescendant(jtree, getsubtree(jtree, 15), getsubtree(jtree, 10))
@test !isdescendant(jtree, getsubtree(jtree, 1), getsubtree(jtree, 1))
@test getwidth(jtree) == 4

# Figure 4.7 (left)
jtree = JunctionTree(graph, order, MaximalSupernode())

@test getresidual.([jtree], getsubtree.([jtree], 1:17)) == [
    [1, 3, 4],
    [2],
    [1, 3, 4],
    [1, 3, 4],
    [5, 9],
    [6],
    [7, 8],
    [7, 8],
    [5, 9],
    [10, 11],
    [10, 11],
    [12, 13, 14, 16, 17],
    [12, 13, 14, 16, 17],
    [12, 13, 14, 16, 17],
    [15],
    [12, 13, 14, 16, 17],
    [12, 13, 14, 16, 17]]

@test getseperator.([jtree], getsubtree.([jtree], 1:17)) == [
    [5, 15],
    [3, 4],
    [5, 15],
    [5, 15],
    [15, 16],
    [9, 16],
    [9, 15],
    [9, 15],
    [15, 16],
    [13, 14, 17],
    [13, 14, 17],
    [],
    [],
    [],
    [16, 17],
    [],
    []]

@test getlevel.([jtree], getsubtree.([jtree], 1:17)) == [
    3,
    4,
    3,
    3,
    2,
    3,
    3,
    3,
    2,
    1,
    1,
    0,
    0,
    0,
    1,
    0,
    0 ]

@test isdescendant(jtree, getsubtree(jtree, 5), getsubtree(jtree, 15))
@test !isdescendant(jtree, getsubtree(jtree, 15), getsubtree(jtree, 5))
@test !isdescendant(jtree, getsubtree(jtree, 10), getsubtree(jtree, 15))
@test !isdescendant(jtree, getsubtree(jtree, 15), getsubtree(jtree, 10))
@test !isdescendant(jtree, getsubtree(jtree, 1), getsubtree(jtree, 1))
@test getwidth(jtree) == 4

# Figure 4.9
jtree = JunctionTree(graph, order, FundamentalSupernode())

@test getresidual.([jtree], getsubtree.([jtree], 1:17)) == [
    [1],
    [2],
    [3, 4],
    [3, 4],
    [5],
    [6],
    [7, 8],
    [7, 8],
    [9],
    [10, 11],
    [10, 11],
    [12],
    [13, 14],
    [13, 14],
    [15],
    [16, 17],
    [16, 17]]

@test getseperator.([jtree], getsubtree.([jtree], 1:17)) == [
    [3, 4, 5, 15],
    [3, 4],
    [5, 15],
    [5, 15],
    [9, 15, 16],
    [9, 16],
    [9, 15],
    [9, 15],
    [15, 16],
    [13, 14, 17],
    [13, 14, 17],
    [13, 14, 16, 17],
    [16, 17],
    [16, 17],
    [16, 17],
    [],
    []]

@test isdescendant(jtree, getsubtree(jtree, 5), getsubtree(jtree, 15))
@test !isdescendant(jtree, getsubtree(jtree, 15), getsubtree(jtree, 5))
@test !isdescendant(jtree, getsubtree(jtree, 10), getsubtree(jtree, 15))
@test !isdescendant(jtree, getsubtree(jtree, 15), getsubtree(jtree, 10))
@test !isdescendant(jtree, getsubtree(jtree, 1), getsubtree(jtree, 1))
@test getwidth(jtree) == 4
=#
