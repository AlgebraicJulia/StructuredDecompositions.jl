# This file examines the relationship between run time, bag size, and number of bags for a fixed 40 node graph.
# Note since bags must overlap, the colorings_fixed will need to deal with more than 40 total nodes due to overlap.
# We will compare the results against running the graph coloring algorithm on Graphs.jl

# We will consider cases with 2(even), 4(even), 8(uneven), 16(uneven) bags
# We will first benchmark the coloring algorithm in graphs.jl


colors = (2, 3)


struct Coloring
    n::Int    
end


function K(n::Integer)
    complete_graph(Graph, n)
end


function (coloring::Coloring)(graph::Graph)
    FinSet(homomorphisms(graph, K(coloring.n); alg=HomomorphismQuery()))
end


function (coloring::Coloring)(f::ACSetTransformation)  
    FinFunction(λ -> compose(f, λ), coloring(codom(f)), coloring(dom(f)))
end
 
 
function skeletal_coloring(n::Integer)
    skeleton ∘ Coloring(n)
end 


graph = @acset Graph begin
    V = 40
    E = 63
    src = [1, 1, 2, 2, 3, 4, 4, 5, 5, 6, 7, 7, 8, 8, 9, 10, 10, 10, 10, 10, 10, 11, 11, 12, 13, 13, 14, 14, 15, 16, 16, 17, 17, 18, 19,
          21, 21, 22, 23, 23, 24, 24, 25, 26, 26, 27, 27, 28, 29, 31, 31, 32, 33, 33, 34, 34, 35, 36, 36, 37, 37, 38, 39]
    tgt = [2, 3, 3, 4, 4, 5, 6, 6, 7, 7, 8, 9, 9, 10, 10, 11, 12, 21, 22, 31, 32, 12, 13, 13, 14, 15, 15, 16, 16, 17, 18, 18, 19, 19, 20,
          22, 23, 23, 24, 25, 25, 26, 26, 27, 28, 28, 29, 29, 30, 32, 33, 33, 34, 35, 35, 36, 36, 37, 38, 38, 39, 39, 40]
end


for i in colors
    target = K(i)
    SUITE["graph coloring fixed"]["$i coloring"]["Catlab"] = @benchmarkable is_homomorphic($graph, $target)
end


simple = UndirectedGraph(adjacency_matrix(graph))


for i in colors
    SUITE["graph coloring fixed"]["$i coloring"]["SimpleGraphAlgorithms"] = @benchmarkable chromatic_number($simple)
end


# 1 bag case
# 40 nodes, no adhesions

# bag 1
H01 = graph

H02 = @acset Graph begin
    V = 1
    E = 0
    src = []
    tgt = []
end

H01_02 = @acset Graph begin
    V = 1
    E = 0
    src = []
    tgt = []
end

G = @acset Graph begin
    V = 2
    E = 1
    src = [1]
    tgt = [1]
end

# transformations
Γ = FinDomFunctor(
    Dict(1 => H01, 2 => H02, 3 => H01_02),
    Dict(
        1 => ACSetTransformation(H01_02, H01, V=[1]),
        2 => ACSetTransformation(H01_02, H02, V=[1])),
    ∫(G))


decomp = StrDecomp(G, Γ)

for i in colors
    SUITE["graph coloring fixed"]["$i coloring"]["StructuredDecompositions"]["1 bags"] = @benchmarkable decide_sheaf_tree_shape($(skeletal_coloring(i)), $decomp)
end

# 2 bag case 
# 21 nodes each, 1 node adhesion

# bag 1
H01 = @acset Graph begin
    V = 21
    E = 31
    src = [1, 1, 2, 2, 3, 4, 4, 5, 5, 6, 7, 7, 8, 8, 9, 10, 10, 11, 11, 12, 13, 13, 14, 14, 15, 16, 16, 17, 17, 18, 19]
    tgt = [2, 3, 3, 4, 4, 5, 6, 6, 7, 7, 8, 9, 9, 10, 10, 11, 12, 12, 13, 13, 14, 15, 15, 16, 16, 17, 18, 18, 19, 19, 20]
end

# adhesion 1, 2
H01_02 = @acset Graph begin
    V = 1
end

# bag 2
H02 = @acset Graph begin
    V = 21
    E = 31
    src = [1, 1, 2, 2, 3, 4, 4, 5, 5, 6, 7, 7, 8, 8, 9, 10, 10, 11, 11, 12, 13, 13, 14, 14, 15, 16, 16, 17, 17, 18, 19]
    tgt = [2, 3, 3, 4, 4, 5, 6, 6, 7, 7, 8, 9, 9, 10, 10, 11, 12, 12, 13, 13, 14, 15, 15, 16, 16, 17, 18, 18, 19, 19, 20]
end

G = @acset Graph begin
    V = 2
    E = 1
    src = [1]
    tgt = [1]
end

# transformations
Γ = FinDomFunctor(
    Dict(1 => H01, 2 => H02, 3 => H01_02),
    Dict(
        1 => ACSetTransformation(H01_02, H01, V=[1]),
        2 => ACSetTransformation(H01_02, H02, V=[1])),
    ∫(G))


decomp = StrDecomp(G, Γ)

for i in colors
    SUITE["graph coloring fixed"]["$i coloring"]["StructuredDecompositions"]["2 bags"] = @benchmarkable decide_sheaf_tree_shape($(skeletal_coloring(i)), $decomp)
end


# 4 bag case
# 1 bag w/ 10 nodes, 3 bags w/ 11 nodes, 1 node adhesion

# bag 1
H01 = @acset Graph begin
    V = 10
    E = 15
    src = [1, 1, 2, 2, 3, 4, 4, 5, 5, 6, 7, 7, 8, 8, 9]
    tgt = [2, 3, 3, 4, 4, 5, 6, 6, 7, 7, 8, 9, 9, 10, 10]
end

# adhesion 1, 2
H01_02 = @acset Graph begin
    V = 1
end

# bag 2
H02 = @acset Graph begin
    V = 11
    E = 16
    src = [1, 1, 2, 2, 3, 4, 4, 5, 5, 6, 7, 7, 8, 8, 9, 10]
    tgt = [2, 3, 3, 4, 4, 5, 6, 6, 7, 7, 8, 9, 9, 10, 10, 11]
end

# adhesion 3, 4
H02_03 = @acset Graph begin
    V = 1
end

# bag 3
H03 = @acset Graph begin
    V = 11
    E = 16
    src = [1, 1, 2, 2, 3, 4, 4, 5, 5, 6, 7, 7, 8, 8, 9, 10]
    tgt = [2, 3, 3, 4, 4, 5, 6, 6, 7, 7, 8, 9, 9, 10, 10, 11]
end

# adhesion 3, 4
H03_04 = @acset Graph begin
    V = 1
end

# bag 4
H04 = @acset Graph begin
    V = 11
    E = 16
    src = [1, 1, 2, 2, 3, 4, 4, 5, 5, 6, 7, 7, 8, 8, 9, 10]
    tgt = [2, 3, 3, 4, 4, 5, 6, 6, 7, 7, 8, 9, 9, 10, 10, 11]
end

G = @acset Graph begin
    V = 4
    E = 3
    src = [1, 2, 3]
    tgt = [2, 3, 4]
end

# transformations
Γ = FinDomFunctor(
    Dict(1 => H01, 2 => H02, 3 => H03, 4 => H04, 5 => H01_02, 6 => H02_03, 7 => H03_04),
    Dict(
        1 => ACSetTransformation(H01_02, H01, V=[1]),
        2 => ACSetTransformation(H01_02, H02, V=[1]),
        3 => ACSetTransformation(H02_03, H02, V=[1]),
        4 => ACSetTransformation(H02_03, H03, V=[1]),
        5 => ACSetTransformation(H03_04, H03, V=[1]),
        6 => ACSetTransformation(H03_04, H04, V=[1])),
    ∫(G))


decomp = StrDecomp(G, Γ)

for i in colors
    SUITE["graph coloring fixed"]["$i coloring"]["StructuredDecompositions"]["4 bags"] = @benchmarkable decide_sheaf_tree_shape($(skeletal_coloring(i)), $decomp)
end

# 8 bag case 
# bags between 4-7 nodes, 47 total nodes, 7 total adhesion nodes)

# bag 1
H01 = @acset Graph begin
    V = 4
    E = 5
    src = [1, 1, 2, 2, 3]
    tgt = [2, 3, 3, 4, 4]
end

# adhesion 1, 2
H01_02 = @acset Graph begin
    V = 1
end

# bag 2
H02 = @acset Graph begin
    V = 7
    E = 10
    src = [1, 1, 2, 2, 3, 4, 4, 5, 5, 6]
    tgt = [2, 3, 3, 4, 4, 5, 6, 6, 7, 7]
end

# adhesion 3, 4
H02_03 = @acset Graph begin
    V = 1
end

# bag 3
H03 = @acset Graph begin
    V = 7
    E = 10
    src = [1, 1, 2, 2, 3, 4, 4, 5, 5, 6]
    tgt = [2, 3, 3, 4, 4, 5, 6, 6, 7, 7]
end

# adhesion 3, 4
H03_04 = @acset Graph begin
    V = 1
end

# bag 4
H04 = @acset Graph begin
    V = 5
    E = 6
    src = [1, 1, 2, 2, 3, 4]
    tgt = [2, 3, 3, 4, 4, 5]
end

# adhesion 2, 5
H02_05 = @acset Graph begin
    V = 1
end

# bag 5
H05 = @acset Graph begin
    V = 7
    E = 10
    src = [1, 1, 2, 2, 3, 4, 4, 5, 5, 6]
    tgt = [2, 3, 3, 4, 4, 5, 6, 6, 7, 7]
end

# adhesion 5, 6
H05_06 = @acset Graph begin
    V = 1
end

# bag 6
H06 = @acset Graph begin
    V = 5
    E = 6
    src = [1, 1, 2, 2, 3, 4]
    tgt = [2, 3, 3, 4, 4, 5]
end

# adhesion 2, 7
H02_07 = @acset Graph begin
    V = 1
end

# bag 7
H07 = @acset Graph begin
    V = 7
    E = 10
    src = [1, 1, 2, 2, 3, 4, 4, 5, 5, 6]
    tgt = [2, 3, 3, 4, 4, 5, 6, 6, 7, 7]
end

# adhesion 7, 8
H07_08 = @acset Graph begin
    V = 1
end

# bag 8
H08 = @acset Graph begin
    V = 5
    E = 6
    src = [1, 1, 2, 2, 3, 4]
    tgt = [2, 3, 3, 4, 4, 5]
end

G = @acset Graph begin
    V = 8
    E = 7
    src = [1, 2, 2, 2, 3, 5, 7]
    tgt = [2, 3, 5, 7, 4, 6, 8]
end

# transformations
Γ = FinDomFunctor(
    Dict(1 => H01, 2 => H02, 3 => H03, 4 => H04, 5 => H05, 6 => H06, 7 => H07, 8 => H08,
         9 => H01_02, 10 => H02_03, 11 => H03_04, 12 => H02_05, 13 => H05_06, 14 => H02_07, 15 => H07_08),
    Dict(
        1  => ACSetTransformation(H01_02, H01, V=[4]),
        2  => ACSetTransformation(H01_02, H02, V=[1]),
        3  => ACSetTransformation(H02_03, H02, V=[7]),
        4  => ACSetTransformation(H02_03, H03, V=[1]),
        5  => ACSetTransformation(H03_04, H03, V=[7]),
        6  => ACSetTransformation(H03_04, H04, V=[1]),
        7  => ACSetTransformation(H02_05, H02, V=[7]),
        8  => ACSetTransformation(H02_05, H05, V=[1]),
        9  => ACSetTransformation(H05_06, H05, V=[7]),
        10 => ACSetTransformation(H05_06, H06, V=[1]),
        11 => ACSetTransformation(H02_07, H02, V=[7]),
        12 => ACSetTransformation(H02_07, H07, V=[1]),
        13 => ACSetTransformation(H07_08, H07, V=[7]),
        14 => ACSetTransformation(H07_08, H08, V=[1])),
    ∫(G))

decomp = StrDecomp(G, Γ)

for i in colors
    SUITE["graph coloring fixed"]["$i coloring"]["StructuredDecompositions"]["8 bags"] = @benchmarkable decide_sheaf_tree_shape($(skeletal_coloring(i)), $decomp)
end

# 12 bag case 
# 4-5 nodes per bag, 51 total nodes, 11 total adhesion nodes

# bag 1
H01 = @acset Graph begin
    V = 4
    E = 5
    src = [1, 1, 2, 2, 3]
    tgt = [2, 3, 3, 4, 4]
end

# adhesion 1, 2
H01_02 = @acset Graph begin
    V = 1
end

# bag 2
H02 = @acset Graph begin
    V = 4
    E = 5
    src = [1, 1, 2, 2, 3]
    tgt = [2, 3, 3, 4, 4]
end

# adhesion 3, 4
H02_03 = @acset Graph begin
    V = 1
end

# bag 3
H03 = @acset Graph begin
    V = 4
    E = 5
    src = [1, 1, 2, 2, 3]
    tgt = [2, 3, 3, 4, 4]
end

# adhesion 3, 4
H03_04 = @acset Graph begin
    V = 1
end

# bag 4
H04 = @acset Graph begin
    V = 4
    E = 5
    src = [1, 1, 2, 2, 3]
    tgt = [2, 3, 3, 4, 4]
end

# adhesion 4, 5
H04_05 = @acset Graph begin
    V = 1
end

# bag 5
H05 = @acset Graph begin
    V = 4
    E = 5
    src = [1, 1, 2, 2, 3]
    tgt = [2, 3, 3, 4, 4]
end

# adhesion 5, 6
H05_06 = @acset Graph begin
    V = 1
end

# bag 6
H06 = @acset Graph begin
    V = 5
    E = 6
    src = [1, 1, 2, 2, 3, 4]
    tgt = [2, 3, 3, 4, 4, 5]
end

# adhesion 3, 7
H03_07 = @acset Graph begin
    V = 1
end

# bag 7
H07 = @acset Graph begin
    V = 4
    E = 5
    src = [1, 1, 2, 2, 3]
    tgt = [2, 3, 3, 4, 4]
end

# adhesion 7, 8
H07_08 = @acset Graph begin
    V = 1
end

# bag 8
H08 = @acset Graph begin
    V = 4
    E = 5
    src = [1, 1, 2, 2, 3]
    tgt = [2, 3, 3, 4, 4]
end

# adhesion 8, 9
H08_09 = @acset Graph begin
    V = 1
end

# bag 9
H09 = @acset Graph begin
    V = 5
    E = 6
    src = [1, 1, 2, 2, 3, 4]
    tgt = [2, 3, 3, 4, 4, 5]
end

# adhesion 3, 10
H03_10   = @acset Graph begin
    V = 1
end

# bag 10
H10 = @acset Graph begin
    V = 4
    E = 5
    src = [1, 1, 2, 2, 3]
    tgt = [2, 3, 3, 4, 4]
end

# adhesion 10, 11
H10_11   = @acset Graph begin
    V = 1
end

# bag 11
H11 = @acset Graph begin
    V = 4
    E = 5
    src = [1, 1, 2, 2, 3]
    tgt = [2, 3, 3, 4, 4]
end

# adhesion 11, 12
H11_12   = @acset Graph begin
    V = 1
end

# bag 12
H12 = @acset Graph begin
    V = 5
    E = 6
    src = [1, 1, 2, 2, 3, 4]
    tgt = [2, 3, 3, 4, 4, 5]
end

G = @acset Graph begin
    V = 12
    E = 11
    src = [1, 2, 3, 3, 3, 4, 5, 7, 8, 10, 11]
    tgt = [2, 3, 4, 7, 10, 5, 6, 8, 9, 11, 12]
end

# transformations
Γ = FinDomFunctor(
    Dict(1 => H01, 2 => H02, 3 => H03, 4 => H04, 5 => H05, 6 => H06, 7 => H07, 8 => H08, 9 => H09, 10 => H10, 11 => H11, 12 => H12, 
         13 => H01_02, 14 => H02_03, 15 => H03_04, 16 => H04_05, 17 => H05_06, 18 => H03_07, 19 => H07_08, 20 => H08_09, 21 => H03_10, 22 => H10_11, 23 => H11_12),
    Dict(
        1  => ACSetTransformation(H01_02, H01, V=[4]),
        2  => ACSetTransformation(H01_02, H02, V=[1]),
        3  => ACSetTransformation(H02_03, H02, V=[4]),
        4  => ACSetTransformation(H02_03, H03, V=[1]),
        5  => ACSetTransformation(H03_04, H03, V=[4]),
        6  => ACSetTransformation(H03_04, H04, V=[1]),
        7  => ACSetTransformation(H04_05, H04, V=[4]),
        8  => ACSetTransformation(H04_05, H05, V=[1]),
        9  => ACSetTransformation(H05_06, H05, V=[4]),
        10 => ACSetTransformation(H05_06, H06, V=[1]),
        11 => ACSetTransformation(H03_07, H03, V=[4]),
        12 => ACSetTransformation(H03_07, H07, V=[1]),
        13 => ACSetTransformation(H07_08, H07, V=[4]),
        14 => ACSetTransformation(H07_08, H08, V=[1]),
        15 => ACSetTransformation(H08_09, H08, V=[4]),
        16 => ACSetTransformation(H08_09, H09, V=[1]),
        17 => ACSetTransformation(H03_10, H03, V=[4]),
        18 => ACSetTransformation(H03_10, H10, V=[1]),
        19 => ACSetTransformation(H10_11, H10, V=[4]),
        20 => ACSetTransformation(H10_11, H11, V=[1]),
        21 => ACSetTransformation(H11_12, H11, V=[4]),
        22 => ACSetTransformation(H11_12, H12, V=[1])),
    ∫(G))

decomp = StrDecomp(G, Γ)

for i in colors
    SUITE["graph coloring fixed"]["$i coloring"]["StructuredDecompositions"]["12 bags"] = @benchmarkable decide_sheaf_tree_shape($(skeletal_coloring(i)), $decomp)
end

