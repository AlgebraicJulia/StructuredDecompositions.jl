# This file benchmarks the StructuredDecompositions method of graph coloring

# RESULTS = sheaf tree skeletal coloring is n linear for false returns, worse than n exponential for true returns
    # shouldn't be worse than exponential for true returns so I will take a look at it after considering k

# The following structures and functions were pulled from DecidingSheaves.jl 


function K(n::Integer)
    complete_graph(Graph, n)
end


struct Coloring
    n::Int    
end


function (coloring::Coloring)(graph::Graph)
    FinSet(homomorphisms(graph, K(coloring.n)))
end


function (coloring::Coloring)(f::ACSetTransformation)  
    FinFunction(λ -> compose(f, λ), coloring(codom(f)), coloring(dom(f)))
end
 
 
function skeletal_coloring(n::Integer)
    skeleton ∘ Coloring(n)
end 


# Decomposition 1
# small n per bag(4) and small bags k(2), 3 coloring

# bag 1
H1 = @acset Graph begin
    V = 4
    E = 4
    src = [1, 2, 3, 4]
    tgt = [2, 3, 4, 1]
end

# adhesion 1, 2
H12 = @acset Graph begin
    V = 2
end

# bag 2
H2 = @acset Graph begin
    V = 4
    E = 5
    src = [1, 2, 3, 4, 1]
    tgt = [2, 3, 4, 1, 3]
end

# decomposition shape
G = @acset Graph begin
    V = 2
    E = 1
    src = [1]
    tgt = [2]
end

# transformations
Γ = FinDomFunctor(
    Dict(1 => H1, 2 => H2, 3 => H12),
    Dict(
        1 => ACSetTransformation(H12, H1, V=[3, 4]),
        2 => ACSetTransformation(H12, H2, V=[2, 1])),
    ∫(G))

decomp = StrDecomp(G, Γ)
graph = ob(colimit(decomp))

for i in (2, 3, 4)
    SUITE["GraphColoring"]["Decomposition 1"]["$i Coloring"]["HomSearch"] = @benchmarkable is_homomorphic($graph, $(K(i)))
    SUITE["GraphColoring"]["Decomposition 1"]["$i Coloring"]["StructuredDecompositions"] = @benchmarkable decide_sheaf_tree_shape($(skeletal_coloring(i)), $decomp)
end

# Decomposition 2
# medium n per bag(10) and small bags k(2), 3 coloring

# bag 1
H1 = @acset Graph begin
    V = 10
    E = 16
    src = [1, 1, 1, 2, 2, 3, 4, 4, 4, 6, 6, 7, 7, 7, 8, 9]
    tgt = [2, 5, 6, 3, 8, 4, 5, 9, 10, 7, 10, 8, 9, 10, 9, 10]
end

# adhesion 1, 2
H12 = @acset Graph begin
    V = 2
end

# bag 2
H2 = @acset Graph begin
    V = 10
    E = 18
    src = [1, 1, 2, 2, 2, 2, 3, 4, 4, 5, 5, 6, 6, 7, 7, 7, 8, 9]
    tgt = [2, 9, 3, 4, 9, 10, 4, 5, 10, 6, 10, 7, 10, 8, 9, 10, 9, 10]
end

# decomposition shape
G = @acset Graph begin
    V = 2
    E = 1
    src = [1]
    tgt = [2]
end

# transformations
Γ = FinDomFunctor(
    Dict(1 => H1, 2 => H2, 3 => H12),
    Dict(
        1 => ACSetTransformation(H12, H1, V=[4, 2]),
        2 => ACSetTransformation(H12, H2, V=[1, 3])),
    ∫(G))

decomp = StrDecomp(G, Γ)
graph = ob(colimit(decomp))

for i in (2, 3, 4)
    SUITE["GraphColoring"]["Decomposition 2"]["$i Coloring"]["HomSearch"] = @benchmarkable is_homomorphic($graph, $(K(i)))
    SUITE["GraphColoring"]["Decomposition 2"]["$i Coloring"]["StructuredDecompositions"] = @benchmarkable decide_sheaf_tree_shape($(skeletal_coloring(i)), $decomp)
end


#=
# Decomposition 3 
# large n per bag(20) and small bags k(2), 3 coloring

# bag 1
H1 = @acset Graph begin
    V = 20
    E = 34
    src = [1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13, 14, 15, 15, 16, 16, 17, 17, 18]
    tgt = [2, 3, 4, 5, 4, 6, 5, 7, 8, 9, 7, 10, 8, 11, 9, 12, 13, 14, 11, 15, 12, 15, 13, 16, 14, 17, 18, 16, 19, 17, 19, 18, 19, 20]
end
  
# adhesion 1, 2
H12 = @acset Graph begin
    V = 2
end
  
# bag 2
H2 = @acset Graph begin
    V = 20
    E = 34
    src = [1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13, 14, 15, 15, 16, 16, 17, 17, 18]
    tgt = [2, 3, 4, 5, 4, 6, 5, 7, 8, 9, 7, 10, 8, 11, 9, 12, 13, 14, 11, 15, 12, 15, 13, 16, 14, 17, 18, 16, 19, 17, 19, 18, 19, 20]
end
  
G = @acset Graph begin
    V = 2
    E = 1
    src = [1]
    tgt = [2]
end
  
# transformations
Γ = FinDomFunctor(
    Dict(1 => H1, 2 => H2, 3 => H12),
    Dict(
        1 => ACSetTransformation(H12, H1, V=[1, 2]),
        2 => ACSetTransformation(H12, H2, V=[1, 2])),
  ∫(G))

decomp = StrDecomp(G, Γ)
graph = ob(colimit(decomp))

for i in (2, 3, 4)
    SUITE["GraphColoring"]["Decomposition 3"]["$i Coloring"]["HomSearch"] = @benchmarkable is_homomorphic($graph, $(K(i)))
    SUITE["GraphColoring"]["Decomposition 3"]["$i Coloring"]["StructuredDecompositions"] = @benchmarkable decide_sheaf_tree_shape($(skeletal_coloring(i)), $decomp)
end
=#


# Decomposition 4 
# small n per bag(5) with medium bags k(4), 3 coloring

# bag 1
H1 = @acset Graph begin
    V = 5
    E = 7
    src = [1, 2, 2, 2, 3, 4, 5]
    tgt = [2, 3, 4, 5, 4, 5, 1]
end

# adhesion 1, 2
H12 = @acset Graph begin
    V = 1
end

# bag 2
H2 = @acset Graph begin
    V = 5
    E = 7
    src = [1, 2, 2, 2, 3, 4, 5]
    tgt = [2, 3, 4, 5, 4, 5, 1]
end

# adhesion 3, 4
H23 = @acset Graph begin
    V = 1
end

# bag 3
H3 = @acset Graph begin
    V = 5
    E = 7
    src = [1, 2, 2, 2, 3, 4, 5]
    tgt = [2, 3, 4, 5, 4, 5, 1]
end

# adhesion 3, 4
H34 = @acset Graph begin
    V = 1
end

# bag 4
H4 = @acset Graph begin
    V = 5
    E = 7
    src = [1, 2, 2, 2, 3, 4, 5]
    tgt = [2, 3, 4, 5, 4, 5, 1]
end

# decomposition shape
G = @acset Graph begin
    V = 4
    E = 3
    src = [1, 2, 3]
    tgt = [2, 3, 4]
end

# transformations
Γ = FinDomFunctor(
    Dict(1 => H1, 2 => H2, 3 => H3, 4 => H4, 5 => H12, 6 => H23, 7 => H34),
    Dict(
        1 => ACSetTransformation(H12, H1, V=[1]),
        2 => ACSetTransformation(H12, H2, V=[1]),
        3 => ACSetTransformation(H23, H2, V=[1]),
        4 => ACSetTransformation(H23, H3, V=[1]),
        5 => ACSetTransformation(H34, H3, V=[1]),
        6 => ACSetTransformation(H34, H4, V=[1])),
    ∫(G))

decomp = StrDecomp(G, Γ)
graph = ob(colimit(decomp))

for i in (2, 3, 4)
    SUITE["GraphColoring"]["Decomposition 4"]["$i Coloring"]["HomSearch"] = @benchmarkable is_homomorphic($graph, $(K(i)))
    SUITE["GraphColoring"]["Decomposition 4"]["$i Coloring"]["StructuredDecompositions"] = @benchmarkable decide_sheaf_tree_shape($(skeletal_coloring(i)), $decomp)
end


# Decomposition 5 
# medium n per bag(10) with medium bags k(4), 3 coloring

# bag 1
H1 = @acset Graph begin
    V = 10
    E = 16
    src = [1, 1, 1, 2, 2, 3, 4, 4, 4, 6, 6, 7, 7, 7, 8, 9]
    tgt = [2, 5, 6, 3, 8, 4, 5, 9, 10, 7, 10, 8, 9, 10, 9, 10]
end

# adhesion 1, 2
H12 = @acset Graph begin
    V = 1
end

# bag 2
H2 = @acset Graph begin
    V = 10
    E = 18
    src = [1, 1, 2, 2, 2, 2, 3, 4, 4, 5, 5, 6, 6, 7, 7, 7, 8, 9]
    tgt = [2, 9, 3, 4, 9, 10, 4, 5, 10, 6, 10, 7, 10, 8, 9, 10, 9, 10]
end

# adhesion 3, 4
H23 = @acset Graph begin
    V = 1
end

# bag 3
H3 = @acset Graph begin
    V = 10
    E = 16
    src = [1, 1, 1, 2, 2, 3, 4, 4, 4, 6, 6, 7, 7, 7, 8, 9]
    tgt = [2, 5, 6, 3, 8, 4, 5, 9, 10, 7, 10, 8, 9, 10, 9, 10]
end

# adhesion 3, 4
H34 = @acset Graph begin
    V = 1
end

# bag 4
H4 = @acset Graph begin
    V = 10
    E = 18
    src = [1, 1, 2, 2, 2, 2, 3, 4, 4, 5, 5, 6, 6, 7, 7, 7, 8, 9]
    tgt = [2, 9, 3, 4, 9, 10, 4, 5, 10, 6, 10, 7, 10, 8, 9, 10, 9, 10]
end

G = @acset Graph begin
    V = 4
    E = 3
    src = [1, 3, 4]
    tgt = [2, 2, 3]
end

# transformations
Γ = FinDomFunctor(
    Dict(1 => H1, 2 => H2, 3 => H3, 4 => H4, 5 => H12, 6 => H23, 7 => H34),
    Dict(
        1 => ACSetTransformation(H12, H1, V=[4]),
        2 => ACSetTransformation(H12, H2, V=[1]),
        3 => ACSetTransformation(H23, H2, V=[1]),
        4 => ACSetTransformation(H23, H3, V=[1]),
        5 => ACSetTransformation(H34, H3, V=[1]),
        6 => ACSetTransformation(H34, H4, V=[4])),
    ∫(G))

decomp = StrDecomp(G, Γ)
graph = ob(colimit(decomp))

for i in (2, 3, 4)
    SUITE["GraphColoring"]["Decomposition 5"]["$i Coloring"]["HomSearch"] = @benchmarkable is_homomorphic($graph, $(K(i)))
    SUITE["GraphColoring"]["Decomposition 5"]["$i Coloring"]["StructuredDecompositions"] = @benchmarkable decide_sheaf_tree_shape($(skeletal_coloring(i)), $decomp)
end

# Decomposition 6 
# large n per bag(20) with medium bags k(4), 3 coloring

# bag 1
H1 = @acset Graph begin
    V = 20
    E = 34
    src = [1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13, 14, 15, 15, 16, 16, 17, 17, 18]
    tgt = [2, 3, 4, 5, 4, 6, 5, 7, 8, 9, 7, 10, 8, 11, 9, 12, 13, 14, 11, 15, 12, 15, 13, 16, 14, 17, 18, 16, 19, 17, 19, 18, 19, 20]
end

# adhesion 1, 2
H12 = @acset Graph begin
    V = 1
end

#   bag 2
H2 = @acset Graph begin
    V = 20
    E = 34
    src = [1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13, 14, 15, 15, 16, 16, 17, 17, 18]
    tgt = [2, 3, 4, 5, 4, 6, 5, 7, 8, 9, 7, 10, 8, 11, 9, 12, 13, 14, 11, 15, 12, 15, 13, 16, 14, 17, 18, 16, 19, 17, 19, 18, 19, 20]
end

# adhesion 3, 4
H23 = @acset Graph begin
    V = 1
end

#bag 3
H₃ = @acset Graph begin
    V = 20
    E = 34
    src = [1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13, 14, 15, 15, 16, 16, 17, 17, 18]
    tgt = [2, 3, 4, 5, 4, 6, 5, 7, 8, 9, 7, 10, 8, 11, 9, 12, 13, 14, 11, 15, 12, 15, 13, 16, 14, 17, 18, 16, 19, 17, 19, 18, 19, 20]
end

# adhesion 3, 4
H34 = @acset Graph begin
  V = 1
end

# bag 4
H4 = @acset Graph begin
    V = 20
    E = 34
    src = [1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13, 14, 15, 15, 16, 16, 17, 17, 18]
    tgt = [2, 3, 4, 5, 4, 6, 5, 7, 8, 9, 7, 10, 8, 11, 9, 12, 13, 14, 11, 15, 12, 15, 13, 16, 14, 17, 18, 16, 19, 17, 19, 18, 19, 20]
end

G = @acset Graph begin
    V = 4
    E = 3
    src = [1, 2, 3]
    tgt = [2, 3, 4]
end

# transformations
Γ = FinDomFunctor(
    Dict(1 => H1, 2 => H2, 3 => H3, 4 => H4, 5 => H12, 6 => H23, 7 => H34),
    Dict(
        1 => ACSetTransformation(H12, H1, V=[1]),
        2 => ACSetTransformation(H12, H2, V=[1]),
        3 => ACSetTransformation(H23, H2, V=[1]),
        4 => ACSetTransformation(H23, H3, V=[1]),
        5 => ACSetTransformation(H34, H3, V=[1]),
        6 => ACSetTransformation(H34, H4, V=[1])),
    ∫(G))

decomp = StrDecomp(G, Γ)
graph = ob(colimit(decomp))

for i in (2, 3, 4)
    SUITE["GraphColoring"]["Decomposition 6"]["$i Coloring"]["HomSearch"] = @benchmarkable is_homomorphic($graph, $(K(i)))
    SUITE["GraphColoring"]["Decomposition 6"]["$i Coloring"]["StructuredDecompositions"] = @benchmarkable decide_sheaf_tree_shape($(skeletal_coloring(i)), $decomp)
end

# Decomposition 7 
# small n per bag(5) with large bags k(10), 3 coloring

# bag 1
H1 = @acset Graph begin
    V = 5
    E = 7
    src = [1, 2, 2, 2, 3, 4, 5]
    tgt = [2, 3, 4, 5, 4, 5, 1]
end

# adhesion 1, 2
H12 = @acset Graph begin
    V = 1
end

# bag 2
H2 = @acset Graph begin
    V = 5
    E = 7
    src = [1, 2, 2, 2, 3, 4, 5]
    tgt = [2, 3, 4, 5, 4, 5, 1]
end

# adhesion 3, 4
H23 = @acset Graph begin
    V = 1
end

# bag 3
H3 = @acset Graph begin
    V = 5
    E = 7
    src = [1, 2, 2, 2, 3, 4, 5]
    tgt = [2, 3, 4, 5, 4, 5, 1]
end

# adhesion 3, 4
H34 = @acset Graph begin
    V = 1
end

# bag 4
H4 = @acset Graph begin
    V = 5
    E = 7
    src = [1, 2, 2, 2, 3, 4, 5]
    tgt = [2, 3, 4, 5, 4, 5, 1]
end

# adhesion 4, 5
H45 = @acset Graph begin
    V = 1
end

# bag 5
H5 = @acset Graph begin
    V = 5
    E = 7
    src = [1, 2, 2, 2, 3, 4, 5]
    tgt = [2, 3, 4, 5, 4, 5, 1]
end

# adhesion 5, 6
H56 = @acset Graph begin
    V = 1
end

# bag 6
H6 = @acset Graph begin
    V = 5
    E = 7
    src = [1, 2, 2, 2, 3, 4, 5]
    tgt = [2, 3, 4, 5, 4, 5, 1]
end

# adhesion 6, 7
H67 = @acset Graph begin
    V = 1
end

# bag 7
H7 = @acset Graph begin
    V = 5
    E = 7
    src = [1, 2, 2, 2, 3, 4, 5]
    tgt = [2, 3, 4, 5, 4, 5, 1]
end

# adhesion 7, 8
H78 = @acset Graph begin
    V = 1
end

# bag 8
H8 = @acset Graph begin
    V = 5
    E = 7
    src = [1, 2, 2, 2, 3, 4, 5]
    tgt = [2, 3, 4, 5, 4, 5, 1]
end

# adhesion 8, 9
H89 = @acset Graph begin
    V = 1
end

# bag 9
H9 = @acset Graph begin
    V = 5
    E = 7
    src = [1, 2, 2, 2, 3, 4, 5]
    tgt = [2, 3, 4, 5, 4, 5, 1]
end

# adhesion 9, 10
H90   = @acset Graph begin
    V = 1
end

# bag 10
H0 = @acset Graph begin
    V = 5
    E = 7
    src = [1, 2, 2, 2, 3, 4, 5]
    tgt = [2, 3, 4, 5, 4, 5, 1]
end

# decomposition shape
G = @acset Graph begin
    V = 10
    E = 9
    src = [1, 2, 3, 4, 5, 6, 7, 8, 9]
    tgt = [2, 3, 4, 5, 6, 7, 8, 9, 10]
end

# transformations
Γ = FinDomFunctor(
    Dict(1 => H1, 2 => H2, 3 => H3, 4 => H4, 5 => H5, 6 => H6, 7 => H7, 8 => H8, 9 => H9,
        10 => H0, 11 => H12, 12 => H23, 13 => H34, 14 => H45, 15 => H56, 16 => H67, 17 => H78, 18 => H89, 19 => H90),
    Dict(
        1  => ACSetTransformation(H12, H1, V=[1]),
        2  => ACSetTransformation(H12, H2, V=[1]),
        3  => ACSetTransformation(H23, H2, V=[1]),
        4  => ACSetTransformation(H23, H3, V=[1]),
        5  => ACSetTransformation(H34, H3, V=[1]),
        6  => ACSetTransformation(H34, H4, V=[1]),
        7  => ACSetTransformation(H45, H4, V=[1]),
        8  => ACSetTransformation(H45, H5, V=[1]),
        9  => ACSetTransformation(H56, H5, V=[1]),
        10 => ACSetTransformation(H56, H6, V=[1]),
        11 => ACSetTransformation(H67, H6, V=[1]),
        12 => ACSetTransformation(H67, H7, V=[1]),
        13 => ACSetTransformation(H78, H7, V=[1]),
        14 => ACSetTransformation(H78, H8, V=[1]),
        15 => ACSetTransformation(H89, H8, V=[1]),
        16 => ACSetTransformation(H89, H9, V=[1]),
        17 => ACSetTransformation(H90, H9, V=[1]),
        18 => ACSetTransformation(H90, H0, V=[1])),
    ∫(G))

decomp = StrDecomp(G, Γ)
graph = ob(colimit(decomp))

for i in (2, 3, 4)
    SUITE["GraphColoring"]["Decomposition 6"]["$i Coloring"]["HomSearch"] = @benchmarkable is_homomorphic($graph, $(K(i)))
    SUITE["GraphColoring"]["Decomposition 6"]["$i Coloring"]["StructuredDecompositions"] = @benchmarkable decide_sheaf_tree_shape($(skeletal_coloring(i)), $decomp)
end
