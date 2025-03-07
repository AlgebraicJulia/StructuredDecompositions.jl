module TestDecompositions

using Test
using PartialFunctions

using StructuredDecompositions.Decompositions 
using StructuredDecompositions.FunctorUtils

using Catlab.Graphics
using Catlab.Graphs
using Catlab.ACSetInterface
using Catlab.CategoricalAlgebra


# Define the instance #######################
# bag 1
H1 = @acset Graph begin
    V = 3
    E = 2
    src = [1, 2]
    tgt = [2, 3]
end

# adhesion 1, 2
H12 = @acset Graph begin
    V = 2
end

# bag 2
H2 = @acset Graph begin
    V = 4
    E = 3
    src = [1, 2, 3]
    tgt = [2, 3, 4]
end

# adhesion 2, 3
H23 = @acset Graph begin
    V = 1
end

# bag 3
H3 = @acset Graph begin
    V = 2
    E = 1
    src = [1]
    tgt = [2]
end

# Make the decomp ###########
# The shape
G = @acset Graph begin
    V = 3
    E = 2
    src = [1, 2]
    tgt = [2, 3]
end

# the functor
Γ = FinDomFunctor(
    Dict(1 => H1, 2 => H2, 3 => H3, 4 => H12, 5 => H23),
    Dict(
        1 => ACSetTransformation(H12, H1, V=[1, 3]),
        2 => ACSetTransformation(H23, H2, V=[1]),
        3 => ACSetTransformation(H12, H2, V=[4, 1]),
        4 => ACSetTransformation(H23, H3, V=[1])),
    ∫(G))

# the decomposition
bigdecomp = StrDecomp(G, Γ)

#f = ACSetTransformation(Γ₀[4], Γ₀[1], V=[1, 3])

@test H1 ∈ bags(bigdecomp) && H2 ∈ bags(bigdecomp) && !(H12 ∈ bags(bigdecomp))

v𝐃 = 𝐃 $ vs

bigdecomp_to_sets = v𝐃(bigdecomp)
@test all(s -> dom(s[1]) == dom(s[2]), adhesionSpans(bigdecomp_to_sets))

s𝐃 = 𝐃 $ skeleton

bigdecomp_skeleton = s𝐃(bigdecomp_to_sets)

@test bags(bigdecomp_skeleton) == map(FinSet, [3,4,2])
@test adhesions(bigdecomp_skeleton) == map(FinSet, [2,1])
@test all(s -> dom(s[1]) == dom(s[2]), adhesionSpans(bigdecomp_skeleton))


##################################
# Integration with JunctionTrees #
##################################


graph = SymmetricGraph(17)

add_edges!(graph,
    [1, 1, 1, 1,  2, 2, 5, 5,  6, 6,  7, 7, 7,  10, 10, 10, 10, 12, 12, 12, 12, 15],
    [3, 4, 5, 15, 3, 4, 9, 16, 9, 16, 8, 9, 15, 11, 13, 14, 17, 13, 14, 16, 17, 17])

decomposition = StrDecomp(graph; alg=1:17)

@test decomposition.decomp_shape == @acset Graph begin
    V = 8
    E = 7
    tgt = [1, 2, 3, 4, 5, 6, 7]
    src = [8, 3, 6, 6, 6, 7, 8]
end

@test map(i -> ob_map(decomposition.diagram, i), 1:15) == [
    induced_subgraph(graph, [10, 11, 13, 14, 17]), # j k m n q
    induced_subgraph(graph, [2, 3, 4]),            # b c d
    induced_subgraph(graph, [1, 3, 4, 5, 15]),     # a c d e o
    induced_subgraph(graph, [6, 9, 16]),           # f i p
    induced_subgraph(graph, [7, 8, 9, 15]),        # g h i o
    induced_subgraph(graph, [5, 9, 15, 16]),       # e i o p
    induced_subgraph(graph, [15, 16, 17]),         # o p q
    induced_subgraph(graph, [12, 13, 14, 16, 17]), # l m n p q
    induced_subgraph(graph, [13, 14, 17]),         # m n q
    induced_subgraph(graph, [3, 4]),               # c d
    induced_subgraph(graph, [5, 15]),              # e o
    induced_subgraph(graph, [9, 16]),              # i p
    induced_subgraph(graph, [9, 15]),              # i o
    induced_subgraph(graph, [15, 16]),             # o p
    induced_subgraph(graph, [16, 17]),             # p q
]

@test map(i -> hom_map(decomposition.diagram, i), 1:14) == [
    ACSetTransformation(induced_subgraph(graph, [13, 14, 17]), induced_subgraph(graph, [12, 13, 14, 16, 17]), V=[2, 3, 5], E=Int[]), # m n q → l m n p q
    ACSetTransformation(induced_subgraph(graph, [3, 4]),       induced_subgraph(graph, [1, 3, 4, 5, 15]),     V=[2, 3],    E=Int[]), # c d   → a c d e o
    ACSetTransformation(induced_subgraph(graph, [5, 15]),      induced_subgraph(graph, [5, 9, 15, 16]),       V=[1, 3],    E=Int[]), # e o   → e i o p
    ACSetTransformation(induced_subgraph(graph, [9, 16]),      induced_subgraph(graph, [5, 9, 15, 16]),       V=[2, 4],    E=Int[]), # i p   → e i o p
    ACSetTransformation(induced_subgraph(graph, [9, 15]),      induced_subgraph(graph, [5, 9, 15, 16]),       V=[2, 3],    E=Int[]), # i o   → e i o p
    ACSetTransformation(induced_subgraph(graph, [15, 16]),     induced_subgraph(graph, [15, 16, 17]),         V=[1, 2],    E=Int[]), # o p   → o p q
    ACSetTransformation(induced_subgraph(graph, [16, 17]),     induced_subgraph(graph, [12, 13, 14, 16, 17]), V=[4, 5],    E=Int[]), # p q   → l m n p q
    ACSetTransformation(induced_subgraph(graph, [13, 14, 17]), induced_subgraph(graph, [10, 11, 13, 14, 17]), V=[3, 4, 5], E=Int[]), # m n q → j k m n q
    ACSetTransformation(induced_subgraph(graph, [3, 4]),       induced_subgraph(graph, [2, 3, 4]),            V=[2, 3],    E=Int[]), # c d   → b c d
    ACSetTransformation(induced_subgraph(graph, [5, 15]),      induced_subgraph(graph, [1, 3, 4, 5, 15]),     V=[4, 5],    E=Int[]), # e o   → a c d e o
    ACSetTransformation(induced_subgraph(graph, [9, 16]),      induced_subgraph(graph, [6, 9, 16]),           V=[2, 3],    E=Int[]), # i p   → f i p
    ACSetTransformation(induced_subgraph(graph, [9, 15]),      induced_subgraph(graph, [7, 8, 9, 15]),        V=[3, 4],    E=Int[]), # i o   → g h i o
    ACSetTransformation(induced_subgraph(graph, [15, 16]),     induced_subgraph(graph, [5, 9, 15, 16]),       V=[3, 4],    E=Int[]), # o p   → e i o p
    ACSetTransformation(induced_subgraph(graph, [16, 17]),     induced_subgraph(graph, [15, 16, 17]),         V=[2, 3],    E=Int[]), # p q   → o p q
]


end
