using LinearAlgebra
using SparseArrays
using StructuredDecompositions
using StructuredDecompositions.JunctionTrees: LinkedList
using Test


@testset "linked lists" begin
    @testset "singly linked list" begin
        list = SinglyLinkedList(zeros(Int), Vector{Int}(undef, 3))
        @test isempty(list)
        @test collect(list) == []

        @test pushfirst!(list, 1) === list
        @test !isempty(list)
        @test collect(list) == [1]

        @test pushfirst!(list, 2) === list
        @test !isempty(list)
        @test collect(list) == [2, 1]

        @test pushfirst!(list, 3) === list
        @test !isempty(list)
        @test collect(list) == [3, 2, 1]

        @test popfirst!(list) == 3
        @test !isempty(list)
        @test collect(list) == [2, 1]

        @test popfirst!(list) == 2
        @test !isempty(list)
        @test collect(list) == [1]

        @test popfirst!(list) == 1        
        @test isempty(list)
        @test collect(list) == []
    end

    @testset "doubly linked list" begin
        list = LinkedList(zeros(Int), Vector{Int}(undef, 3), Vector{Int}(undef, 3))
        @test isempty(list)
        @test collect(list) == []

        @test pushfirst!(list, 1) === list
        @test !isempty(list)
        @test collect(list) == [1]

        @test pushfirst!(list, 2) === list
        @test !isempty(list)
        @test collect(list) == [2, 1]

        @test pushfirst!(list, 3) === list
        @test !isempty(list)
        @test collect(list) == [3, 2, 1]

        @test delete!(list, 2) === list
        @test !isempty(list)
        @test collect(list) == [3, 1]

        @test delete!(list, 3) === list
        @test !isempty(list)
        @test collect(list) == [1]

        @test popfirst!(list) == 1       
        @test isempty(list)
        @test collect(list) == []
    end
end


@testset "extension" begin
    matrix = zeros(0, 0)
    @test_throws "Laplacians"      permutation(matrix; alg=Spectral())
    @test_throws "Metis"           permutation(matrix; alg=NodeND())
    @test_throws "TreeWidthSolver" permutation(matrix; alg=BT())
end


import Laplacians
import Metis
import TreeWidthSolver


@testset "representation" begin
    @test isa(repr("text/plain", MCS()), String)
    @test isa(repr("text/plain", RCM()), String)
    @test isa(repr("text/plain", AMD()), String)
    @test isa(repr("text/plain", SymAMD()), String)
    @test isa(repr("text/plain", MMD()), String)
    @test isa(repr("text/plain", NodeND()), String)
    @test isa(repr("text/plain", BT()), String)

    list = SinglyLinkedList(ones(Int), [2, 0])
    @test isa(repr("text/plain", list), String)
    list = LinkedList(ones(Int), [2, 0], [0, 1])
    @test isa(repr("text/plain", list), String)

    matrix = ones(2, 2)
    label, tree = eliminationtree(matrix)
    @test isa(repr("text/plain", tree), String)
    label, tree = supernodetree(matrix)
    @test isa(repr("text/plain", tree), String)
    label, tree = junctiontree(matrix)
    @test isa(repr("text/plain", tree), String)
end


@testset "null graph" begin
    matrix = zeros(0, 0)
    @test ischordal(matrix)
    @test isfilled(matrix)
    @test iszero(treewidth(matrix))

    @test permutation(matrix; alg=MCS())      == ([], [])
    @test permutation(matrix; alg=RCM())      == ([], [])
    @test permutation(matrix; alg=AMD())      == ([], [])
    @test permutation(matrix; alg=SymAMD())   == ([], [])
    @test permutation(matrix; alg=MMD())      == ([], []) skip=true
    @test permutation(matrix; alg=NodeND())   == ([], []) skip=true
    @test permutation(matrix; alg=Spectral()) == ([], []) skip=true
    @test permutation(matrix; alg=BT())       == ([], [])

    label, tree = junctiontree(matrix; snd=Nodal())
    @test iszero(length(tree))
    @test isnothing(rootindex(tree))
    @test iszero(treewidth(tree))
    @test iszero(nnz(tree))

    label, tree = junctiontree(matrix; snd=Maximal())
    @test iszero(length(tree))
    @test isnothing(rootindex(tree))
    @test iszero(treewidth(tree))
    @test iszero(nnz(tree))

    label, tree = junctiontree(matrix; snd=Fundamental())
    @test iszero(length(tree))
    @test isnothing(rootindex(tree))
    @test iszero(treewidth(tree))
    @test iszero(nnz(tree))
end


@testset "singleton graph" begin
    matrix = zeros(1, 1)
    @test ischordal(matrix)
    @test isfilled(matrix)
    @test iszero(treewidth(matrix))

    @test permutation(matrix; alg=MCS())      == ([1], [1])
    @test permutation(matrix; alg=RCM())      == ([1], [1])
    @test permutation(matrix; alg=AMD())      == ([1], [1])
    @test permutation(matrix; alg=SymAMD())   == ([1], [1])
    @test permutation(matrix; alg=MMD())      == ([1], [1])
    @test permutation(matrix; alg=NodeND())   == ([1], [1])
    @test permutation(matrix; alg=Spectral()) == ([1], [1]) skip=true
    @test permutation(matrix; alg=BT())       == ([1], [1])

    label, tree = junctiontree(matrix; snd=Nodal())
    @test isone(length(tree))
    @test isone(rootindex(tree))
    @test iszero(treewidth(tree))
    @test iszero(nnz(tree))
    @test isnothing(parentindex(tree, 1))
    @test isempty(childindices(tree, 1))
    @test isempty(separator(tree, 1))
    @test isempty(relative(tree, 1))
    @test isone(only(residual(tree, 1)))
    @test isone(only(tree[1]))

    label, tree = junctiontree(matrix; snd=Maximal())
    @test isone(length(tree))
    @test isone(length(tree))
    @test isone(rootindex(tree))
    @test iszero(treewidth(tree))
    @test iszero(nnz(tree))
    @test isnothing(parentindex(tree, 1))
    @test isempty(childindices(tree, 1))
    @test isempty(separator(tree, 1))
    @test isempty(relative(tree, 1))
    @test isone(only(residual(tree, 1)))
    @test isone(only(tree[1]))

    label, tree = junctiontree(matrix; snd=Fundamental())
    @test isone(length(tree))
    @test isone(rootindex(tree))
    @test iszero(treewidth(tree))
    @test iszero(nnz(tree))
    @test isnothing(parentindex(tree, 1))
    @test isempty(childindices(tree, 1))
    @test isempty(separator(tree, 1))
    @test isempty(relative(tree, 1))
    @test isone(only(residual(tree, 1)))
    @test isone(only(tree[1]))
end


@testset "vandenberghe and andersen" begin
    # Chordal Graphs and Semidefinite Optimization
    # Vandenberghe and Andersen
    matrix = [
        0  0  1  1  1  0  0  0  0  0  0  0  0  0  1  0  0
        0  0  1  1  0  0  0  0  0  0  0  0  0  0  0  0  0
        1  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
        1  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
        1  0  0  0  0  0  0  0  1  0  0  0  0  0  0  1  0
        0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  1  0
        0  0  0  0  0  0  0  1  1  0  0  0  0  0  1  0  0
        0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0
        0  0  0  0  1  1  1  0  0  0  0  0  0  0  0  0  0
        0  0  0  0  0  0  0  0  0  0  1  0  1  1  0  0  1
        0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0
        0  0  0  0  0  0  0  0  0  0  0  0  1  1  0  1  1
        0  0  0  0  0  0  0  0  0  1  0  1  0  0  0  0  0
        0  0  0  0  0  0  0  0  0  1  0  1  0  0  0  0  0
        1  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  1
        0  0  0  0  1  1  0  0  0  0  0  1  0  0  0  0  0
        0  0  0  0  0  0  0  0  0  1  0  1  0  0  1  0  0
    ]

    # Figure 4.2
    extension = [
        0  0  1  1  1  0  0  0  0  0  0  0  0  0  1  0  0
        0  0  1  1  0  0  0  0  0  0  0  0  0  0  0  0  0
        1  1  0  1  1  0  0  0  0  0  0  0  0  0  1  0  0
        1  1  1  0  1  0  0  0  0  0  0  0  0  0  1  0  0
        1  0  1  1  0  0  0  0  1  0  0  0  0  0  1  1  0
        0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  1  0
        0  0  0  0  0  0  0  1  1  0  0  0  0  0  1  0  0
        0  0  0  0  0  0  1  0  1  0  0  0  0  0  1  0  0
        0  0  0  0  1  1  1  1  0  0  0  0  0  0  1  1  0
        0  0  0  0  0  0  0  0  0  0  1  0  1  1  0  0  1
        0  0  0  0  0  0  0  0  0  1  0  0  1  1  0  0  1
        0  0  0  0  0  0  0  0  0  0  0  0  1  1  0  1  1
        0  0  0  0  0  0  0  0  0  1  1  1  0  1  0  1  1
        0  0  0  0  0  0  0  0  0  1  1  1  1  0  0  1  1
        1  0  1  1  1  0  1  1  1  0  0  0  0  0  0  1  1
        0  0  0  0  1  1  0  0  1  0  0  1  1  1  1  0  1
        0  0  0  0  0  0  0  0  0  1  1  1  1  1  1  1  0
    ]

    @test !ischordal(matrix)
    @test ischordal(extension)
    @test isfilled(tril(extension))
    @test treewidth(matrix; alg=1:17) == 4
    @test treewidth(extension; alg=1:17) == 4

    order, index = permutation(matrix; alg=MCS())
    @test isa(order, Vector{Int})
    @test isa(index, Vector{Int})
    @test length(order) == 17
    @test order[index] == 1:17

    order, index = permutation(matrix; alg=RCM())
    @test isa(order, Vector{Int})
    @test isa(index, Vector{Int})
    @test length(order) == 17
    @test order[index] == 1:17

    order, index = permutation(matrix; alg=AMD())
    @test isa(order, Vector{Int})
    @test isa(index, Vector{Int})
    @test length(order) == 17
    @test order[index] == 1:17

    order, index = permutation(matrix; alg=SymAMD())
    @test isa(order, Vector{Int})
    @test isa(index, Vector{Int})
    @test length(order) == 17
    @test order[index] == 1:17

    order, index = permutation(matrix; alg=MMD())
    @test isa(order, Vector{Int})
    @test isa(index, Vector{Int})
    @test length(order) == 17
    @test order[index] == 1:17

    order, index = permutation(matrix; alg=NodeND())
    @test isa(order, Vector{Int})
    @test isa(index, Vector{Int})
    @test length(order) == 17
    @test order[index] == 1:17

    order, index = permutation(matrix; alg=Spectral())
    @test isa(order, Vector{Int})
    @test isa(index, Vector{Int})
    @test length(order) == 17
    @test order[index] == 1:17

    order, index = permutation(matrix; alg=BT())
    @test isa(order, Vector{Int})
    @test isa(index, Vector{Int})
    @test length(order) == 17
    @test order[index] == 1:17

    # Figure 4.3
    label, tree = junctiontree(matrix; alg=1:17, snd=Nodal())
    @test length(tree) == 17
    @test rootindex(tree) == 17
    @test treewidth(tree) == 4
    @test nnz(tree) == sum(extension) ÷ 2
    @test filledgraph(tree) == tril(extension[label, label])

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

    @test map(i -> view(label, residual(tree, i)), 1:17) == [
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

    @test map(i -> view(label, separator(tree, i)), 1:17) == [
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

    @test all(1:17) do i
        tree[i] == [residual(tree, i); separator(tree, i)]
    end

    @test all(1:16) do i
        j = parentindex(tree, i)
        tree[j][relative(tree, i)] == separator(tree, i)
    end

    @test relative(tree, 17) == []


    # Figure 4.7 (left)
    label, tree = junctiontree(matrix; alg=1:17, snd=Maximal())
    @test length(tree) == 8
    @test rootindex(tree) == 8
    @test treewidth(tree) == 4
    @test nnz(tree) == sum(extension) ÷ 2
    @test filledgraph(tree) == tril(extension[label, label])

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

    @test map(i -> view(label, residual(tree, i)), 1:8) == [
        [10, 11],             # j k
        [2],                  # b
        [1, 3, 4],            # a c d
        [6],                  # f
        [7, 8],               # g h
        [5, 9],               # e i
        [15],                 # o
        [12, 13, 14, 16, 17], # l m n p q
    ]

    @test map(i -> view(label, separator(tree, i)), 1:8) == [
        [13, 14, 17], # m n q
        [3, 4],       # c d
        [5, 15],      # e o
        [9, 16],      # i p
        [9, 15],      # i o
        [15, 16],     # o p
        [16, 17],     # p q
        [],           #
    ]

    @test all(1:8) do i
        tree[i] == [residual(tree, i); separator(tree, i)]
    end

    @test all(1:7) do i
        j = parentindex(tree, i)
        tree[j][relative(tree, i)] == separator(tree, i)
    end

    @test relative(tree, 8) == []


    # Figure 4.9
    label, tree = junctiontree(matrix; alg=1:17, snd=Fundamental())
    @test length(tree) == 12
    @test rootindex(tree) == 12
    @test treewidth(tree) == 4
    @test nnz(tree) == sum(extension) ÷ 2
    @test filledgraph(tree) == tril(extension[label, label])

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

    @test map(i -> view(label, residual(tree, i)), 1:12) == [
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

    @test map(i -> view(label, separator(tree, i)), 1:12) == [
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

    @test all(1:12) do i
        tree[i] == [residual(tree, i); separator(tree, i)]
    end

    @test all(1:11) do i
        j = parentindex(tree, i)
        tree[j][relative(tree, i)] == separator(tree, i)
    end

    @test relative(tree, 12) == []
end
