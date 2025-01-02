# A junction tree.
# This type implements the abstract graph and abstract tree interfaces.
struct JunctionTree <: AbstractTree
    tree::Tree
    sndptr::Vector{Int}
    sepptr::Vector{Int}
    sepval::Vector{Int}
    relval::Vector{Int}
end


function JunctionTree()
    JunctionTree(Tree(), Int[1], Int[1], Int[], Int[])
end


# Construct a junction tree. 
function jtree(graph, ealg::EliminationAlgorithm=DEFAULT_ELIMINATION_ALGORITHM, stype::SupernodeType=DEFAULT_SUPERNODE_TYPE)
    graph = adjacency_matrix(graph)
    order = Permutation(graph, ealg)
    order, jtree!(order, permute!(OrderedGraph(graph), order), stype)
end


# Construct a junction tree. 
function jtree(graph, order::Permutation, stype::SupernodeType=DEFAULT_SUPERNODE_TYPE)
    jtree(adjacency_matrix(graph), order, stype)
end


# Construct a junction tree. 
function jtree(graph::SparseMatrixCSC, order::Permutation, stype::SupernodeType)
    order = deepcopy(order)
    order, jtree!(order, permute!(OrderedGraph(graph), order), stype)
end


# Construct a junction tree. 
function jtree!(order::Permutation, graph::OrderedGraph, stype::SupernodeType)
    tree, sndptr, sepptr = stree!(order, graph, stype)
    sepval = sepvals(graph, tree, sndptr, sepptr)
    relval = relvals(tree, sndptr, sepptr, sepval)
    JunctionTree(tree, sndptr, sepptr, sepval, relval)
end


# Get the seperators of every node of a supernodal elimination tree.
function sepvals(graph::OrderedGraph, stree::Tree, sndptr::AbstractVector, sepptr::AbstractVector)
    temporary = zeros(Int, nv(graph))
    sepval = Vector{Int}(undef, last(sepptr) - 1)
    p = 1

    for j in vertices(stree)
        for v in inneighbors(graph, sndptr[j])
            if sndptr[j + 1] <= v
                sepval[p] = v
                temporary[v] = j
                p += 1
            end
        end

        for i in childindices(stree, j), v in view(sepval, sepptr[i]:sepptr[i + 1] - 1)
            if sndptr[j + 1] <= v && temporary[v] != j
                sepval[p] = v
                temporary[v] = j
                p += 1
            end
        end

        sort!(view(sepval, sepptr[j]:sepptr[j + 1] - 1))
    end

    sepval
end


# Get the relative indices of every node of a supernodal elimination tree.
function relvals(stree::Tree, sndptr::AbstractVector, sepptr::AbstractVector, sepval::AbstractVector)
    relval = Vector{Int}(undef, length(sepval))
    p = 1

    for edge in Graphs.edges(stree)
        i = Graphs.dst(edge)
        j = Graphs.src(edge)

        while p < sepptr[i + 1] && sepval[p] < sndptr[j + 1]
            relval[p] = sepval[p] - sndptr[j] + 1
            p += 1
        end

        q = sepptr[j]

        while p < sepptr[i + 1]
            q += searchsortedfirst(view(sepval, q:sepptr[j + 1] - 1), sepval[p])
            relval[p] = q - sepptr[j] - sndptr[j] + sndptr[j + 1]
            p += 1
        end
    end

    relval
end


# Get the clique at node i.
function clique(tree::JunctionTree, i::Integer)
    SumVector(residual(tree, i), seperator(tree, i))
end


# Compute the width of a junction tree
function treewidth(tree::JunctionTree)
    maximum(vertices(tree)) do i
        length(residual(tree, i)) + length(seperator(tree, i)) - 1
    end
end


# Get the residual at node i.
function residual(jtree::JunctionTree, i::Integer)
    jtree.sndptr[i]:jtree.sndptr[i + 1] - 1
end


# Get the seperator at node i.
function seperator(jtree::JunctionTree, i::Integer)
    view(jtree.sepval, jtree.sepptr[i]:jtree.sepptr[i + 1] - 1)
end


# Get the relative indices at node i.
function relative(jtree::JunctionTree, i::Integer)
    view(jtree.relval, jtree.sepptr[i]:jtree.sepptr[i + 1] - 1)
end


##########################
# Indexed Tree Interface #
##########################


function firstchildindex(tree::JunctionTree, i::Integer)
    firstchildindex(tree.tree, i)
end


function AbstractTrees.rootindex(tree::JunctionTree)
    rootindex(tree.tree)
end


function AbstractTrees.parentindex(tree::JunctionTree, i::Integer)
    parentindex(tree.tree, i)
end


function AbstractTrees.nextsiblingindex(tree::JunctionTree, i::Integer)
    nextsiblingindex(tree.tree, i)
end


function AbstractTrees.nodevalue(tree::JunctionTree, i::Integer)
    clique(tree, i)
end


############################
# Abstract Graph Interface #
############################


function Graphs.nv(tree::JunctionTree)
    nv(tree.tree)
end
