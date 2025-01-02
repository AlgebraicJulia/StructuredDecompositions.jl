const AbstractGraph = Union{Graphs.AbstractGraph, Catlab.HasGraph, SparseMatrixCSC}


# Construct an order using the reverse Cuthill-McKee algorithm. Uses CuthillMcKee.jl.
function Permutation(graph::AbstractGraph, ealg::CuthillMcKeeJL_RCM)
    order = CuthillMcKee.symrcm(adjacency_matrix(graph))
    Permutation(order)
end


# Construct an order using the reverse Cuthill-McKee algorithm. Uses SYMRCMJL_RCM.jl.
function Permutation(graph::AbstractGraph, ealg::SymRCMJL_RCM)
    order = SymRCM.symrcm(adjacency_matrix(graph))
    Permutation(order)
end


# Construct an order using the approximate minimum degree algorithm. Uses AMD.jl.
function Permutation(graph::AbstractGraph, ealg::AMDJL_AMD)
    order = AMD.amd(adjacency_matrix(graph), ealg.meta)
    Permutation(order)
end


# Construct an order using the SYMAMD algorithm. Uses AMD.jl.
function Permutation(graph::AbstractGraph, ealg::AMDJL_SYMAMD)
    order = AMD.symamd(adjacency_matrix(graph), ealg.meta)
    Permutation(order)
end


# Construct an order using the nested dissection heuristic. Uses Metis.jl.
function Permutation(graph::AbstractGraph, ealg::MetisJL_ND)
    order, index = Metis.permutation(adjacency_matrix(graph))
    Permutation(order, index)
end


# Construct an order using the Bouchitte-Todinca algorithm. Uses TreeWidthSolver.jl.
function Permutation(graph::AbstractGraph, ealg::TreeWidthSolverJL_BT)
    order = reverse!(reduce(vcat, TreeWidthSolver.elimination_order(simple(graph))))
    Permutation(order)
end


# Construct an order using the maximum cardinality search algorithm.
function Permutation(graph::AbstractGraph, ealg::MCS)
    order = mcs(graph)
    Permutation(order)
end


# Construct an order using a descent-first search.
function Permutation(graph::AbstractGraph, ealg::DFS)
    order = dfs(graph)
    Permutation(order)
end


function is_directed(::Type{T}) where T <: Graphs.AbstractGraph
    Graphs.is_directed(T)
end


function is_directed(::Type{T}) where T <: Catlab.HasGraph
    Catlab.is_directed(T)
end


function is_directed(::Type{T}) where T <: SparseMatrixCSC
    false
end


function is_directed(graph::T) where T <: AbstractGraph
    is_directed(T)
end


function nv(graph::Graphs.AbstractGraph)
    Graphs.nv(graph)
end


function nv(graph::Catlab.HasGraph)
    Catlab.nv(graph)
end


function nv(graph::SparseMatrixCSC)
    size(graph, 2)
end


function ne(graph::Graphs.AbstractGraph)
    Graphs.ne(graph)
end


function ne(graph::Catlab.HasGraph, is_directed::Val{true}=Val(is_directed(graph)))
    Catlab.ne(graph)
end


function ne(graph::Catlab.HasGraph, is_directed::Val{false})
    Catlab.ne(graph) ÷ 2
end


function ne(graph::SparseMatrixCSC)
    nnz(graph) ÷ 2
end


function vertices(graph::Graphs.AbstractGraph)
    Graphs.vertices(graph)
end


function vertices(graph::Catlab.HasGraph)
    Catlab.vertices(graph)
end


function vertices(graph::SparseMatrixCSC)
    axes(graph, 2)
end


function inneighbors(graph::Graphs.AbstractGraph, v::Integer)
    Graphs.inneighbors(graph, v)
end


function inneighbors(graph::Catlab.HasGraph, v::Integer)
    view(view(graph, :src), Catlab.incident(graph, v, :tgt))
end


function inneighbors(graph::SparseMatrixCSC, v::Integer)
    outneighbors(graph, v)
end


function outneighbors(graph::Graphs.AbstractGraph, v::Integer)
    Graphs.outneighbors(graph, v)
end


function outneighbors(graph::Catlab.HasGraph, v::Integer)
    view(view(graph, :tgt), Catlab.incident(graph, v, :src))
end


function outneighbors(graph::SparseMatrixCSC, v::Integer)
    view(rowvals(graph), nzrange(graph, v))
end


function all_neighbors(graph::AbstractGraph, v::Integer, is_directed::Val{true}=Val(is_directed(graph)))
    SumVector(inneighbors(graph, v), outneighbors(graph, v))
end


function all_neighbors(graph::AbstractGraph, v::Integer, is_directed::Val{false})
    outneighbors(graph, v)
end


# Construct the adjacency matrix of a graph.
function adjacency_matrix(graph::AbstractGraph)
    rowval = Vector{Int}(undef, 2ne(graph))
    colptr = Vector{Int}(undef, nv(graph) + 1)
    colptr[1] = 1

    for i in vertices(graph)
        p = colptr[i]

        for j in all_neighbors(graph, i)
            rowval[p] = j
            p += 1
        end

        colptr[i + 1] = p
        sort!(view(rowval, colptr[i]:colptr[i + 1] - 1))
    end

    nzval = ones(Bool, 2ne(graph))
    SparseMatrixCSC(nv(graph), nv(graph), colptr, rowval, nzval)
end


# Construct the adjacency matrix of a graph.
function adjacency_matrix(graph::SparseMatrixCSC)
    graph
end


function simple(graph::AbstractGraph)
    fadjlist = Vector{Vector{Int}}(undef, nv(graph))

    for i in vertices(graph)
        fadjlist[i] = sort(all_neighbors(graph, i))
    end

    Graphs.Graph(ne(graph), fadjlist)
end


function simple(graph::Graphs.Graph)
    graph
end


# Simple Linear-Time Algorithms to Test Chordality of Graphs, Test Acyclicity of Hypergraphs, and Selectively Reduce Acyclic Hypergraphs
# Tarjan and Yannakakis
# Maximum Cardinality Search
function mcs(graph::AbstractGraph)
    n = nv(graph)
    α = Vector{Int}(undef, n)
    size = Vector{Int}(undef, n)
    set = Vector{LinkedList{Int}}(undef, n)
    pointer = Vector{ListNode{Int}}(undef, n)

    for i in 1:n
        size[i] = 1
        set[i] = LinkedList{Int}()
        pointer[i] = push!(set[1], i)
    end

    i = n
    j = 1

    while i >= 1
        v = first(set[j])
        deleteat!(set[j], pointer[v])        
        α[i] = v
        size[v] = 0

        for w in all_neighbors(graph, v)
            if size[w] >= 1
                deleteat!(set[size[w]], pointer[w])
                size[w] += 1
                pointer[w] = push!(set[size[w]], w)
            end
        end

        i -= 1
        j += 1

        while j >= 1 && isempty(set[j])
            j -= 1
        end
    end

    α
end


# Find a topological ordering of a directed graph by performing a descent-first search.
function dfs(graph::AbstractGraph)
    order = FixedStack{Int}(nv(graph))
    stack = FixedStack{Int}(nv(graph))
    color = zeros(Int, nv(graph))

    for v in vertices(graph)
        if color[v] == 0
            push!(stack, v)

            while !isempty(stack)
                u = last(stack)

                if color[u] == 0
                    color[u] = 1

                    for n in outneighbors(graph, u)
                        if color[n] == 0
                            push!(stack, n)
                        elseif color[n] == 1
                            error()
                        end
                    end
                else
                    pop!(stack)

                    if color[u] == 1
                        color[u] = 2
                        push!(order, u)
                    end
                end
            end
        end
    end

    order.items
end


# A Compact Row Storage Scheme for Cholesky Factors Using Elimination Trees
# Liu
# Algorithm 4.2: Elimination Tree by Path Compression.
function etree(graph::AbstractGraph)
    n = nv(graph)
    parent = Vector{Int}(undef, n)
    ancestor = Vector{Int}(undef, n)

    for i in 1:n
        parent[i] = 0
        ancestor[i] = 0

        for k in outneighbors(graph, i)
            r = k

            while !iszero(ancestor[r]) && ancestor[r] != i
                t = ancestor[r]
                ancestor[r] = i
                r = t
            end

            if iszero(ancestor[r])
                ancestor[r] = i
                parent[r] = i
            end
        end
    end

    Tree(parent)
end


# An Efficient Algorithm to Compute Row and Column Counts for Sparse Cholesky Factorization
# Gilbert, Ng, and Peyton
# Figure 3: Implementation of algorithm to compute row and column counts.
function supcnt(graph::AbstractGraph, tree::AbstractTree, level::AbstractVector=levels(tree), fdesc::AbstractVector=firstdescendants(tree))
    n = nv(tree)
    sets = DisjointSets(n)
    prev_p = zeros(Int, n)
    prev_nbr = zeros(Int, n)
    rc = ones(Int, n)
    wt = ones(Int, n)

    for p in 1:n - 1
        wt[parentindex(tree, p)] = 0
    end
    
    for p in 1:n - 1
        wt[parentindex(tree, p)] -= 1

        for u in inneighbors(graph, p)
            if fdesc[p] > prev_nbr[u]
                wt[p] += 1
                pp = prev_p[u]
                
                if iszero(pp)
                    rc[u] += level[p] - level[u]
                else
                    q = find!(sets, pp)
                    rc[u] += level[p] - level[q]
                    wt[q] -= 1
                end
    
                prev_p[u] = p
            end

            prev_nbr[u] = p
        end

        union!(sets, p, parentindex(tree, p))
    end

    cc = wt

    for p in 1:n - 1
        cc[parentindex(tree, p)] += cc[p]
    end

    rc, cc
end
