# A junction tree.
# This type implements the abstract graph and abstract tree interfaces.
struct JunctionTree <: AbstractVector{Bag}
    tree::Tree
    sndptr::Vector{Int}
    sepptr::Vector{Int}
    sepval::Vector{Int}
    relval::Vector{Int}
end


function Bag(tree::JunctionTree, i::Integer)
    Bag(residual(tree, i), separator(tree, i))
end


# Construct a junction tree. 
function junctiontree(matrix::SparseMatrixCSC, alg::Union{AbstractVector, EliminationAlgorithm}=DEFAULT_ELIMINATION_ALGORITHM, type::SupernodeType=DEFAULT_SUPERNODE_TYPE)
    label, index = permutation(matrix, alg)
    upper = triu(matrix)
    label, junctiontree!(upper, label, invsymperm(upper, index), type)
end


# Construct a junction tree. 
function junctiontree!(temporary::SparseMatrixCSC, label::AbstractVector, upper::SparseMatrixCSC, type::SupernodeType)
    lower, tree, sndptr, sepptr = supernodetree!(temporary, label, upper, type)
    sepval = sepvals(lower, tree, sndptr, sepptr)
    relval = relvals(tree, sndptr, sepptr, sepval)
    JunctionTree(tree, sndptr, sepptr, sepval, relval)
end


# Construct a postordered elimination tree.
function eliminationtree!(temporary::SparseMatrixCSC, label::AbstractVector, upper::SparseMatrixCSC)
    tree = etree(upper)
    index = dfs(tree)
    invpermute!(label, index)
    transpose!(upper, invsymperm!(temporary, upper, index)), invpermute!(tree, index)
end


# Construct a postordered supernodal elimination tree.
function supernodetree!(temporary::SparseMatrixCSC, label::AbstractVector, upper::SparseMatrixCSC, type::SupernodeType)
    lower, etree = eliminationtree!(temporary, label, upper)
    _, colcount = supcnt(lower, etree)
    new, ancestor, tree = stree(etree, colcount, type)
    index = dfs(tree)

    sndptr = Vector{Int}(undef, length(tree) + 1)
    sepptr = Vector{Int}(undef, length(tree) + 1)
    sndval = Vector{Int}(undef, size(lower, 1))
    sndptr[1] = sepptr[1] = 1

    for (i, j) in enumerate(invperm(index))
        v = new[j]
        p = sndptr[i]

        while !isnothing(v) && v != ancestor[j]
            sndval[v] = p
            v = parentindex(etree, v)
            p += 1
        end

        sndptr[i + 1] = p
        sepptr[i + 1] = sndptr[i] + sepptr[i] + colcount[new[j]] - p
    end

    invpermute!(label, sndval)
    invsymperm!(temporary, lower, sndval, ReverseOrdering()), invpermute!(tree, index), sndptr, sepptr
end


# Construct a postordered supernodal elimination tree.
function supernodetree!(temporary::SparseMatrixCSC, label::AbstractVector, upper::SparseMatrixCSC, type::Node)
    lower, tree = eliminationtree!(temporary, label, upper)
    _, colcount = supcnt(lower, tree)
    sndptr = OneTo(size(lower, 1) + 1)
    sepptr = cumsum(vcat(1, colcount .- 1))
    lower, tree, sndptr, sepptr
end


# Get the separators of every node of a supernodal elimination tree.
function sepvals(lower::SparseMatrixCSC, tree::Tree, sndptr::AbstractVector, sepptr::AbstractVector)
    temporary = zeros(Int, size(lower, 1))
    sepval = Vector{Int}(undef, last(sepptr) - 1)
    p = 1

    for j in tree
        for v in view(rowvals(lower), nzrange(lower, sndptr[j]))
            if sndptr[j + 1] <= v
                sepval[p] = v
                temporary[v] = j
                p += 1
            end
        end

        for i in childindices(tree, j), v in view(sepval, sepptr[i]:sepptr[i + 1] - 1)
            if sndptr[j + 1] <= v && temporary[v] != j
                sepval[p] = v
                temporary[v] = j
                p += 1
            end
        end

        sort!(sepval, sepptr[j], sepptr[j + 1] - 1, DEFAULT_STABLE, ForwardOrdering())
    end

    sepval
end


# Get the relative indices of every node of a supernodal elimination tree.
function relvals(tree::Tree, sndptr::AbstractVector, sepptr::AbstractVector, sepval::AbstractVector)
    relval = Vector{Int}(undef, length(sepval))
    p = 1

    for i in tree[1:end - 1]
        j = parentindex(tree, i)

        while p < sepptr[i + 1] && sepval[p] < sndptr[j + 1]
            relval[p] = sepval[p] - sndptr[j] + 1
            p += 1
        end

        q = sepptr[j]

        while p < sepptr[i + 1]
            q = searchsortedfirst(sepval, sepval[p], q, sepptr[j + 1] - 1, ForwardOrdering()) + 1
            relval[p] = q - sepptr[j] - sndptr[j] + sndptr[j + 1]
            p += 1
        end
    end

    relval
end


# Compute the width of a junction tree
function treewidth(tree::JunctionTree)
    maximum(length, tree) - 1
end


# Get the residual at node i.
function residual(tree::JunctionTree, i::Integer)
    tree.sndptr[i]:tree.sndptr[i + 1] - 1
end


# Get the separator at node i.
function separator(tree::JunctionTree, i::Integer)
    view(tree.sepval, tree.sepptr[i]:tree.sepptr[i + 1] - 1)
end


# Get the relative indices at node i.
function relative(tree::JunctionTree, i::Integer)
    view(tree.relval, tree.sepptr[i]:tree.sepptr[i + 1] - 1)
end


###########################
# Abstract Tree Interface #
###########################


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


function AbstractTrees.childindices(tree::JunctionTree, i::Integer)
    childindices(tree.tree, i)
end


function AbstractTrees.ParentLinks(::Type{IndexNode{JunctionTree, Int}})
    StoredParents()
end


function AbstractTrees.SiblingLinks(::Type{IndexNode{JunctionTree, Int}})
    StoredSiblings()
end


function AbstractTrees.NodeType(::Type{IndexNode{JunctionTree, Int}})
    HasNodeType()
end


function AbstractTrees.nodetype(::Type{IndexNode{JunctionTree, Int}})
    IndexNode{JunctionTree, Int}
end


#############################
# Abstract Vector Interface #
#############################


function Base.getindex(tree::JunctionTree, i::Integer)
    Bag(tree, i)
end


function Base.IndexStyle(::Type{JunctionTree})
    IndexLinear()
end


function Base.size(tree::JunctionTree)
    size(tree.tree)
end
