"""
    NestedUWD{T, B, V}

An undirected wiring diagram, represented as a nested collected of undirected wiring
diagrams.
"""
struct NestedUWD{T, B, V}
    diagram::TypedUnnamedRelationDiagram{T, B, V}
    jtree::JunctionTree
    assignmentlist::Vector{Int}
    assignmentindex::Vector{Vector{Int}}
end


function NestedUWD(
    diagram::D,
    jtree::JunctionTree,
    assignmentlist::AbstractVector,
    assignmentindex::AbstractVector) where D <: UndirectedWiringDiagram

    T, B, V = getattributetypes(D)
    relation = TypedUnnamedRelationDiagram{T, B, V}()
    copy_parts!(relation, diagram)
    NestedUWD{T, B, V}(relation, jtree, assignmentlist, assignmentindex)
end


function NestedUWD(diagram::UndirectedWiringDiagram, jtree::JunctionTree)
    n = nparts(diagram, :Box)
    m = length(jtree)
    assignmentlist = Vector{Int}(undef, n)
    assignmentindex = Vector{Vector{Int}}(undef, m)
    assignmentindex .= [[]]
    
    for b in 1:n
        i = getsubtree(jtree, diagram[incident(diagram, b, :box), :junction])
        assignmentlist[b] = i
        push!(assignmentindex[i], b)
    end

    NestedUWD(diagram, jtree, assignmentlist, assignmentindex)
end 


"""
    NestedUWD(
        diagram::UndirectedWiringDiagram,
        [, algorithm::Union{Order, EliminationAlgorithm}]
        [, supernode::Supernode])

Construct a nested undirected wiring diagram.
"""
function NestedUWD(
    diagram::UndirectedWiringDiagram,
    algorithm::Union{Order, EliminationAlgorithm}=DEFAULT_ELIMINATION_ALGORITHM,
    supernode::Supernode=DEFAULT_SUPERNODE)

    jtree = JunctionTree(diagram, algorithm, supernode)
    NestedUWD(diagram, jtree)
end


# Construct a tree decomposition of the line graph of an undirected wiring diagram.
function JunctionTree(
    diagram::UndirectedWiringDiagram,
    algorithm::Union{Order, EliminationAlgorithm}=DEFAULT_ELIMINATION_ALGORITHM,
    supernode::Supernode=DEFAULT_SUPERNODE)

    graph = makegraph(diagram)
    jtree = JunctionTree(graph, algorithm, supernode)

    query = diagram[:outer_junction]
    JunctionTree(getsubtree(jtree, query), jtree)
end


# Construct the line graph of an undirected wiring diagram.
function makegraph(diagram::UndirectedWiringDiagram)
    n = nparts(diagram, :Junction)
    m = nparts(diagram, :Box)
    graph = SymmetricGraph(n)

    for i in 1:m
        junctions = diagram[incident(diagram, i, :box), :junction]
        l = length(junctions)
        
        for i₁ in 1:l - 1
            j₁ = junctions[i₁]

            for i₂ in i₁ + 1:l
                j₂ = junctions[i₂]

                if !has_edge(graph, j₁, j₂)
                    add_edge!(graph, j₁, j₂)
                end
            end
        end
    end

    junctions = diagram[:, :outer_junction]
    l = length(junctions)
    
    for i₁ in 1:l - 1
        j₁ = junctions[i₁]

        for i₂ in i₁ + 1:l
            j₂ = junctions[i₂]

            if !has_edge(graph, j₁, j₂)
                add_edge!(graph, j₁, j₂)
            end
        end
    end

    graph
end


"""
    makeschedule(nuwd::NestedUWD)

Construct a directed wiring diagram that represents the nesting structure of a nested UWD.
"""
function makeschedule(nuwd::NestedUWD{<:Any, T}) where T
    m = length(nuwd.assignmentlist)
    n = length(nuwd.jtree)

    parents = map(1:n - 1) do i
        parentindex(nuwd.jtree, i)
    end

    costs = map(1:n) do i
        length(getresidual(nuwd.jtree, i)) + length(getseperator(nuwd.jtree, i))
    end

    schedule = WiringDiagramACSet{T, Nothing, Union{Int, AbstractBox}, DataType}()

    add_parts!(schedule, :Box, n)
    add_parts!(schedule, :Wire, n - 1)
    add_parts!(schedule, :InPort, m + n - 1)
    add_parts!(schedule, :InWire, m)
    add_parts!(schedule, :OutPort, n)
    add_parts!(schedule, :OutWire, 1)
    add_parts!(schedule, :OuterInPort, m)
    add_parts!(schedule, :OuterOutPort, 1)

    schedule[:, :src] = 1:n - 1
    schedule[:, :tgt] = m + 1:m + n - 1
    schedule[:, :in_src] = 1:m
    schedule[:, :in_tgt] = 1:m
    schedule[:, :out_src] = n:n
    schedule[:, :out_tgt] = 1:1
    schedule[:, :in_port_box] = [nuwd.assignmentlist; parents]
    schedule[:, :out_port_box] = 1:n
    
    schedule[:, :value] = costs
    schedule[:, :box_type] = Box{Int}
    schedule[:, :outer_in_port_type] = nuwd.diagram[:, :name]

    Theory = ThSymmetricMonoidalCategory.Meta.T
    WiringDiagram{Theory, T, Nothing, Int}(schedule, nothing)
end


"""
    function evalschedule(
        f,
        nuwd::NestedUWD,
        generators::Union{AbstractVector, AbstractDict}
        [, operations::AbstractVector])

Evaluate an undirected wiring diagrams given a set of generators for the boxes. The
optional first argument `f` should be callable with the signature
```
    f(diagram, generators)
```
where `diagram` is an undirected wiring diagram, and `generators` is a vector. If `f` is not
specified, then it defaults to `oapply`.
"""
function evalschedule(
    f,
    nuwd::NestedUWD,
    generators::AbstractVector{T},
    operations::AbstractVector=makeoperations(nuwd)) where T

    n = length(nuwd.jtree)
    mailboxes = Vector{T}(undef, n)

    for i in 1:n
        g₁ = generators[nuwd.assignmentindex[i]]
        g₂ = mailboxes[childindices(nuwd.jtree, i)]
        mailboxes[i] = f(operations[i], [g₁; g₂])
    end
    
    mailboxes[n]
end


function evalschedule(
    f,
    nuwd::NestedUWD,
    generators::AbstractDict{<:Any, T},
    operations::AbstractVector=makeoperations(nuwd)) where T

    g = generators
    n = nparts(nuwd.diagram, :Box)
    generators = Vector{T}(undef, n)

    for i in 1:n
        generators[i] = g[nuwd.diagram[i, :name]]
    end

    evalschedule(f, nuwd, generators, operations)
end


function evalschedule(
    nuwd::NestedUWD,
    generators::Union{AbstractVector, AbstractDict},
    operations::AbstractVector=makeoperations(nuwd))

    evalschedule(oapply, nuwd, generators, operations)
end


# For each node i of a nested UWD, construct the undirected wiring diagram corresponding to i.
function makeoperations(nuwd::NestedUWD)
    m = length(nuwd.jtree)

    map(1:m) do i
        makeoperation(nuwd, i)
    end
end


# Construct the undirected wiring diagram corresponding to node i of a nested UWD.
function makeoperation(nuwd::NestedUWD{T, B, V}, i::Integer) where {T, B, V}
    function findjunction(j::Integer)
        v = nuwd.jtree.stree.order.index[j]
        v₁ = nuwd.jtree.stree.firstsupernodelist[i]
        v₂ = nuwd.jtree.stree.lastsupernodelist[i]

        if v <= v₂
            v - v₁ + 1
        else
            v₂ - v₁ + 1 + searchsortedfirst(nuwd.jtree.seperatorlist[i], v)
        end
    end

    residual = getresidual(nuwd.jtree, i)
    seperator = getseperator(nuwd.jtree, i)
    m = length(residual)
    n = length(seperator)

    operation = TypedUnnamedRelationDiagram{T, B, V}()
    add_parts!(operation, :Junction, m + n)

    operation[1:m, :junction_type] = nuwd.diagram[residual, :junction_type]
    operation[1:m, :variable] = nuwd.diagram[residual, :variable]
    operation[m + 1:m + n, :junction_type] = nuwd.diagram[seperator, :junction_type]
    operation[m + 1:m + n, :variable] = nuwd.diagram[seperator, :variable]

    if i < length(nuwd.jtree)
        for j in seperator
            p′ = add_part!(operation, :OuterPort)
            operation[p′, :outer_junction] = m + p′
            operation[p′, :outer_port_type] = nuwd.diagram[j, :junction_type]
        end
    else
        for j in nuwd.diagram[:outer_junction]
            p′ = add_part!(operation, :OuterPort)
            operation[p′, :outer_junction] = findjunction(j)
            operation[p′, :outer_port_type] = nuwd.diagram[j, :junction_type]
        end
    end

    for b in nuwd.assignmentindex[i]
        b′ = add_part!(operation, :Box)
        operation[b′, :name] = nuwd.diagram[b, :name]

        for j in nuwd.diagram[incident(nuwd.diagram, b, :box), :junction]
            p′ = add_part!(operation, :Port)
            operation[p′, :box] = b′
            operation[p′, :junction] = findjunction(j)
            operation[p′, :port_type] = nuwd.diagram[j, :junction_type]
        end
    end

    for b in childindices(nuwd.jtree, i)
        b′ = add_part!(operation, :Box)

        for j in getseperator(nuwd.jtree, b)
            p′ = add_part!(operation, :Port)
            operation[p′, :box] = b′
            operation[p′, :junction] = findjunction(j)
            operation[p′, :port_type] = nuwd.diagram[j, :junction_type]
        end
    end

    operation
end


# Get the attribute types of an undirected wiring diagram.
function getattributetypes(::Type{<:UntypedRelationDiagram{B, V}}) where {B, V}
    Nothing, B, V
end


function getattributetypes(::Type{<:TypedRelationDiagram{T, B, V}}) where {T, B, V}
    T, B, V
end
