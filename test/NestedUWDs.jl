using Catlab.RelationalPrograms
using Catlab.UndirectedWiringDiagrams
using LinearAlgebra
using StructuredDecompositions.NestedUWDs
using Test


# CategoricalTensorNetworks.jl
# https://github.com/AlgebraicJulia/CategoricalTensorNetworks.jl/
function contract_tensor_network(d::UndirectedWiringDiagram,
                                 tensors::AbstractVector{<:AbstractArray})
    @assert nboxes(d) == length(tensors)
    juncs = [junction(d, ports(d, b)) for b in boxes(d)]
    j_out = junction(d, ports(d, outer=true), outer=true)
    contract_tensor_network(tensors, juncs, j_out)
end


function contract_tensor_network(tensors::AbstractVector{<:AbstractArray{T}},
                                 juncs::AbstractVector, j_out) where T
    # Handle important binary case with specialized code.
    if length(tensors) == 2 && length(juncs) == 2
        return contract_tensor_network(Tuple(tensors), Tuple(juncs), j_out)
    end

    jsizes = Tuple(infer_junction_sizes(tensors, juncs, j_out))
    juncs, j_out = map(Tuple, juncs), Tuple(j_out)
    C = zeros(T, Tuple(jsizes[j] for j in j_out))
    for index in CartesianIndices(jsizes)
        x = one(T)
        for (A, junc) in zip(tensors, juncs)
            x *= A[(index[j] for j in junc)...]
        end
        C[(index[j] for j in j_out)...] += x
    end
    return C
end


function contract_tensor_network( # Binary case.
    (A, B)::Tuple{<:AbstractArray{T},<:AbstractArray{T}},
    (jA, jB), j_out) where T
    jsizes = Tuple(infer_junction_sizes((A, B), (jA, jB), j_out))
    jA, jB, j_out = Tuple(jA), Tuple(jB), Tuple(j_out)
    C = zeros(T, Tuple(jsizes[j] for j in j_out))
    for index in CartesianIndices(jsizes)
        C[(index[j] for j in j_out)...] +=
            A[(index[j] for j in jA)...] * B[(index[j] for j in jB)...]
    end
    return C
end


function infer_junction_sizes(tensors, juncs, j_out)
    @assert length(tensors) == length(juncs)
    njunc = maximum(Iterators.flatten((Iterators.flatten(juncs), j_out)))
    jsizes = fill(-1, njunc)
    for (A, junc) in zip(tensors, juncs)
        for (i, j) in enumerate(junc)
            if jsizes[j] == -1
                jsizes[j] = size(A, i)
            else
                @assert jsizes[j] == size(A, i)
            end
        end
    end
    @assert all(s >= 0 for s in jsizes)
    jsizes
end


# out[v,z] = A[v,w] * B[w,x] * C[x,y] * D[y,z]    
diagram = @relation (v, z) begin
    A(v, w)
    B(w, x)
    C(x, y)
    D(y, z)
end

nuwd = NestedUWD(diagram)
A, B, C, D = map(randn, [(3, 4), (4, 5), (5, 6), (6, 7)])
generators = Dict{Symbol, Array{Float64}}(:A => A, :B => B, :C => C, :D => D)
out = evalschedule(contract_tensor_network, nuwd, generators)
@test out ≈ A * B * C * D

# out[] = A[w,x] * B[x,y] * C[y,z] * D[z,w]
diagram = @relation () begin
    A(w, x)
    B(x, y)
    C(y, z)
    D(z, w)
end

nuwd = NestedUWD(diagram)
A, B, C, D = map(randn, [(10, 5), (5, 5), (5, 5), (5, 10)])
generators = Dict{Symbol, Array{Float64}}(:A => A, :B => B, :C => C, :D => D)
out = evalschedule(contract_tensor_network, nuwd, generators)
@test out[] ≈ tr(A * B * C * D)

# out[w,x,y,z] = A[w,x] * B[y,z]
diagram = @relation (w, x, y, z) begin
    A(w, x)
    B(y, z)
end

nuwd = NestedUWD(diagram)
A, B = map(randn, [(3, 4), (5, 6)])
generators = Dict{Symbol, Array{Float64}}(:A => A, :B => B)
out = evalschedule(contract_tensor_network, nuwd, generators)
@test out ≈ (reshape(A, (3, 4, 1, 1)) .* reshape(B, (1, 1, 5, 6)))

# out[] = A[x,y] * B[x,y]
diagram = @relation () begin
    A(x, y)
    B(x, y)
end

nuwd = NestedUWD(diagram)
A, B = map(randn, [(5, 5), (5, 5)])
generators = Dict{Symbol, Array{Float64}}(:A => A, :B => B)
out = evalschedule(contract_tensor_network, nuwd, generators)
@test out[] ≈ dot(vec(A), vec(B))
