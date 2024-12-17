using MLStyle

#using Catlab.Syntax
using Catlab.Graphs, Catlab.Graphics
using Catlab.Theories, Catlab.CategoricalAlgebra

#Reverse direction of the edges
function to_graph(g::ReflexiveGraph)::Graph
    F = FinFunctor(Dict(:V => :V, :E => :E), Dict(:src => :src, :tgt => :tgt), SchGraph, SchReflexiveGraph)
    ΔF = DeltaMigration(F)
    return migrate(Graph, g, ΔF)
  end

# subobject classifier of graphs
Ω = to_graph(complete_graph(ReflexiveGraph, 2))
add_edge!(Ω, 1,1)

one = to_graph(ReflexiveGraph(1))

t = homomorphisms(one, Ω)[1]

#complete graphs
k(n) = to_graph(complete_graph(ReflexiveGraph, n))


#find all characteristic maps
sub_characters(g) = homomorphisms(g, Ω)

#find all subgraphs using the subobject classifier
function sub(g)
    map(f -> legs(pullback(t, f))[2], sub_characters(g))
end

function epis(g::Graph, h::Graph)
    homomorphisms(g,h, epic=true)
end

function epis(gs::Vector{Graph})
    ggs = vec(collect(Base.product(gs,gs)))
    the_tuples = filter(tup -> length(tup[2]) > 0 , [((g,h), epis(g,h)) for (g,h) in ggs])
    return Dict(the_tuples)
end


#=
    In: a tree T and an integer w
    Out: all codecompositions with shape T and width w.
    
    Alg: 
        for each node v of T choose a graph Gᵥ on at most w vertices
        for each edge e=tu in T, make a list of all cospans (using function epis above) of the form 
            Gₜ -->> ? <<-- Gᵤ 
        then any combination of such choices of elements in these lists determine a codecompostion of 
        width w and shape T.  
=#

#=
    In: integer n, width 
    Out: all codecompositions of width w whose shape is a tree on n nodes. 

    Alg: enumerate trees on n nodes, pass them into previous function
=#

K₃ = k(3)
subK₃ = sub(K₃)

epi_dict = epis(map(dom, subK₃))


J₃ = to_graph(ReflexiveGraph(3))
add_edge!(J₃, 1,2)
add_edge!(J₃, 2,3)
to_graphviz(J₃)

K₂ = to_graph(ReflexiveGraph(2))


hs = homomorphisms 