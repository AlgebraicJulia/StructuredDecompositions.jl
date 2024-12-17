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

K₃ = k(3)
subK₃ = sub(K₃)

epi_dict = epis(map(dom, subK₃))


J₃ = to_graph(ReflexiveGraph(3))
add_edge!(J₃, 1,2)
add_edge!(J₃, 2,3)
to_graphviz(J₃)

K₂ = to_graph(ReflexiveGraph(2))


hs = homomorphisms 