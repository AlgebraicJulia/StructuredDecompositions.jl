using ACSets, ACSets.ADTs
using Catlab
# using StructuredDecompositions

# specify a problem
struct Coloring <: AbstractFunctor
    n::Int
    f::Function
    function Coloring(n::Int)
        new(n, g -> homomorphisms(g, K(n)))
    end
end

# coloring is a representable functor Graph(-, K(n))
(c::Coloring)(x::Graph) = FinSet(c.f(x))
(c::Coloring)(f::ACSetTransformation) = FinFunction(λ -> f ∘ λ, c(codom(f)), c(dom(f)))
# notice the contravariance in action

# specify a decomposition

macro graph(head, body)
    v = 0;
    parsebody = @λ begin
        Expr(:block, args...) => parsebody.(args)
        Expr(:call, :E, s, t) => begin
            v = maximum([v, s, t])
            Expr(:call, :E, s, t)
        end
        :: LineNumberNode => nothing
        s => s
    end
    edges = parsebody(body)
    # result = quote
    #     construct(Graph, acsetspec(:Graph,
            quote
                $(Expr(:tuple, [Expr(:call, :V) for i in v])...)
                $(edges...)
            end
            # ))
    # end
end

G1 = construct(Graph, acsetspec(:(Graph), 
    quote
        V(); V(); V(); V(); V();
        E(src=1,tgt=2)
        E(src=2,tgt=3)
        E(src=3,tgt=4)
        E(src=3,tgt=5)
        E(src=4,tgt=5)
end))

macro see(head, body)
        dump(body)
end

# Bag 1
G1 = @acset Graph begin
    V=5; E=5
    src=[1,2,3,3,4]
    tgt=[2,3,4,5,5]
end
# A 12
G12 = @acset Graph begin V=1; end
# Bag 2
G2 = @acset Graph begin
    V=7
    E=12
    src=[1,2,3,4,5,6,1,7,7,7,7,2]
    tgt=[2,3,4,5,6,1,5,1,6,5,3,4]
end
# A 23
G23 = @acset Graph begin
    V=3
    E=2
    src=[1,2]
    tgt=[2,3]
end
# Bag 3
G3 = @acset Graph begin
    V=7
    E=6
    src=[1,2,3,5,4,4]
    tgt=[3,3,4,4,6,7]
end
# A 31
G31 = @acset Graph begin V=2; end

shape = @acset Graph begin
    V=3; E=3
    src=[1,2,3]
    tgt=[2,3,1]
end

Γ=FinDomFunctor(
    Dict(1=>G1,2=>G2,3=>G3,4=>G12,5=>G23,6=>G31),
    Dict(1=>ACSetTransformation(G12, G1, V=[1]),
         2=>ACSetTransformation(G12, G2, V=[6]),
         3=>ACSetTransformation(G23, G2, V=[1,6,5]), # error
         4=>ACSetTransformation(G23, G3, V=[6,4,5]), # error
         5=>ACSetTransformation(G31, G3, V=[1,2]),
         6=>ACSetTransformation(G31, G3, V=[3,4])
    ),
    ∫(shape));

