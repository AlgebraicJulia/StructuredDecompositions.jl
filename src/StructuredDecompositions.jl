module StructuredDecompositions

using Reexport

include("Decompositions.jl")
# include("FunctorUtils.jl") XXX
include("DecidingSheaves.jl")

@reexport using .Decompositions
@reexport using .DecidingSheaves

include("JunctionTrees.jl")
# include("StrTreeDecomp.jl")

end 
