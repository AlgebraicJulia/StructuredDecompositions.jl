module StructuredDecompositions


using Reexport


include("junction_trees/JunctionTrees.jl")
include("decompositions/Decompositions.jl")
include("functor_utils/FunctorUtils.jl")
include("deciding_sheaves/DecidingSheaves.jl")


@reexport using .JunctionTrees


end 
