module StructuredDecompositions

using Reexport

include("Decompositions.jl")
include("FunctorUtils.jl")
include("DecidingSheaves.jl")
include("junction_trees/JunctionTrees.jl")
include("nested_uwds/NestedUWDs.jl")

@reexport using .Decompositions
@reexport using .FunctorUtils
@reexport using .DecidingSheaves

end 
