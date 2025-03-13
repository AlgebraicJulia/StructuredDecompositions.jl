module StructuredDecompositions


using Reexport


include("decompositions/Decompositions.jl")
include("functor_utils/FunctorUtils.jl")
include("deciding_sheaves/DecidingSheaves.jl")


@reexport using .Decompositions
@reexport using .FunctorUtils
@reexport using .DecidingSheaves


end 
