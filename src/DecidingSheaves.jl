module DecidingSheaves

export Presheaf, Sheaf, decide

using ..Decompositions

#####################
#   COMPUTATION
#####################
"""(Pre)-Sheaves
"""
abstract type Presheaf end

abstract type Sheaf <: Presheaf end

"""
Filtering algorithm. 
INPUT: a Finset^{op}-valued structured decomposition d : FG → Finset^{op} 
       and an edge e = xy in G
OUTPUT: a structured decomposition obtained by replacing the span de in d 
        by the span obtained by projecting the pullback of de (i.e. taking images)
"""
function adhesion_filter(d::StructuredDecomposition, x,y)
  
end

#given a decomposition and a list Λ of edges of the decomposition shape
#run adhesion_filter for each edge in the list
function adhesion_filter(d::StructuredDecomposition, Λ)
end

"""Solve the decision problem encoded by a sheaf. 
The algorithm is as follows: 
  compute on each bag,
  compute composites on edges, 
  project back down to bags
  answer 
    "no" if there is an empty bag; 
    "yes" otherwise.
"""
function decide(f::Sheaf, d::StructuredDecomposition)::Bool
end

end