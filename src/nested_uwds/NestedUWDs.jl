module NestedUWDs


using AbstractTrees
using Catlab.ACSetInterface
using Catlab.BasicGraphs
using Catlab.DirectedWiringDiagrams
using Catlab.DirectedWiringDiagrams: WiringDiagramACSet
using Catlab.MonoidalUndirectedWiringDiagrams
using Catlab.MonoidalUndirectedWiringDiagrams: UntypedHypergraphDiagram
using Catlab.RelationalPrograms
using Catlab.RelationalPrograms: TypedUnnamedRelationDiagram
using Catlab.Theories
using Catlab.UndirectedWiringDiagrams
using Catlab.WiringDiagramAlgebras

using ..JunctionTrees
using ..JunctionTrees: insertsorted!, DEFAULT_ELIMINATION_ALGORITHM, DEFAULT_SUPERNODE

# Elimination Algorithms
export EliminationAlgorithm, AMDJL_AMD, CuthillMcKeeJL_RCM, MetisJL_ND, MCS

# Supernodes
export Supernode, Node, MaximalSupernode, FundamentalSupernode

# Nested UWDs
export NestedUWD
export evalschedule, makeschedule, makeoperations


include("nested_uwds.jl")


end
