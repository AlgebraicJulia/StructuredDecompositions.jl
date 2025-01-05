"""
    EliminationAlgorithm

A graph elimination algorithm. The options are
- [`CuthillMcKeeJL_RCM`](@ref)
- [`SymRCMJL_RCM`](@ref)
- [`AMDJL_AMD`](@ref)
- [`AMDJL_SYMAMD`](@ref)
- [`MetisJL_ND`](@ref)
- [`TreeWidthSolverJL_BT`](@ref)
- [`MCS`](@ref)
"""
abstract type EliminationAlgorithm end


"""
    CuthillMcKeeJL_RCM <: EliminationAlgorithm

The reverse Cuthill-McKee algorithm. Uses CuthillMckee.jl.
"""
struct CuthillMcKeeJL_RCM <: EliminationAlgorithm end


"""
    SymRCMJL_RCM <: EliminationAlgorithm

The reverse Cuthill-McKee algorithm. Uses SymRCM.jl.
"""
struct SymRCMJL_RCM <: EliminationAlgorithm end


"""
    AMDJL_AMD <: EliminationAlgorithm

The approximate minimum degree algorithm. Uses AMD.jl.
"""
struct AMDJL_AMD <: EliminationAlgorithm
    meta::AMD.Amd
end

"""
    AMDJL_SYMAMD <: EliminationAlgorithm

The SYMAMD algorithm. Uses AMD.jl.
"""
struct AMDJL_SYMAMD{T} <: EliminationAlgorithm
    meta::AMD.Colamd{T}
end


"""
    MetisJL_ND <: EliminationAlgorithm

The nested dissection heuristic. Uses Metis.jl.
"""
struct MetisJL_ND <: EliminationAlgorithm end


"""
    TreeWidthSolverJL_BT <: EliminationAlgorithm

The Bouchitte-Todinca algorithm. Uses TreeWidthSolver.jl.
"""
struct TreeWidthSolverJL_BT <: EliminationAlgorithm end


"""
    MCS <: EliminationAlgorithm

The maximum cardinality search algorithm.
"""
struct MCS <: EliminationAlgorithm end


function AMDJL_AMD()
    AMDJL_AMD(AMD.Amd())
end


function AMDJL_SYMAMD{T}() where T
    AMDJL_SYMAMD(AMD.Colamd{T}())
end


function AMDJL_SYMAMD()
    AMDJL_SYMAMD{Int}()
end


const DEFAULT_ELIMINATION_ALGORITHM = AMDJL_AMD()
