module SONiA

__precompile__(true)

# Include Modules:
include("PreProcessing.jl")
include("MyPlots.jl")
include("BCs.jl")
include("Solvers.jl")
include("PostProcessing.jl")

# Read and manage Input File:
using .PreProcessing
export readInput, read_LS_PrePost, separate_conn, geoCheckTris, geoCheckQuads

# Plot functions:
using .MyPlots
export plotMesh, plotBC, plotField, PL, PL_FIELD

# BCs
using .BCs
export BC_box, BC_Dirichlet, BC_Neumann, arch, line

# Problem Type:
using .Solvers
export elastSolver

# Post PROCESSING
using . PostProcessing
export splitU, defCoord, uTot, stressCalcGP_tris, stressCalcGP_quads, avarageStress, 
vmStress, nodalStresses, factorCalc

end
