module SONiA

__precompile__(true)

# Include Modules:
include("InputManager.jl")
include("MyPlots.jl")
include("BCs.jl")
include("Solvers.jl")
include("PostProcessing.jl")

# Read and manage Input File:
using .InputManager
export readInput

# Plot functions:
using .MyPlots
export plotMesh, plotBC, plotField, PL

# BCs
using .BCs
export BC_box, BC_Dirichlet, BC_Neumann

# Problem Type:
using .Solvers
export elastSolver

# Post PROCESSING
using . PostProcessing
export splitU, defCoord, uTot, stressCalcGP, avarageStress

end
