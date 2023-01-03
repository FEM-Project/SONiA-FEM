import Pkg
Pkg.activate(".")

using SONiA
using GLMakie


######################################################################
# GEOMETRY
######################################################################

# There are 3 files for this geometry (*_2el, *_36el, *_400el)
GEOM_NAME = "quarterPlate_36el"
conn_tris, conn_quads, coord = readInput(GEOM_NAME)


######################################################################
# MATERIAL
######################################################################

# Define Material and Problem Type:
MAT_NAME = "Custom"
PROB_TYPE = "Plane_Stress"

# Show Geometru:
PL("Geometry", conn_tris, conn_quads, coord)


######################################################################
# BCs
######################################################################

## Define BCs
# Dirichlet:
# Defining Volume-box 1 (first and last index must be the same)
xv = [ -0.1 0.1 0.1 -0.1 -0.1]
yv = [ 0. 0. 210. 210. 0.]

# Define first Dirichelet BC
BC_type = 10
BC_Dirichlet_1 = BC_Dirichlet(BC_box(xv,yv,coord),BC_type,0.0,0.0)
plotBC(BC_box(xv,yv,coord)[:,2], BC_box(xv,yv,coord)[:,3], :red)

# Defining Volume-box 2 (first and last index must be the same)
xv = [ 0.0 210.0 210.0 0.0 0.0]
yv = [ -0.1 -0.1 0.1 0.1 -0.1]

# Define second Dirichelet BC
BC_type = 1
BC_Dirichlet_2 = BC_Dirichlet(BC_box(xv,yv,coord),BC_type,0.0,0.0)
plotBC(BC_box(xv,yv,coord)[:,2],BC_box(xv,yv,coord)[:,3],:red)

# Assemble all BC Dirichlet
all_BC_Dirichlet = vcat(BC_Dirichlet_1, BC_Dirichlet_2)

# BC Neumann
# Define type of BC and values
FORCE_TYPE = "Distributed"
NORMAL_FORCE = 1/200
TANGENTIAL_FORCE = 0

# Defining Volume-box 3 (first and last index must be the same)
xv = [ 200.0-0.1 200.0+0.1 200.0+0.1 200.0-0.1 200.0-0.1]
yv = [ -0.1 -0.1 210. 210. -0.1]

# Define first Neumann BC
NeumannNodes_1 = BC_box(xv,yv,coord)
plotBC(NeumannNodes_1[:,2],NeumannNodes_1[:,3],:blue)

# Calculate Force vector
FF_1 = BC_Neumann(coord, conn_quads, NeumannNodes_1, NORMAL_FORCE, TANGENTIAL_FORCE)

# Assemble all BC Neumann
FF_all = vcat(FF_1)


######################################################################
# SOLVE SYSTEM
######################################################################

## SOLVE SYSTEM
U = elastSolver(conn_tris, conn_quads, coord, MAT_NAME, PROB_TYPE, all_BC_Dirichlet, FF_all)


######################################################################
# POST-PROCESSING
######################################################################

FACTOR = 10

# Recover Displacements
ux,uy = splitU(U,coord)
utot = uTot(ux, uy)

# Recover stresses
sx, sy, txy = stressCalcGP(U, conn_quads, coord, MAT_NAME, PROB_TYPE)
avrsx = avarageStress(sx, conn_quads, coord)
avrsy = avarageStress(sy, conn_quads, coord)
avrtxy = avarageStress(txy, conn_quads, coord)

# Plot Results 
PL("Displacement", utot, conn_tris, conn_quads, coord) 
PL("Sigma X", avrsx, conn_tris, conn_quads, coord)  
PL("Sigma Y", avrsy, conn_tris, conn_quads, coord)  
PL("Tau XY", avrtxy, conn_tris, conn_quads, coord)  

println("..END")