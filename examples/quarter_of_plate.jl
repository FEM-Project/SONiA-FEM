import Pkg
Pkg.activate(".")

using SONiA
using GLMakie
using Revise


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

# Show Geometry:
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

# Assemble all BC Neumann (can be: FF_all = FF_1 + FF_2 + ...)
FF_all = FF_1


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

# Recover stresses at Gauss Points:
sx_t, sy_t, txy_t = stressCalcGP_tris(U, conn_tris, coord, MAT_NAME, PROB_TYPE)
sx_q, sy_q, txy_q = stressCalcGP_quads(U, conn_quads, coord, MAT_NAME, PROB_TYPE)

# Recover stresses at nodes Triangles:
sx_n_t = nodalStresses(sx_t, conn_tris, coord)
sy_n_t = nodalStresses(sy_t, conn_tris, coord)
txy_n_t = nodalStresses(txy_t, conn_tris, coord)

# Recover stresses at nodes Quadrangles:
sx_n_q = nodalStresses(sx_q, conn_quads, coord)
sy_n_q = nodalStresses(sy_q, conn_quads, coord)
txy_n_q = nodalStresses(txy_q, conn_quads, coord)

# Avarage stresses
avrsx = avarageStress(sx_n_t, sx_n_q, conn_tris, conn_quads, coord)
avrsy = avarageStress(sy_n_t, sy_n_q, conn_tris, conn_quads, coord)
avrtxy = avarageStress(txy_n_t, txy_n_q, conn_tris, conn_quads, coord)

# Von Mises
vm = vmStress(avrsx,avrsy,avrtxy)

# Plot Results 
PL("Deformed Geometry", conn_tris, conn_quads, defCoord(coord, ux, uy, FACTOR))
PL("Displacement", utot, conn_tris, conn_quads, defCoord(coord, ux, uy, FACTOR)) 
PL("Sigma X", avrsx, conn_tris, conn_quads, coord)  
PL("Sigma Y", avrsy, conn_tris, conn_quads, coord)  
PL("Tau XY", avrtxy, conn_tris, conn_quads, coord)  
PL("Von Mises", vm, conn_tris, conn_quads, coord)  

println("..END")