import Pkg
Pkg.activate(".")

using SONiA
using GLMakie
using Revise
using TickTock

# Read .k file of LS-PrePost
GEOM_NAME = "triangle.k"
nel, nnode, conn, coord = read_LS_PrePost(GEOM_NAME)
conn_tris, conn_quads = separate_conn(conn)


# Define Material and Problem Type:
MAT_NAME = "Custom"
PROB_TYPE = "Plane_Stress"

# Show Geometry:
PL("Geometry", conn_tris, conn_quads, coord, 0, false)

# Geometry Check
geoCheckQuads(coord, conn_quads)
geoCheckTris(coord, conn_tris)


# Defining Volume-box 1 (first and last index must be the same)
xv = [-10 210 210 -10 -10]
yv = [-10 -10 10 10 -10]

# Define first Dirichelet BC
BC_type = 11
BC_Dirichlet_1 = BC_Dirichlet(BC_box(xv,yv,coord),BC_type,0.0,0.0)
plotBC(BC_box(xv,yv,coord)[:,2], BC_box(xv,yv,coord)[:,3], :red)
lines!(vec(xv),vec(yv),linestyle=:solid,color = :red)

# Assemble all BC Dirichlet
all_BC_Dirichlet = vcat(BC_Dirichlet_1)


# BC Neumann
# Define type of BC and values
FORCE_TYPE = "Distributed"
NORMAL_FORCE = 1/100
TANGENTIAL_FORCE = 0

# Defining Volume-box 3 (first and last index must be the same)
xv = [90 190 210 100 90]
yv = [140 -10 -10 160 140]

# Define first Neumann BC
NeumannNodes_1 = BC_box(xv,yv,coord)

#Plot Arrows
neumann_edges_q = BC_Neumann_edges(conn_quads, NeumannNodes_1)
neumann_edges_t = BC_Neumann_edges(conn_tris, NeumannNodes_1)
neumann_edges = [neumann_edges_q; neumann_edges_t]
plotDistributed(neumann_edges, coord, NORMAL_FORCE*100)

# Plot Nnode Distributed
# plotBC(NeumannNodes_1[:,2],NeumannNodes_1[:,3],:blue)

# PLot Box
# lines!(vec(xv),vec(yv),linestyle=:solid,color = :blue)

# Calculate Force vector
FF_1_tris = BC_Neumann(coord, conn_tris, NeumannNodes_1, NORMAL_FORCE, TANGENTIAL_FORCE)
FF_1_quads = BC_Neumann(coord, conn_quads, NeumannNodes_1, NORMAL_FORCE, TANGENTIAL_FORCE)

# Assemble all BC Neumann (can be: FF_all = FF_1 + FF_2 + ...)
FF_all = FF_1_tris + FF_1_quads

## SOLVE SYSTEM
tick()
U, Ktot, Kstar, FF = elastSolver(conn_tris, conn_quads, coord, MAT_NAME, PROB_TYPE, all_BC_Dirichlet, FF_all)
tock()


######################################################################
# POST-PROCESSING
######################################################################

FACTOR = factorCalc(U)
println("FACTOR: ",FACTOR)

# Recover Displacements
ux,uy = splitU(U,coord)
utot = uTot(ux, uy)

# PL("Displacement", utot, conn_tris, conn_quads, defCoord(coord, ux, uy, FACTOR)) 

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

def = true
# Plot Results 
PL("Deformed Geometry", conn_tris, conn_quads, coord, defCoord(coord, ux, uy, FACTOR), true)
PL_FIELD("Displacement", utot, conn_tris, conn_quads, coord, defCoord(coord, ux, uy, FACTOR), def) 
PL_FIELD("Sigma X", avrsx, conn_tris, conn_quads, coord, defCoord(coord, ux, uy, FACTOR), def)  
PL_FIELD("Sigma Y", avrsy, conn_tris, conn_quads, coord, defCoord(coord, ux, uy, FACTOR), def)  
PL_FIELD("Tau XY", avrtxy, conn_tris, conn_quads, coord, defCoord(coord, ux, uy, FACTOR), def)  
PL_FIELD("Von Mises", vm, conn_tris, conn_quads, coord, defCoord(coord, ux, uy, FACTOR), def) 

println("..END")