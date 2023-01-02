import Pkg
Pkg.activate(".")

using Revise
using SONiA
using GLMakie
using TickTock

GEOM_NAME = "quarterPlate_36el"
conn_tris, conn_quads, coord = readInput(GEOM_NAME)

# Define Material and Problem Type:
MAT_NAME = "Custom"
PROB_TYPE = "Plane_Stress"

# Create new figures:
PL("Geometry", conn_tris, conn_quads, coord)

## Define BCs
# Dirichlet:
# Defining Volume-box 1 (first and last index must be the same)
# Get nodes from Volume Box 1
xv = [ -0.1 0.1 0.1 -0.1 -0.1]
yv = [ 0. 0. 210. 210. 0.]

# Define BC_D 1
BC_type = 10
BC_Dirichlet_1 = BC_Dirichlet(BC_box(xv,yv,coord),BC_type,0.0,0.0)
plotBC(BC_box(xv,yv,coord)[:,2], BC_box(xv,yv,coord)[:,3], :red)

# Defining Volume-box 2 (first and last index must be the same)
# Get nodes from Volume Box 2
xv = [ 0.0 210.0 210.0 0.0 0.0]
yv = [ -0.1 -0.1 0.1 0.1 -0.1]
# Define BC_D 2
BC_type = 1
BC_Dirichlet_2 = BC_Dirichlet(BC_box(xv,yv,coord),BC_type,0.0,0.0)
plotBC(BC_box(xv,yv,coord)[:,2],BC_box(xv,yv,coord)[:,3],:red)

# All BC Dirichlet
all_BC_Dirichlet = vcat(BC_Dirichlet_1, BC_Dirichlet_2)

# BC Neumann
xv = [ 200.0-0.1 200.0+0.1 200.0+0.1 200.0-0.1 200.0-0.1]
yv = [ -0.1 -0.1 210. 210. -0.1]
NeumannNodes_1 = BC_box(xv,yv,coord)
plotBC(NeumannNodes_1[:,2],NeumannNodes_1[:,3],:blue)

FORCE_TYPE = "Distributed"
NORMAL_FORCE = 1/200
TANGENTIAL_FORCE = 0
FF_1 = BC_Neumann(coord, conn_quads, NeumannNodes_1, NORMAL_FORCE, TANGENTIAL_FORCE)

FF_all = vcat(FF_1)

tick()
U = elastSolver(conn_tris, conn_quads, coord, MAT_NAME, PROB_TYPE, all_BC_Dirichlet, FF_all)
tock()

## POST PROCESSING
ux,uy = splitU(U,coord)

FACTOR = 10
plotMesh(conn_tris, defCoord(coord, ux, uy, FACTOR), :solid, :black)
plotMesh(conn_quads, defCoord(coord, ux, uy, FACTOR), :solid, :black)

utot = uTot(ux, uy)


# Define new plot  
PL("Displacement", utot, conn_tris, conn_quads, coord)  
# p2 = Figure()
# ax = Axis(p2[1, 1], aspect=DataAspect(),title = "Utot")
# plotField(defCoord(coord, ux, uy, FACTOR), conn_quads, utot)
# plotMesh(conn_tris, defCoord(coord, ux, uy, FACTOR), :solid, :black)
# plotMesh(conn_quads, defCoord(coord, ux, uy, FACTOR), :solid, :black)
# GLMakie.activate!()
# display(GLMakie.Screen(), p2)


sx, sy, txy = stressCalcGP(U, conn_quads, coord, MAT_NAME, PROB_TYPE)

avrsx = avarageStress(sx, conn_quads, coord)
avrsy = avarageStress(sy, conn_quads, coord)
avrtxy = avarageStress(txy, conn_quads, coord)

PL("Sigma X", avrsx, conn_tris, conn_quads, coord)  
PL("Sigma Y", avrsy, conn_tris, conn_quads, coord)  
PL("Tau XY", avrtxy, conn_tris, conn_quads, coord)  

# Define new plot    
# p3 = Figure()
# ax = Axis(p3[1, 1], aspect=DataAspect(),title = "Sigma x")
# plotField(defCoord(coord, ux, uy, FACTOR), conn_quads, avrsx)
# plotMesh(conn_tris, defCoord(coord, ux, uy, FACTOR), :solid, :black)
# plotMesh(conn_quads, defCoord(coord, ux, uy, FACTOR), :solid, :black)
# GLMakie.activate!()
# display(GLMakie.Screen(), p3)

# avrsy = avarageStress(sy, conn_quads, coord)
# # Define new plot    
# p4 = Figure()
# ax = Axis(p4[1, 1], aspect=DataAspect(),title = "Sigma y")
# plotField(defCoord(coord, ux, uy, FACTOR), conn_quads, avrsy)
# plotMesh(conn_tris, defCoord(coord, ux, uy, FACTOR), :solid, :black)
# plotMesh(conn_quads, defCoord(coord, ux, uy, FACTOR), :solid, :black)
# GLMakie.activate!()
# display(GLMakie.Screen(), p4)

# avrtxy = avarageStress(txy, conn_quads, coord)
# # Define new plot    
# p5 = Figure()
# ax = Axis(p5[1, 1], aspect=DataAspect(),title = "Tau xy")
# plotField(defCoord(coord, ux, uy, FACTOR), conn_quads, avrtxy)
# plotMesh(conn_tris, defCoord(coord, ux, uy, FACTOR), :solid, :black)
# plotMesh(conn_quads, defCoord(coord, ux, uy, FACTOR), :solid, :black)
# GLMakie.activate!()
# display(GLMakie.Screen(), p5)


println("..END")