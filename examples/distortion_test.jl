import Pkg
Pkg.activate(".")

using SONiA
using GLMakie
using Revise


######################################################################
# GEOMETRY
######################################################################

# There are 3 files for this geometry (*_2el, *_36el, *_400el)
GEOM_NAME = "distortion_test"
conn_tris, conn_quads, coord = readInput(GEOM_NAME)


######################################################################
# MATERIAL
######################################################################

# Define Material and Problem Type:
MAT_NAME = "Custom"
PROB_TYPE = "Plane_Stress"

# Show Geometry:
PL("Geometry", conn_tris, conn_quads, coord, coord, false)

# Geometry Check
geoCheckQuads(coord, conn_quads)