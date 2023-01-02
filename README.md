# SONiA

SONiA is a SOlver for Numerical Aplications written in Julia Programming Language to enhance knowledge on this topic

## How to create/read an example file?

#### 0. Write a geometry file as follow:
* first of all the number of nodes
* you need a connectivity matrix where you define the element number, the number of nodes and the nodes' order
* now you need node numbering and coordinates for nodes
**EXAMPLE:**
1    4   1     1     2     3     4  
2    4   1     4     3     6     5  
1         50.0         0.0    
2        200.0         0.0    
3        200.0       200.0    
4         35.3553     35.3553  
5          0.0        50.0   
6          0.0       200.0  

#### 1. Create a new julia file and activate SONiA package
   
```julia
import Pkg
Pkg.activate(".")
   
using SONiA 
```

#### 2. Select and read your geometry file

```julia
GEOM_NAME = "foo"
conn_tris, conn_quads, coord = readInput(GEOM_NAME) 
```

you now have connectivity for triangles (conn_tris), for quadrangles (conn_quads) and the coordinates matrix (coord)

#### 3. Define Material and problem type

```julia
MAT_NAME = "Custom"
PROB_TYPE = "Plane_Stress"
```

you can both create your material using "Custom" option and defining a linear-elastic material or use a meteria in MatLib.jl

#### 4. Dirichlet BCs
* create a box to get the nodes you want to constaint
```julia
xv = [ -0.1 0.1 0.1 -0.1 -0.1]
yv = [ 0. 0. 210. 210. 0.]
```

* define the type of constraint
```julia
BC_type = 10
```

* Apply BC
```julia
BC_Dirichlet_1 = BC_Dirichlet(BC_box(xv,yv,coord),BC_type,0.0,0.0)
```
here you can also define a prescribed displacement (in the example is set to 0.0)

* Get all Dirichlet BCs
```julia
all_BC_Dirichlet = vcat(BC_Dirichlet_1, BC_Dirichlet_2,...)
```

#### 5. Repeat step **4.** for all Dirichlet BCs

#### 6. Neumann BCs (here things becomes hard!)
* Define Nodal Forces
...

* Define Distributed Forces
```julia
FORCE_TYPE = "Distributed"
NORMAL_FORCE = 1
TANGENTIAL_FORCE = 1
```

* create a box to get the nodes you want to constaint
```julia
xv = [ 200.0-0.1 200.0+0.1 200.0+0.1 200.0-0.1 200.0-0.1]
yv = [ -0.1 -0.1 210. 210. -0.1]
```

* Apply BC
```julia
FF_1 = BC_Neumann(coord, conn_quads, BC_box(xv,yv,coord), NORMAL_FORCE, TANGENTIAL_FORCE)
```

* Add to total Neumann BCs
```julia
FF_all = vcat(FF_1)
``` 

#### 7. Plot the situation till here
* ...

#### 8. SOLVE
* Call the solver you need to get the resutl
```julia
U = elastSolver(conn_tris, conn_quads, coord, MAT_NAME, PROB_TYPE, all_BC_Dirichlet, FF_all)
```

---
---

#### 9. Post Processing
* Define multiplier factor
```julia
FACTOR = 10
```