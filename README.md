# SONiA-FEM

SONiA is a SOlver for Numerical Aplications written in Julia Programming Language to enhance knowledge on this topic

<img title="" src="doc//Fig.png" alt="" width="735">

## How to use it

* Download and install Julia Programming language for your platform: https://julialang.org/
* Download and install Visual Studio Code for your platform: https://code.visualstudio.com/
* Follow the steps to set up Julia in VS Code: https://code.visualstudio.com/docs/languages/julia

## How to create/read an example file?

#### 0. Write a geometry file as follow:

**CASE 1** *- from scratches*

* first of all the number of elements

* you need a connectivity matrix where you define the element number (1,2,..,n), the number of nodes (3 or 4), the associated material (1,2,..,n) and the nodes' order.
  
  * Right now the solver can handle just 1 material, but the connectivity matrix is generalized to implement multiple materials

* now you need node numbering and coordinates for nodes
  **EXAMPLE:**
  
  ```
  2
  1 4 1 1 2 3 4  
  2 4 1 4 3 6 5  
  1 50.0 0.0    
  2 200.0 0.0    
  3 200.0 200.0    
  4 35.3553 35.3553  
  5 0.0 50.0   
  6 0.0 200.0  
  ```

**CASE 2** *- using LS-PREPOST*

* Simply use LS-PREPOST to create a *.k* file

#### 1. Create a new Julia file and activate SONiA package

```julia
import Pkg
Pkg.activate(".")

using SONiA 
```

#### 2. Select and read your geometry file

**CASE 1** *- your geometry*

```julia
GEOM_NAME = "foo"
conn_tris, conn_quads, coord = readInput(GEOM_NAME) 
```

**CASE 2** *- LS-PREPOST geometry*

```julia
GEOM_NAME = "foo.k"
nel, nnode, conn, coord = read_LS_PrePost(GEOM_NAME)
conn_tris, conn_quads = separate_conn(conn)
```

you now have connectivity for triangles (conn_tris), for quadrangles (conn_quads) and the coordinates matrix (coord)

why not? Plot your geometry using PL function

```julia
PL("Geometry", conn_tris, conn_quads, coord, 0, false)
```

#### 3. Define Material and problem type

```julia
MAT_NAME = "Steel"
PROB_TYPE = "Plane_Stress"
```

All materials are defined in the folder "materials" and are linear-elastic. You can choose to use custom linear-elastic material modifying the "Custom.txt" file or adding new material files.
IMPORTANT: only linear-elastic material work!

#### 4. Dirichlet BCs

* create a box to get the nodes you want to constraint
  
  ```julia
  xv = [ -0.1 0.1 0.1 -0.1 -0.1]
  yv = [ 0. 0. 210. 210. 0.]
  ```

* define the type of constraint
  
  * 00 = free
  
  * 01 = free x, constraint y
  
  * 10 = constraint x, free y
  
  * 11 = constraint x and y
  
  ```julia
  BC_type = 10
  ```

* Apply BC
  
  ```julia
  BC_Dirichlet_1 = BC_Dirichlet(BC_box(xv,yv,coord),BC_type,0.0,0.0)
  ```
  
  here you can also define a prescribed displacement (in the example is set to 0.0) 
  If you like you can plot a different color (:red) for the selected nodes using
  
  ```julia
  plotBC(BC_box(xv,yv,coord)[:,2], BC_box(xv,yv,coord)[:,3], :red)
  ```

* Get all Dirichlet BCs
  
  ```julia
  all_BC_Dirichlet = vcat(BC_Dirichlet_1, BC_Dirichlet_2,...)
  ```

* Repeat step **4.** for all Dirichlet BCs

#### 5. Neumann BCs (here things becomes hard!)

* Define Nodal Forces
  ...

* Define Distributed Forces
  
  ```julia
  FORCE_TYPE = "Distributed"
  NORMAL_FORCE = 1
  TANGENTIAL_FORCE = 1
  ```

* create a box to get the nodes you want to constraint
  
  ```julia
  xv = [ -0.1 0.1 0.1 -0.1 -0.1]
  yv = [ 0. 0. 210. 210. 0.]
  ```

* Get Neumann nodes
  
  ```julia
  NeumannNodes_1 = BC_box(xv,yv,coord)
  ```

* Get Neumann edges for quadrangles (q) and triangles (t)
  
  ```julia
  neumann_edges_q = BC_Neumann_edges(conn_quads, NeumannNodes_1)
  neumann_edges_t = BC_Neumann_edges(conn_tris, NeumannNodes_1)
  ```

* Plot arrows on edges for distributed forces
  
  ```julia
  neumann_edges = [neumann_edges_q; neumann_edges_t]
  plotDistributed(neumann_edges, coord, NORMAL_FORCE*100)
  ```

* Calculate forces vector for both triangles and quadrangle and assemble total **Force Vector**
  
  ```julia
  FF_1_tris = BC_Neumann(coord, conn_tris, NeumannNodes_1, NORMAL_FORCE, TANGENTIAL_FORCE)
  FF_1_quads = BC_Neumann(coord, conn_quads, NeumannNodes_1, NORMAL_FORCE, TANGENTIAL_FORCE)
  FF_all = FF_1_tris + FF_1_quads + ...
  ```

#### 6. SOLVE LINEAR SYSTEM

* Call the solver you need to get the results
  
  ```julia
  U, Ktot, Kstar, FF = elastSolver(conn_tris, conn_quads, coord, MAT_NAME, PROB_TYPE, all_BC_Dirichlet, FF_all)
  ```

---

---

#### 7. Post Processing

* Automatically define multiplier factor
  
  ```julia
  FACTOR = factorCalc(U)
  ```

* Recover fields to plot
  
  * separate ux and uy using the function splitU
    
    ```julia
    ux,uy = splitU(U,coord)
    ```
  
  * calculate Utot
    
    ```julia
    utot = uTot(ux, uy)
    ```
  
  * calculate stresses at GPs positions for triangles and quadrangles
    
    ```julia
    sx_t, sy_t, txy_t = stressCalcGP_tris(U, conn_tris, coord, MAT_NAME, PROB_TYPE)
    sx_q, sy_q, txy_q = stressCalcGP_quads(U, conn_quads, coord, MAT_NAME, PROB_TYPE)
    ```
  
  * calculate stresses at nodes positions for triangles and quadrangles
    
    ```julia
    sx_n_t = nodalStresses(sx_t, conn_tris, coord)
    sy_n_t = nodalStresses(sy_t, conn_tris, coord)
    txy_n_t = nodalStresses(txy_t, conn_tris, coord)
    ```
    
    ```julia
    sx_n_q = nodalStresses(sx_q, conn_quads, coord)
    sy_n_q = nodalStresses(sy_q, conn_quads, coord)
    txy_n_q = nodalStresses(txy_q, conn_quads, coord)
    ```
  
  * calculate average stresses
    
    ```julia
    avrsx = avarageStress(sx_n_t, sx_n_q, conn_tris, conn_quads, coord)
    avrsy = avarageStress(sy_n_t, sy_n_q, conn_tris, conn_quads, coord)
    avrtxy = avarageStress(txy_n_t, txy_n_q, conn_tris, conn_quads, coord)
    ```
  
  * calculate Von Mises stresses
    
    ```julia
    vm = vmStress(avrsx,avrsy,avrtxy)
    ```

* Plot the result of the deformed geometry using PL
  
  ```julia
  PL("Deformed Geometry", conn_tris, conn_quads, coord, defCoord(coord, ux, uy, FACTOR), true)
  ```

* Plot the fields you want using PL_FIELD
  
  ```julia
  PL_FIELD("Displacement", utot, conn_tris, conn_quads, coord, defCoord(coord, ux, uy, FACTOR), def) 
  PL_FIELD("Sigma X", avrsx, conn_tris, conn_quads, coord, defCoord(coord, ux, uy, FACTOR), def)  
  PL_FIELD("Sigma Y", avrsy, conn_tris, conn_quads, coord, defCoord(coord, ux, uy, FACTOR), def)  
  PL_FIELD("Tau XY", avrtxy, conn_tris, conn_quads, coord, defCoord(coord, ux, uy, FACTOR), def)  
  PL_FIELD("Von Mises", vm, conn_tris, conn_quads, coord, defCoord(coord, ux, uy, FACTOR), def)
  ```
  
  Deformed coords can be calculated using defCoord
  
  ```julia
  defCoord(coord, ux, uy, FACTOR)
  ```

---

---

See "quarter_of_plate.jl" example to use the program.