# To Do
## Important step to make the code available for Git:
* ~~Learn how to write a module in Julia...~~
* ~~Create modules:~~
    * ~~**InputManager**: read Input file and return conn and coord~~
    * ~~**FemElements**: is the library of All FEM Elements (insert TL1 and QL1 elements as default)~~
    * ~~**MyPlots**: put here all the function necessaire to plot geometries and results~~
        * Plot stressess
    * ~~**Solvers**: is the library of all solvers (Thermal and Elastic as default)~~
    * ~~**BCs**: put here stuff for BCs calculation (Neumann and Dirichlet)~~
    * ~~**MatLib**: Material library~~
    * ~~**ProblemType**: put here problem types (Plane Stress, Plane Strain, etc.)~~
    * ~~**StiffenssElements**: put here all stiffness calculations (not sure it is necessaire)~~
    * ~~**PostPorcessing**: put here post processing stuff~~
        * Calculate stressess
* Generalise the code (not only QL1 but also TL1). This task is due to have a more usable code not for a specific purpose (right now some routine should be generalized)
* Define an interactive way to select edges and apply BCs

## Less important task:
* Maybe B-Matrix can be created in an other function
* Create a Thermo-Mechanical Solver
* Modal analysis may be interesting!
* Mortar methods for non conformal meshes?
* Write Documentation
* Cite stuff

## Ideas:
* Implement NEFEM and non isoparametric elements
* Non-Conformal mapping?
* Contact Solver?
* Inverse-FEM is possible? Start from result to get informations upon contraint...
