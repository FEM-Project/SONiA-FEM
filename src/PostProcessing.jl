module PostProcessing

__precompile__(true)

export splitU, defCoord, uTot, stressCalcGP, avarageStress, vmStress

include("legzo.jl")
include("FemElements.jl")
include("MatLib.jl")
include("ProblemType.jl")

import .FemElements
using LinearAlgebra, .MatLib, .ProblemType, Statistics

function splitU(U, coord)

    # Define Necessaire Quantities:
    n_dof = 2
    nnode = size(coord,1)
    ngdof = nnode * n_dof
    
    ux = zeros(trunc(Int,ngdof/2),1)
    uy = zeros(trunc(Int,ngdof/2),1)
    j = 1
    k = 1
    for i in 1:ngdof
        if i%2 == 0
            uy[j] = U[i]
            j += 1
        else
            ux[k] = U[i]
            k += 1
        end
    end
    return ux, uy
end

function defCoord(coord, ux, uy, f)
    nnode = size(coord,1)
    coord_def = zeros(nnode,3)
    coord_def[:,1] = coord[:,1]
    coord_def[:,2] = coord[:,2].+ux.*f
    coord_def[:,3] = coord[:,3].+uy.*f
    return coord_def
end

function uTot(ux, uy)
    utot = Array{Float64}(undef,0)
    for i in 1:lastindex(ux)
    append!(utot, norm(ux[i],uy[i]))
    end
    return utot
end

function stressCalcGP(U, conn, coord, MAT_NAME, PROBLEM_TYPE)
    # Iterate over all elements
    nel = size(conn,1)

    mat = material(MAT_NAME)
    if PROBLEM_TYPE == "Plane_Stress"
        D = Plane_Stress(mat.E, mat.nu, mat.thick)
    elseif PROBLEM_TYPE == "Plane_Strain"

    elseif PROBLEM_TYPE == "Axisymmetric"

    end

    sigmaXel = Array{Float64}(undef, 0, 4)
    sigmaYel = Array{Float64}(undef, 0, 4)
    tauXYel = Array{Float64}(undef, 0, 4)

    for i in 1:nel

        # Nodal displacements
        nodes = conn[i,4:end]
        
        sigmaX = Array{Float64}(undef, 0)
        sigmaY = Array{Float64}(undef, 0)
        tauXY = Array{Float64}(undef, 0)

        nodes_gdl = Array{Int64}(undef, 0)
        for j in 1:lastindex(nodes)
            append!(nodes_gdl, nodes[j]*2-1)
            append!(nodes_gdl, nodes[j]*2)
        end

        Delta = U[nodes_gdl]

        x = coord[nodes, 2]
        y = coord[nodes, 3]

        # Loop over GP
        n_gauss = 2

        J::Matrix = zeros(2,2)
        BB::Matrix = zeros(2,4)
        B::Matrix = zeros(3,8)


        # Gauss Points Calculation:
        pg,w = legzo(n_gauss,1,-1)
 
        for i in 1:lastindex(pg) 
            xi = pg[i]
        
            for j in 1:lastindex(pg)
                eta = pg[j]

                # Calculate derivative of shape fucntion over (xi,eta)
                dNxi, dNeta = FemElements.QL1_dN(xi,eta)
    
                # Calculate Jacobian Matrix components using dot() fucntion in order to return ascalar value
                J[1,1] = dot(dNxi,x)
                J[2,1] = dot(dNeta,x)
                J[1,2] = dot(dNxi,y)
                J[2,2] = dot(dNeta,y)
                # Calculate BB Matrix (Derivative of Shape Functions in cartesian space)
                BB = inv(J)*[dNxi' ; dNeta']
            
                # Assemble B Matrix (Structural Analysis)
                B = [BB[1,1] 0       BB[2,1] ; 
                     0       BB[2,1] BB[1,1] ;
                     BB[1,2] 0       BB[2,2] ;
                     0       BB[2,2] BB[1,2] ;
                     BB[1,3] 0       BB[2,3] ;
                     0       BB[2,3] BB[1,3] ;
                     BB[1,4] 0       BB[2,4] ;
                     0       BB[2,4] BB[1,4]]
                
                sigmaNode = D*B'*Delta
                append!(sigmaX, sigmaNode[1])
                append!(sigmaY, sigmaNode[2])
                append!(tauXY, sigmaNode[3])
            end
        end
        
        sigmaXel = vcat(sigmaXel, sigmaX')
        sigmaYel = vcat(sigmaYel, sigmaY')
        tauXYel = vcat(tauXYel, tauXY')
    end
return sigmaXel, sigmaYel, tauXYel
end

function avarageStress(stress, conn, coord)

    nnode = size(coord,1)
    nel = size(stress,1)
    ss = [Vector{Float64}(undef,0) for _ in 1:nnode]
        
    for i in 1:nel
        nodes = conn[i,4:end]
        # get back stresses to nodes form GPs
        T = [1+sqrt(3)/2 -1/2 -1/2 1-sqrt(3)/2 ;
            -1/2 1-sqrt(3)/2 1+sqrt(3)/2 -1/2 ;
            1-sqrt(3)/2 -1/2 -1/2 1+sqrt(3)/2 ;
            -1/2 1+sqrt(3)/2 1-sqrt(3)/2 -1/2 ]
        snode = T*stress[i,:]

        for j in 1:lastindex(nodes)
            append!(ss[nodes[j]],snode[j])
        end
    end

    avrStress = Array{Float64}(undef, 0)
    # create Avarage Vector
    for i in 1:nnode
        append!(avrStress, mean(ss[i]))
    end

    return avrStress
    
end

function vmStress(sx,sy,txy)
    # Pass to this function avarage values of stresses components
    vm = sqrt.(sx.^2-sx.*sy+sy.^2+3*txy.^2)
    return vm
end

end