module PostProcessing

__precompile__(true)

export splitU, defCoord, uTot, stressCalcGP_tris, stressCalcGP_quads, avarageStress, vmStress, nodalStresses

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
        uu = sqrt(ux[i]^2+uy[i]^2)
        append!(utot, uu)
    end
    return utot
end

function stressCalcGP_quads(U, conn, coord, MAT_NAME, PROBLEM_TYPE)
    # Iterate over all elements
    nel = size(conn,1)

    ngauss = 4

    mat = material(MAT_NAME)
    if PROBLEM_TYPE == "Plane_Stress"
        D = Plane_Stress(mat.E, mat.nu, mat.thick)
    elseif PROBLEM_TYPE == "Plane_Strain"
        D = Plane_Strain(mat.E, mat.nu)
    elseif PROBLEM_TYPE == "Axisymmetric"

    end

    sigmaXel = Array{Float64}(undef, 0, ngauss+1)
    sigmaYel = Array{Float64}(undef, 0, ngauss+1)
    tauXYel = Array{Float64}(undef, 0, ngauss+1)

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
 
        append!(sigmaX, conn[i,1])
        append!(sigmaY, conn[i,1])
        append!(tauXY, conn[i,1])

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

function stressCalcGP_tris(U, conn, coord, MAT_NAME, PROBLEM_TYPE)
    # Iterate over all elements
    nel = size(conn,1)

    ngauss = 1

    mat = material(MAT_NAME)
    if PROBLEM_TYPE == "Plane_Stress"
        D = Plane_Stress(mat.E, mat.nu, mat.thick)
    elseif PROBLEM_TYPE == "Plane_Strain"
        D = Plane_Strain(mat.E, mat.nu)
    elseif PROBLEM_TYPE == "Axisymmetric"

    end

    sigmaXel = Array{Float64}(undef, 0, ngauss+1)
    sigmaYel = Array{Float64}(undef, 0, ngauss+1)
    tauXYel = Array{Float64}(undef, 0, ngauss+1)

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

        J::Matrix = zeros(2,2)
        BB::Matrix = zeros(2,4)
        B::Matrix = zeros(3,6)

        # Gauss Points Calculation:
        pg = [1/3]
        w = [0.5]

        append!(sigmaX, conn[i,1])
        append!(sigmaY, conn[i,1])
        append!(tauXY, conn[i,1])
        
        # BE CAREFUL TO EXTEND THIS FUNCTION TO p>1 BECAUSE THE INTEGRATION DOMAIN ISN'T A SQUARE!!
        for i in 1:lastindex(pg) 
            xi = pg[i]
            eta = pg[i]
            ww = w[i]

            # Calculate derivative of shape fucntion over (xi,eta)
            dNxi, dNeta = FemElements.TL1_dN(xi,eta)
    
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
                 0       BB[2,3] BB[1,3] ]
            
            sigmaNode = D*B'*Delta
            append!(sigmaX, sigmaNode[1])
            append!(sigmaY, sigmaNode[2])
            append!(tauXY, sigmaNode[3])
        end

        sigmaXel = vcat(sigmaXel, sigmaX')
        sigmaYel = vcat(sigmaYel, sigmaY')
        tauXYel = vcat(tauXYel, tauXY')
    end
return sigmaXel, sigmaYel, tauXYel
end

function nodalStresses(stress, conn, coord)
    nnode = size(coord,1)
    nel = size(conn,1)
    if nel > 0
        type_el = lastindex(conn[1,4:end])
    else
        type_el = 0
        ss = Array{Float64}(undef, 0, type_el+1)
    end
    
    if type_el == 3
        ss = Array{Float64}(undef, 0, type_el+1)
        for i in 1:nel
            tmp = vec([stress[i,1] stress[i,2] stress[i,2] stress[i,2]])
            ss = vcat(ss, tmp')
        end
    elseif type_el == 4
        ss = Array{Float64}(undef, 0, type_el+1)
        for i in 1:nel
            nodes = conn[i,4:end]
            # get back stresses to nodes form GPs
            T = [1+sqrt(3)/2 -1/2 -1/2 1-sqrt(3)/2 ;
                -1/2 1-sqrt(3)/2 1+sqrt(3)/2 -1/2 ;
                1-sqrt(3)/2 -1/2 -1/2 1+sqrt(3)/2 ;
                -1/2 1+sqrt(3)/2 1-sqrt(3)/2 -1/2 ]
            snode = T*stress[i,2:end]
            tmp = append!([stress[i,1]], snode)
            ss = vcat(ss, tmp')
        end
    end
    return ss
end

function avarageStress(stress_tris, stress_quads, conn_tris, conn_quads, coord)

    nnode = size(coord,1)
    ss = [Vector{Float64}(undef,0) for _ in 1:nnode]

    # Trialgles:
    nel_tris = size(stress_tris,1)
    for i in 1:nel_tris
        nodes = conn_tris[i,4:end]
        for j in 1:lastindex(nodes)
            append!(ss[nodes[j]],stress_tris[i,j+1])
        end
    end

    # Quadrangles:
    nel_quads = size(stress_quads,1)
    for i in 1:nel_quads
        nodes = conn_quads[i,4:end]
        for j in 1:lastindex(nodes)
            append!(ss[nodes[j]],stress_quads[i,j+1])
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