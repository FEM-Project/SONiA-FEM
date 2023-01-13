module Solvers

__precompile__(true)

export elastSolver

include("StiffnessElements.jl")
include("MatLib.jl")
include("ProblemType.jl")

using SparseArrays, .StiffnessElements, .MatLib, .ProblemType


function elastSolver(conn_tris, conn_quads, coord, MAT_NAME, PROBLEM_TYPE, all_BC_Dirichlet, FF)

    # Define Necessaire Quantities:
    n_dof = 2
    nnode = size(coord,1)
    nel_TL = size(conn_tris,1)
    nel_QL = size(conn_quads,1)
    ngdof = nnode * n_dof

    println("\n **ASSEMBLE STIFFNESS MATRIX...**")
    # Initialise Ktot
    ktot = spzeros(size(coord, 1) * n_dof, size(coord, 1) * n_dof)

    mat = material(MAT_NAME)
    if PROBLEM_TYPE == "Plane_Stress"
        D = Plane_Stress(mat.E, mat.nu, mat.thick)
    elseif PROBLEM_TYPE == "Plane_Strain"
        D = Plane_Strain(mat.E, mat.nu)
    elseif PROBLEM_TYPE == "Axisymmetric"

    end


    println("...assemble triangle elements")
    # Calculate Area and Stiffness Matrix for TL elements
    for i in 1:nel_TL
        pts = coord[conn_tris[i, 4:end], 2:end]
        nodes = conn_tris[i, 4:end]

        kel = stiff_TL1(pts, D, mat.thick)

        # ASSEMBLE STIFFNESS MATRIX
        # Defining assemble key
        key = Array{Int32}(undef, 0)
        for j in 1:lastindex(nodes)
            append!(key, [nodes[j] * 2 - 1 nodes[j] * 2])
        end
        # Assemble matrix
        for j in 1:lastindex(key)
            for k in 1:lastindex(key)
                ktot[key[j], key[k]] = ktot[key[j], key[k]] + kel[j, k]
            end
        end
    end

    println("...assemble quadrangle elements")
    # Calculate Area and Stiffness Matrix for QL elements
    for i in 1:nel_QL
        pts = coord[conn_quads[i, 4:end], 2:end]
        nodes = conn_quads[i, 4:end]

        kel = stiff_QL1(pts, D, mat.thick)

        # ASSEMBLE STIFFNESS MATRIX
        # Defining assemble key
        key = Array{Int32}(undef, 0)
        for j in 1:lastindex(nodes)
            append!(key, [nodes[j] * 2 - 1 nodes[j] * 2])
        end
        # Assemble matrix
        for j in 1:lastindex(key)
            for k in 1:lastindex(key)
                ktot[key[j], key[k]] = ktot[key[j], key[k]] + kel[j, k]
            end
        end
    end

    println("\n **APPLYING BOUNDARY CONDITIONS...**")
    
    # Reduce Stiffness Matrix
    
    i_remove = Array{Int32}(undef,0)
    valuedoffix = zeros(ngdof,1)
    
    for i in 1:size(all_BC_Dirichlet,1)
        if all_BC_Dirichlet[i,2] == 0.0
            # Do nothing
        elseif all_BC_Dirichlet[i,2] == 1.0
            # y constraint
            index = trunc(Int, all_BC_Dirichlet[i,1])
            append!(i_remove, index*2) 
            valuedoffix[index*2] = all_BC_Dirichlet[i,4]
        elseif all_BC_Dirichlet[i,2] == 10.0
            # x constraint
            index = trunc(Int, all_BC_Dirichlet[i,1])
            append!(i_remove, index*2-1) 
            valuedoffix[index*2-1] = all_BC_Dirichlet[i,3]
        elseif all_BC_Dirichlet[i,2] == 11.0
            # x-y constraint
            index = trunc(Int, all_BC_Dirichlet[i,1])
            append!(i_remove, [index*2-1 index*2]) 
            valuedoffix[index*2-1] = all_BC_Dirichlet[i,3]
            valuedoffix[index*2] = all_BC_Dirichlet[i,4]
        end
    end
    Kstar = ktot[1:end .∉ [i_remove], 1:end .∉ [i_remove]]
    
    ifdoffix = zeros(ngdof,1)
    for i in 1:ngdof
        if i in i_remove
            ifdoffix[i] = 1
        end
    end

    Fstar = FF[1:end .∉ [i_remove]]
    
    println("\n **SOLVE LINEAR SISTEM: LU DECOMPOSITION...**")
    Ustar = Kstar \ Fstar

    U = zeros(ngdof,1)

    #...assemble full global displacement vector
    icount = 0;
    for igdof = 1:ngdof
        if ifdoffix[igdof]==0.
            icount = icount+1;
            U[igdof]=Ustar[icount];
        else
            U[igdof]=valuedoffix[igdof];
        end
    end

return U
end

# Here you should put some other solvers!

end