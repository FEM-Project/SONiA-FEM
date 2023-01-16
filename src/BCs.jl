module BCs

__precompile__(true)

export BC_box, BC_Dirichlet, BC_Neumann, arch, line

include("legzo.jl")
include("FemElements.jl")
import .FemElements

using StaticArrays, PolygonOps, LinearAlgebra

# Utilities
function BC_box(xv, yv, coord)

    nnode = size(coord, 1)

    # Create Volume Box
    polygon = SVector.(xv, yv)

    # Defining Search points
    points = vec(SVector.(coord[:, 2], coord[:, 3]))

    # Find (true or false) the inside point
    inside = [inpolygon(p, polygon; in=true, on=false, out=false) for p in points]

    # Calculate the number of inside points
    n_BC = count(i -> (i == 1), inside)

    # Generate BC_nodeBox Vector
    BC_nodeBox = zeros(n_BC, 3)
    j = 1
    for i in 1:nnode
        if inside[i]
            BC_nodeBox[j, 1] = i
            BC_nodeBox[j, 2:3] = coord[i, 2:3]
            j = j + 1
        end
    end

    return BC_nodeBox
end

function arch(centre,r,alfa1,alfa2)
    # Angle growing clockwise
    x0 = centre[1]
    y0 = centre[2]
    step = 100
    t = LinRange(alfa1*pi/180,alfa2*pi/180,step)
    x = x0 .+ r*cos.(t)
    y = y0 .+ r*sin.(t)
    return x, y
end

function line(p1,p2)
    x = [p1[1] p2[1]]
    y = [p1[2] p2[2]]
    return x, y
end

# Dirichlet
function BC_Dirichlet(nodes, BC_type, ux, uy)
    # Number Dirichlet nodes per Box
    nnode = size(nodes, 1)

    # Define BC_Dirichlet vector
    Dirichlet = zeros(nnode, 4)
    Dirichlet[:, 1] = nodes[:, 1]
    Dirichlet[:, 2] = fill(BC_type, nnode)
    Dirichlet[:, 3] = fill(ux, nnode)
    Dirichlet[:, 4] = fill(uy, nnode)

    return Dirichlet
end

# Neumann
function BC_Neumann(coord, conn, NeumannNodes, applied_p, applied_t)

    nel = size(conn,1)
    nnode = size(coord,1)
    n_dof = 2
    ngdof = nnode * n_dof

    # Loop over elements
    neumann_edges = Array{Float64}(undef, 0, 2)
    for i in 1:nel
        n_edges = Array{Float64}(undef, 0)
        element_node = conn[i, 4:end]
        count = 0
        for j in 1:lastindex(element_node)
            if element_node[j] in NeumannNodes[:, 1]
                count += 1
                append!(n_edges, element_node[j])
            end
        end
        if count == 2
            #println("element ", i, " is on the edge")
            neumann_edges = vcat(neumann_edges, n_edges')
        end
    end
    println("Neumann edges are: ", neumann_edges)


    Fneumann = zeros(ngdof, 1)

    gp, w = legzo(1, -1, 1)

    for i in 1:size(neumann_edges, 1)

        edge = trunc.(Int, neumann_edges[i, :])
        x = coord[edge, 2]
        y = coord[edge, 3]

        Fneumann_tmp = zeros(ngdof, 1)

        pGxyz = 0
        pGxyz_n = 0
        pGxyz_t = 0
        ll = 0
        llx = 0
        lly = 0
        for j in 1:lastindex(gp)
            xi = gp[j]
            N = FemElements.EM1_N(xi)
            dN = FemElements.EM1_dN(xi)
            # Jacobian
            # NOTA: c'é un valore assoluto perché questa "lunghezza" non
            # deve essere negativa (non avrebbe senso), il segno deve
            # essere deciso dal carico applicato (???)
            Jx = dN' * x
            Jy = dN' * y
            J = norm([Jx Jy], 2)
            xx = N' * x
            yy = N' * y
            ll = ll + J * w[j]
            llx = llx + Jx * w[j]
            lly = lly + Jy * w[j]

            pn = Array{Float64}(undef, 0)
            pt = Array{Float64}(undef, 0)
            for k in 1:lastindex(x)
                # if you like you can define a function to apply
                append!(pn, applied_p)
                append!(pt, applied_t)
            end

            # Load interpolated with shape function for each edge
            pq_n = N' * pn
            pq_t = N' * pt

            # Normal calculation, ref [The finite element method. Vol 1 - Zienkiewicz,
            # Taylor] pag. 211
            # in this case [0 0 1] means positive traction and negative
            # compression, otherwise [0 0 -1] means the opposit
            pGxyz_n = pGxyz_n .+ pq_n * N * w[j] .* cross([Jx; Jy; 0], [0; 0; 1])'
            # Tangential pressure positive in anticlockwise direction
            pGxyz_t = pGxyz_t .+ pq_t * N * w[j] .* [Jx Jy 0]

            pGxyz = pGxyz_n + pGxyz_t

        end
        # llxy = sqrt(llx^2 + lly^2)
        # println("This edge length is: ", llxy)
        Fneumann_tmp[2*edge.-1] = pGxyz[:, 1]
        Fneumann_tmp[2*edge] = pGxyz[:, 2]
        Fneumann = Fneumann + Fneumann_tmp
    end
    return Fneumann
end

end