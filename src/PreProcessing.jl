module PreProcessing

__precompile__(true)

export readInput, readMat, read_LS_PrePost, separate_conn, geoCheckTris, geoCheckQuads

using DelimitedFiles: readdlm, findall
using CairoMakie, LinearAlgebra
using GLMakie.GeometryBasics

function read_LS_PrePost(NAME)

    println("**READING INPUT FILE .k from LS-PrePost...**")
    f = readdlm("geometries//$NAME")
    
    nrow = size(f,1)
    
    for i in 1:nrow
        # Get connectivity raw matrix
        if f[i,1] == "*ELEMENT_SHELL"
            count = 2
            global raw_conn = Array{Int64}(undef, 0, 10)
            while count != 0
                line = f[i+count,:]
                raw_conn = vcat(raw_conn, line[1:end-1]')
                count += 1
                if f[i+count,1] == "*NODE"
                    count = 0
                end
            end
        end
    
        # Get coordinates raw matrix
        if f[i,1] == "*NODE"
            count = 2
            global raw_coord = Array{Float64}(undef, 0, 6)
            while count != 0
                line = f[i+count,:]
                raw_coord = vcat(raw_coord, line[1:6]')
                count += 1
                if f[i+count,1] == "*END"
                    count = 0
                end
            end
        end
    end
    
    nel = size(raw_conn,1)
    nnode = size(raw_coord,1)
    
    conn = zeros(Int64,nel,7)
    
    tris_count = 0
    quads_count = 0
    
    for i in 1:nel
        conn[i,1] = raw_conn[i,1]
        conn[i,3:end] = raw_conn[i,2:6]
        if raw_conn[i,5] == raw_conn[i,6]
            conn[i,2] = 3 # it is a triangle
            tris_count += 1
        else
            conn[i,2] = 4 # it is a quadrangle
            quads_count += 1
        end
    end
    
    coord = raw_coord[:,1:3]
    coord = convert(Array{Float64}, coord)

    println("Number of TL elements: ", tris_count)
    println("Number of QL elements: ", quads_count)
    println("Total number of elements: ", nel)
    println("Total number of nodes: ", nnode)
    
    return nel, nnode, conn, coord
    
    end

function separate_conn(conn)
    conn_tris = Array{Int64}(undef, 0, 6)
    conn_quads = Array{Int64}(undef, 0, 7)
    
    nel = size(conn,1)
    for i in 1:nel
        if conn[i,2] == 3
            conn_tris = vcat(conn_tris, conn[i,1:6]')
        elseif conn[i,2] == 4
            conn_quads = vcat(conn_quads, conn[i,1:7]')
        end
    end
    
    return conn_tris, conn_quads
end

function readInput(NAME)
    println("**READING INPUT FILE...**")
    f = readdlm("geometries//$NAME.txt")
    nel = f[1, 1]
    g = f[2:end, :]
    # Create Connectivity for TL elements
    nel_TL = lastindex(findall(x -> x == 3, f[2:nel+1, 2]))
    if nel_TL > 0
        conn_tris = g[findall(x -> x == 3, f[2:nel+1, 2]), 1:6]
        conn_tris = convert(Array{Int64}, conn_tris)
    else
        conn_tris = Array{Float64}(undef,0)
    end

    # Create Connectivity for QL elements
    nel_QL = lastindex(findall(x -> x == 4, f[2:nel+1, 2]))
    if nel_QL > 0
        conn_quads = g[findall(x -> x == 4, f[2:nel+1, 2]), 1:7]
        conn_quads = convert(Array{Int64}, conn_quads)
    else
        conn_quads = Array{Float64}(undef,0)
    end

    println("Number of TL elements: ", nel_TL)
    println("Number of QL elements: ", nel_QL)
    println("Total number of elements: ", nel)

    # Create Coordinates vector (coord)
    coord = f[nel+2:end, 1:3]
    coord = convert(Array{Float64}, coord)
    nnode = size(coord, 1)
    println("Total number of nodes: ", nnode)

    return conn_tris, conn_quads, coord
end

function readMat(NAME)
    f = readdlm("materials//$NAME.txt")
    E = f[2,1]
    nu = f[2,2]
    thick = f[2,3]

    return E, nu, thick
end

function triangleArea(conn, coord)

    # Calculate length of first edge 1-2:
    x1 = coord[conn[1,4],2]
    x2 = coord[conn[1,5],2]
    y1 = coord[conn[1,4],3]
    y2 = coord[conn[1,5],3]
    a = sqrt((x1-x2)^2 + (y1-y2)^2)

    # Calculate length of second edge 2-3:
    x1 = coord[conn[1,5],2]
    x2 = coord[conn[1,6],2]
    y1 = coord[conn[1,5],3]
    y2 = coord[conn[1,6],3]
    b = sqrt((x1-x2)^2 + (y1-y2)^2)

    # Calculate length of third edge 3-1:
    x1 = coord[conn[1,6],2]
    x2 = coord[conn[1,4],2]
    y1 = coord[conn[1,6],3]
    y2 = coord[conn[1,4],3]
    c = sqrt((x1-x2)^2 + (y1-y2)^2)

    # Calculate half perimeter
    p = (a+b+c)/2

    # Calculate Triangle Area using Erone's Formula:
    A = sqrt(p*(p-a)*(p-b)*(p-c))

    return A
end

function quadrangleArea(conn, coord)

    qA = 0.0

    for i in 1:2

        if i == 1
            one = 4
            two = 5
            three = 6
        else
            one = 4
            two = 6
            three = 7
        end

        # Calculate length of first edge 1-2:
        x1 = coord[conn[1,one],2]
        x2 = coord[conn[1,two],2]
        y1 = coord[conn[1,one],3]
        y2 = coord[conn[1,two],3]
        a = sqrt((x1-x2)^2 + (y1-y2)^2)

        # Calculate length of second edge 2-3:
        x1 = coord[conn[1,two],2]
        x2 = coord[conn[1,three],2]
        y1 = coord[conn[1,two],3]
        y2 = coord[conn[1,three],3]
        b = sqrt((x1-x2)^2 + (y1-y2)^2)

        # Calculate length of third edge 3-1:
        x1 = coord[conn[1,three],2]
        x2 = coord[conn[1,one],2]
        y1 = coord[conn[1,three],3]
        y2 = coord[conn[1,one],3]
        c = sqrt((x1-x2)^2 + (y1-y2)^2)

    
        # Calculate half perimeter
        p = (a+b+c)/2

        # Calculate Triangle Area using Erone's Formula:
        tA = sqrt(p*(p-a)*(p-b)*(p-c))
        # println("tA: ",tA)
        qA = qA + tA

    end

    return qA
end

function distortionQL1(conn, coord)

    A = 4 * quadrangleArea(conn,coord)
    conn123 = [1 1 1 conn[1,4] conn[1,5] conn[1,6]]
    conn124 = [1 1 1 conn[1,4] conn[1,5] conn[1,7]]
    B = 4 * (triangleArea(conn123, coord) - 
             triangleArea(conn124, coord))

    conn134 = [1 1 1 conn[1,4] conn[1,6] conn[1,7]]
    C = 4 * (triangleArea(conn134, coord) - 
             triangleArea(conn124, coord))

    delta = (sqrt(B^2+C^2))/A

    
    return delta
end

function geoCheckQuads(coord, conn)

    nel = size(conn,1)
    for i in 1:nel
        delta = distortionQL1(conn[i,:]', coord)

        one = conn[i,4]
        two = conn[i,5]
        three = conn[i,6]
        four = conn[i,7]

        x1 = coord[one, 2]
        x2 = coord[two, 2]
        x3 = coord[three, 2]
        x4 = coord[four, 2]

        y1 = coord[one, 3]
        y2 = coord[two, 3]
        y3 = coord[three, 3]
        y4 = coord[four, 3]

        p = Polygon(Point2f[(x1, y1), (x2, y2), (x3, y3), (x4, y4)])
        xx = vec([x1 x2 x3 x4 x1]) 
        yy = vec([y1 y2 y3 y4 y1])

        if delta >= 0.2
            println("WARNING: element ",i, " appear to be distorted!")
            poly!(p, color = :red)
            lines!(xx,yy,linestyle=:solid, color = :black)
        else
            poly!(p, color = :green)
            lines!(xx,yy,linestyle=:solid, color = :black)
        end
    end

end

function geoCheckTris(coord, conn)

    nel = size(conn,1)
    for i in 1:nel
        x1 = coord[conn[i,4],2]
        x2 = coord[conn[i,5],2]
        x3 = coord[conn[i,6],2]

        y1 = coord[conn[i,4],3]
        y2 = coord[conn[i,5],3]
        y3 = coord[conn[i,6],3]

        b1 = y3-y2
        b2 = y1-y3
        b3 = y1-y2

        c1 = x2-x3
        c2 = x1-x3
        c3 = x1-x2

        # Calculate Jacobian
        J = [-c3 -b3
            -c2 -b2]

        # Calculate Area
        A = det(J)/2

        # Matrix B
        Bs = [-1 1 0
              -1 0 1]

        # Matrix B in cartesian system
        Bi = J \ Bs 

        # Stiffness Matrix
        K = A*Transpose(Bi)*Bi

        # Calculate eigevalues
        eigv = eigvals(K)

        # Remove 0 eigevalues
        eps = 1e-10
        eigv2 = eigv[findall(x->x>eps, eigv)]

        # Conditioning factor normalized
        cond = (maximum(eigv2)/minimum(eigv2))/3

        p = Polygon(Point2f[(x1, y1), (x2, y2), (x3, y3)])
        xx = vec([x1 x2 x3 x1]) 
        yy = vec([y1 y2 y3 y1])

        if cond >= 2
            println("WARNING: element ",i, " appear to be distorted!")
            poly!(p, color = :red)
            lines!(xx,yy,linestyle=:solid, color = :black)
        else
            poly!(p, color = :green)
            lines!(xx,yy,linestyle=:solid, color = :black)
        end

    end
end

end