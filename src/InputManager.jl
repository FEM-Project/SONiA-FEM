module InputManager

__precompile__(true)

export readInput, readMat, read_LS_PrePost, separate_conn

using DelimitedFiles: readdlm, findall

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
end