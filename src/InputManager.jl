module InputManager

__precompile__(true)

export readInput, readMat

using DelimitedFiles: readdlm, findall

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