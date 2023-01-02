module MyPlots

__precompile__(true)

include("FemElements.jl")
import .FemElements

export plotMesh, plotBC, plotField

using CairoMakie, LinearAlgebra

    # Plot Mesh:
    function plotMesh(conn, coord, linetype, linecolor)
        # Read connectivity and coords and plot elements
        for i in 1:size(conn,1)
            x = coord[conn[i,4:end],2]
            y = coord[conn[i,4:end],3]
            # Add last node
            xx = [x ; x[1]] 
            yy = [y ; y[1]]
            # Display
            lines!(xx,yy,linestyle=linetype,color = linecolor)
            scatter!(x,y, color = :black, markersize = 15px)
        end
    end

    # Plot BC Nodes:
    function plotBC(x, y, col)
        scatter!(x, y, color=col, markersize=12px)
    end

    function plotField(coord, conn, field)

        nel = size(conn,1)

        # PLOT CONTOUR
        xs = LinRange(-1, 1, 4)
        ys = LinRange(-1, 1, 4)
        
        # Loop over elements
        for i in 1:nel
            pts = coord[conn[i,4:end], 2:end]
            x = pts[:, 1]
            y = pts[:, 2]
        
            # Recover nodal values of the field (Displacement/Stresses/etc.)
            field_node = field[conn[i,4:end]]
        
            xg = Array{Float64}(undef,0)
            yg = Array{Float64}(undef,0)
            zg = Array{Float64}(undef,0)
        
            for i in 1:lastindex(xs)
                for j in 1:lastindex(ys)
                    N = FemElements.QL1_N(xs[i], ys[j])
                    append!(xg, N'*x)
                    append!(yg, N'*y)
                    VALUE = dot(N,field_node)
                    append!(zg, VALUE)
                end
            end
            min = minimum(field)
            max = maximum(field)
            lev = 10
            avr = abs(max-min)
            tricontourf!(xg, yg, zg, extendlow = :auto, extendhigh = :auto, mode = :normal, levels = min:avr/lev:max)
        end
    end

end