module MyPlots

__precompile__(true)

include("FemElements.jl")
include("PostProcessing.jl")
import .FemElements, .PostProcessing

export plotMesh, plotBC, plotField, PL

using LinearAlgebra, GLMakie

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
        if nel > 0
            type_el = lastindex(conn[1,4:end])
        else
            type_el = 0
        end

        # PLOT CONTOUR
        if type_el == 3
            # Domain 0 - 1
            xs = LinRange(0, 1, 4)
            ys = LinRange(0, 1, 4)
        elseif type_el == 4
            # Domain -1 - 1
            xs = LinRange(-1, 1, 4)
            ys = LinRange(-1, 1, 4)
        end

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
        
            if type_el == 3
                # Countour in a triangular domain
                for i in 1:lastindex(xs)
                    for j in 1:lastindex(ys)-i+1
                        N = FemElements.TL1_N(xs[i], ys[j])
                        append!(xg, N'*x)
                        append!(yg, N'*y)
                        VALUE = dot(N,field_node)
                        append!(zg, VALUE)
                    end
                end
            elseif type_el == 4
                # Contour in a quadrangular domain
                for i in 1:lastindex(xs)
                    for j in 1:lastindex(ys)
                        N = FemElements.QL1_N(xs[i], ys[j])
                        append!(xg, N'*x)
                        append!(yg, N'*y)
                        VALUE = dot(N,field_node)
                        append!(zg, VALUE)
                    end
                end
            end

            min = minimum(field)
            max = maximum(field)
            lev = 10
            avr = abs(max-min)
            tricontourf!(xg, yg, zg, extendlow = :auto, extendhigh = :auto, mode = :normal, levels = min:avr/lev:max)
        end
    end

    # Try to use multiple dispatch
    function PL(TITLE, conn_tris, conn_quads, coord)
        p = Figure()
        ax = Axis(p[1, 1], aspect=DataAspect(),title=TITLE)
        plotMesh(conn_tris, coord, :dash, :gray)
        plotMesh(conn_quads, coord, :dash, :gray)
        GLMakie.activate!()
        display(GLMakie.Screen(), p)
    end

    function PL(TITLE, FIELD, conn_tris, conn_quads, coord)
        p = Figure()
        ax = Axis(p[1, 1], aspect=DataAspect(),title = TITLE)
        plotField(coord, conn_tris, FIELD)
        plotField(coord, conn_quads, FIELD)
        plotMesh(conn_tris, coord, :solid, :black)
        plotMesh(conn_quads, coord, :solid, :black)
        GLMakie.activate!()
        display(GLMakie.Screen(), p)
    end

end