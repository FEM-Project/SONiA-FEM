module MyPlots

__precompile__(true)

include("FemElements.jl")
include("PostProcessing.jl")
import .FemElements, .PostProcessing

export plotMesh, plotBC, plotField, PL, PL_FIELD, plotDistributed

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
            # scatter!(x,y, color = :black, markersize = 15px)
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
            if avr < 1e-8
                avr = 1
            end
            
            tricontourf!(xg, yg, zg, extendlow = :auto, extendhigh = :auto, 
            mode = :normal, levels = min:avr/lev:max)
        end
    end

    # Try to use multiple dispatch
    function PL(TITLE, conn_tris, conn_quads, coord, coord_def, def)
        p = Figure()
        ax = Axis(p[1, 1], aspect=DataAspect(),title=TITLE)
        if def
            plotMesh(conn_tris, coord, :dash, :gray)
            plotMesh(conn_quads, coord, :dash, :gray)
            plotMesh(conn_tris, coord_def, :solid, :black)
            plotMesh(conn_quads, coord_def, :solid, :black)
        else
            plotMesh(conn_tris, coord, :solid, :black)
            plotMesh(conn_quads, coord, :solid, :black)
        end
        GLMakie.activate!()
        display(GLMakie.Screen(), p)
        return p
    end

    function PL_FIELD(TITLE, FIELD, conn_tris, conn_quads, coord, coord_def, def)
        p = Figure()
        ax = Axis(p[1, 1], aspect=DataAspect(),title = TITLE)
        if def
            plotMesh(conn_tris, coord, :dash, :gray)
            plotMesh(conn_quads, coord, :dash, :gray)

            plotField(coord_def, conn_tris, FIELD)
            plotField(coord_def, conn_quads, FIELD)

            plotMesh(conn_tris, coord_def, :solid, :black)
            plotMesh(conn_quads, coord_def, :solid, :black)
        else
            plotField(coord, conn_tris, FIELD)
            plotField(coord, conn_quads, FIELD)

            plotMesh(conn_tris, coord, :solid, :black)
            plotMesh(conn_quads, coord, :solid, :black)
        end

        min = minimum(FIELD)
        max = maximum(FIELD)
        eps = 1e-8
        if max-min <= eps && max-min != 0
            Colorbar(p[1, 2], label = TITLE, colormap = :viridis, limits = (min*0.999, min*1.001))
        elseif max-min <= eps && max-min == 0
            Colorbar(p[1, 2], label = TITLE, colormap = :viridis, limits = (-0.001, 0.001))
        else
            Colorbar(p[1, 2], label = TITLE, colormap = :viridis, limits = (min, max))
        end
        
        GLMakie.activate!()
        display(GLMakie.Screen(), p)
    end

    function plotDistributed(neumann_edges, coord, NORMAL_FORCE)
        
        ## ISSUES:
            # - the arrow length is proportional to the edge length
            #   and not to the force values
            # - compression is not moved out the geometry
        
        for i in 1:size(neumann_edges, 1)
            edge = trunc.(Int, neumann_edges[i, :])
            x = coord[edge, 2]
            y = coord[edge, 3]

            x1 = x[1]
            x2 = x[2]
            y1 = y[1]
            y2 = y[2]
            
            dx = x2-x1
            dy = y2-y1
            ll = sqrt(dx^2+dy^2)

            s = sign(NORMAL_FORCE)
            f = abs(NORMAL_FORCE)

            narrow = 3

            if dx != 0
                m = dy/dx
                alfa = atan(m)
                qq = (x2*y1 - x1*y2)/(x2-x1)
                xx = LinRange(x1,x2,narrow+2)
                yy = m.*xx.+qq
            else
                yy = LinRange(y1,y2,narrow+2)
                xx = fill(x1, lastindex(yy))
                alfa = pi/2
            end

            xx = xx[2:end-1]
            yy = yy[2:end-1]

            xs = [x1; x2]
            ys = [y1; y2]
            

            #= # primo quadrante:
            if alfa >= 0 && alfa <= pi/2
                ddx = s*f*cos(alfa)
                ddy = s*f*sin(alfa)
            elseif alfa > pi/2 && alfa < pi
                ddx = -s*f*cos(alfa)
                ddy = -s*f*sin(alfa)
            elseif alfa >= pi && alfa <= 3*pi/4
                ddx = s*f*cos(alfa)
                ddy = s*f*sin(alfa)
            else
                ddx = -s*f*cos(alfa)
                ddy = -s*f*sin(alfa)
            end =#

            if s < 0
                ux = fill(-dy, lastindex(xx))
                uy = fill(dx, lastindex(yy))
            elseif s > 0
                ux = fill(dy, lastindex(xx))
                uy = fill(-dx, lastindex(yy))
            end
            arrows!(xx, yy, ux, uy, arrowsize = 20.0, lengthscale = 0.5, color = :blue)
        end
    end

end