using PyPlot
using PyCall
@pyimport matplotlib.animation as animation

function plot_results(grid, solution, nx; savepdf = false, projection = "2d")

    function meshgrid(xin, yin)
        nx  =   length(xin)
        ny  =   length(yin)
        xout=   zeros(ny,nx)
        yout=   zeros(ny,nx)
        for jx = 1:nx
            for ix = 1:ny
                xout[ix,jx] = xin[jx]
                yout[ix,jx] = yin[ix]
            end
        end
        return (x=xout, y=yout)
    end

    n  = length(grid)
    x  = LinRange(0.0, 300.0, nx)

    xx, yy = meshgrid(x, grid)
    T      = zeros(size(xx))
    P      = zeros(size(xx)) 
    

    if typeof(solution) == Vector{Float64}
        for i in 1:nx
            T[:, i]  = solution[1:n]
            P[:, i]  = solution[n+1:end]
        end
    end
    
    if projection == "3d"

        clf()
        fig_1 = figure()

        ax = fig_1.add_subplot(projection="3d")

        plot_surface(yy, xx, T, rstride=1, cstride=1, cmap=PyPlot.cm.coolwarm,
                            linewidth=0, antialiased=false)

        ylabel("x-dimension (-)")
        xlabel("y-dimension (-)")
        zlabel("temperature (-)")
       
        savefig("temperature.pdf")

        clf()
        fig_2 = figure()

        ax = fig_2.add_subplot(1, 2, 2, projection="3d")

        plot_surface(yy, xx, P, rstride=1, cstride=1, cmap=PyPlot.cm.coolwarm,
                            linewidth=0, antialiased=false)

        ylabel("x-dimension (-)")
        xlabel("y-dimension (-)")
        zlabel("pressure (-)")

        savefig("pressure.pdf")

    elseif projection == "2d"
        
        clf()
        fig_1 = figure()

        temp = contourf(xx, yy, T, cmap=PyPlot.cm.coolwarm, antialiased=false)

        ylabel("x-dimension (-)")
        xlabel("y-dimension (-)")

        fig_1.colorbar(temp)
       
        if savepdf
            savefig("temperature.pdf")
        end

        clf()
        fig_2 = figure()

        pres = contourf(xx, yy, P, cmap=PyPlot.cm.coolwarm, antialiased=false)

        ylabel("x-dimension (-)")
        xlabel("y-dimension (-)")
        
        fig_2.colorbar(pres)
        
        if savepdf
            savefig("pressure.pdf")
        end

    end
    return fig_1, fig_2
end


function animate_solution(grid, solution, nx; specie = "temperature", projection = "2d")
    
    function meshgrid(xin, yin)
        nx  =   length(xin)
        ny  =   length(yin)
        xout=   zeros(ny,nx)
        yout=   zeros(ny,nx)
        for jx = 1:nx
            for ix = 1:ny
                xout[ix,jx] = xin[jx]
                yout[ix,jx] = yin[ix]
            end
        end
        return (x=xout, y=yout)
    end

    n  = length(grid)
    x  = LinRange(0.0, 300.0, nx)

    xx, yy = meshgrid(x, grid)
    T      = zeros(size(xx))
    P      = zeros(size(xx)) 
    
    fig = figure()

    function animation_frame(i)
        clf()
        
        for idx in 1:nx
            T[:, idx]  = solution[i][1:n]
            P[:, idx]  = solution[i][n+1:end]
        end

        ax = fig.add_subplot(projection="3d")

        plot_surface(yy, xx, T, rstride=1, cstride=1, cmap=PyPlot.cm.coolwarm,
                            linewidth=0, antialiased=false)

        title("t = " + str(solution.t[i]) + " s")
        ylabel("x-dimension (-)")
        xlabel("y-dimension (-)")
        zlabel("temperature (-)")
    end

    gif = animation.FuncAnimation(fig, animation_frame, frames=100, interval=10, blit= false)
    gif[:save]("animation.gif", bitrate=-1, extra_args=["-vcodec", "libx264", "-pix_fmt", "yuv420p"])
end