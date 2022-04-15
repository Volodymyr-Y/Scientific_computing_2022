using PyPlot
using PyCall
@pyimport matplotlib.animation as animation

function plot_results(grid, solution, nx, filename; savepdf = false, projection = "2d")

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
       
        savefig("T"*filename*".pdf")

        clf()
        fig_2 = figure()

        ax = fig_2.add_subplot(1, 2, 2, projection="3d")

        plot_surface(yy, xx, P, rstride=1, cstride=1, cmap=PyPlot.cm.coolwarm,
                            linewidth=0, antialiased=false)

        ylabel("x-dimension (-)")
        xlabel("y-dimension (-)")
        zlabel("pressure (-)")

        savefig("P"*filename*".pdf")

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

function plot_flux(xx, yy,fluxY, fluxX; scaler = true, save = true)
    clf()
    fig, ax = subplots()
    ax.set_title("water mass flux (q)")
  
    if scaler
      scaler_mult =  vcat(LinRange(1, 2, 40)) .* ones(40, 60)
      Q = ax.quiver(xx[1:25:end, 1:5:end], yy[1:25:end, 1:5:end], fluxX[1:25:end, 1:5:end], (scaler_mult.*fluxY[1:25:end, 1:5:end])./2.0, (scaler_mult.*fluxY[1:25:end, 1:5:end])./2.0, units = "height", cmap = cmap=PyPlot.cm.coolwarm)
  
    else 
      Q = ax.quiver(xx[1:25:end, 1:5:end], yy[1:25:end, 1:5:end], fluxX[1:25:end, 1:5:end], fluxY[1:25:end, 1:5:end], fluxY[1:25:end, 1:5:end], units = "height", cmap = cmap=PyPlot.cm.coolwarm)
    end
  
    # qk = ax.quiverkey(Q, 0.9, 0.9, 2, "water mass flux", labelpos="E",
    #                    coordinates="figure")
  
  
  
    cm = PyPlot.cm.coolwarm
  
    sm = PyPlot.cm.ScalarMappable(cmap=cm)
    sm.set_array([minimum(fluxY), maximum(fluxY)])
  
    colorbar(sm)
  
    # ax.set_aspect("equal")
  
    xlabel("x-dimension (-)")
    ylabel("y-dimension (-)")
  
    if save
      savefig("flux.pdf")
    end
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

