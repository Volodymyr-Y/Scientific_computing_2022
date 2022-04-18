using SparseArrays
using ForwardDiff
using LinearAlgebra
using DifferentialEquations
using Sundials
using NLsolve
using BenchmarkTools
using TimerOutputs
using DelimitedFiles
using PyPlot
using PyCall

pygui(true)

function plot_results(gridx, gridy,T,P, filename; savepdf = false, projection = "2d")
    if projection == "3d"

        clf()
        fig_1 = figure()

        ax = fig_1.add_subplot(projection="3d")

        plot_surface(gridy, gridx, T, cmap=PyPlot.cm.coolwarm,
                            linewidth=0, antialiased=false, alpha=0.8)

        ylabel("x-dimension (-)")
        xlabel("y-dimension (-)")
        zlabel("temperature (-)")
       
        savefig("T"*filename*".pdf")

        clf()
        fig_2 = figure()

        ax = fig_2.add_subplot(1, 2, 2, projection="3d")

        plot_surface(gridy, gridx, P, rstride=1, cstride=1, cmap=PyPlot.cm.coolwarm,
                            linewidth=0, antialiased=false, alpha=0.8)

        ylabel("x-dimension (-)")
        xlabel("y-dimension (-)")
        zlabel("pressure (-)")

        savefig("P"*filename*".pdf")

    elseif projection == "2d"
        
        clf()
        fig_1 = figure()

        temp = contourf(gridx, gridy, T, cmap=PyPlot.cm.coolwarm, antialiased=false)

        ylabel("x-dimension (-)")
        xlabel("y-dimension (-)")

        fig_1.colorbar(temp)
       
        if savepdf
            savefig("temperature.pdf")
        end

        #clf()
        fig_2 = figure()

        pres = contourf(gridx, gridy, P, cmap=PyPlot.cm.coolwarm, antialiased=false)

        ylabel("x-dimension (-)")
        xlabel("y-dimension (-)")
        
        fig_2.colorbar(pres)
        
        if savepdf
            savefig("pressure.pdf")
        end

    end
    return fig_1, fig_2
end

function new_assemble_nonlinear_system(gridX,gridY,bc_bott,bc_top;heating_from_left = false)
    N = (n_y,n_x)
    λ = 0.01
    c = 0.001
    α = 0.01
    k = 100.0
    ρ_ref = 1.0
    ϵ = 10.0^(-6)
    h_x = gridX[:,2:N[2]] - gridX[:,1:N[2]-1]     
    h_y = gridY[1:(N[1]-1),:] - gridY[2:N[1],:]  
    h_left = vcat(h_y[1,:]',h_y) 
    h_right = vcat(h_y,h_y[end,:]')  
    h_bottom = hcat(h_x[:,1], h_x)
    h_top = hcat(h_x,h_x[:, end])  

    function ∇_left(u)
        u1 = zeros(typeof(u[1]),N)
        for j in 2:N[2]
            u1[:,j] = u[:,j-1]-u[:,j]
        end    
        return u1 .* h_left ./ h_bottom
    end

    function ∇_right(u)
        u1 = zeros(typeof(u[1]),N)
        for j in 1:N[2]-1
            u1[:,j] = u[:,j+1]-u[:,j] 
        end
        return u1 .* h_right ./ h_top
    end

    function ∇_bottom(u)
        u1 = zeros(typeof(u[1]),N)
        for i in 2:N[1]
            u1[i,:] = u[i-1,:]-u[i,:] 
        end
        return u1 .* h_bottom ./ h_left
    end

    function ∇_top(u)
        u1 = zeros(typeof(u[1]),N)
        for i in 1:N[1]-1
            u1[i,:] = u[i+1,:]-u[i,:] 
        end
        return u1 .* h_top ./ h_right
    end

    function L(u)
        u1 = zeros(typeof(u[1]),N)
        u1[:,1] = u[:,1]*2.0
        for i in 2:N[1]
            u1[i,:] = u[i-1,:]+u[i,:] 
        end
        return u1 / 2.0 
    end

    function R(u)
        u1 = zeros(typeof(u[1]),N)
        u1[:,end] = u[:,end]*2.0
        for i in 1:N[1]-1
            u1[i,:] = u[i+1,:]+u[i,:] 
        end
        return u1 / 2.0 
    end

    function B(u)
        u1 = zeros(typeof(u[1]),N)
        u1[end,:] = u[end,:]*2.0
        for i in 1:N[1]-1
                u1[i,:] = u[i+1,:] + u[i,:]
        end
        return u1 / 2.0
    end
    
    function North(u)
        u1 = zeros(typeof(u[1]),N)
        u1[:,1] = u[:,1]*2.0
        for i in 2:N[1]
            u1[i,:] = u[i-1,:] + u[i,:]
        end
        return u1 / 2.0
    end

    function system(x)
        T = x[1:N[1]*N[2]]
        P = x[N[1]*N[2]+1:2*N[1]*N[2]]
        T = reshape(T,(N[1],N[2]))
        P = reshape(P,(N[1],N[2]))
        
        F1 = (λ/c)*(∇_left(T)+∇_right(T)+∇_top(T)+∇_bottom(T)) +
        ρ_ref*k*(L(T).*∇_left(P) + R(T).*∇_right(P)+North(T).*∇_top(P) + B(T).*∇_bottom(P))-
        (ρ_ref^2)*k*(h_bottom.*B(T)-h_top.*North(T)) +
        ρ_ref*k*α*(h_bottom.*B(T).^2 - h_top.*North(T).^2) 

        F1[end,:] += (1/ϵ) * T[end,:] .- (1/ϵ) * bc_bott
        F1[1,:] += (1/ϵ) * T[1,:] .- (1/ϵ) * bc_top

        if heating_from_left
            F1[:,1] += (1/ϵ) * T[:,1] .- (1/ϵ) .* (0.5 .- gridy./300)
        end

        F2 = (∇_left(P)+∇_right(P)+∇_top(P)+∇_bottom(P) )+
        α*(h_bottom.*B(T)-h_top.*North(T))

        F2[1,:] += (1/ϵ) * P[1,:] .- (1/ϵ) * bc_top

        F1 = reshape(F1,N[1]*N[2])
        F2 = reshape(F2,N[1]*N[2])
        return(vcat(F1,F2))
    end
end

function newton(A,b,u0; tol=1.0e-8, maxit=100)

    result=DiffResults.JacobianResult(u0)
    history=Float64[]
    u=copy(u0)
    it=1
    while it<maxit
    ForwardDiff.jacobian!(result,(v)->A(v)-b ,u)
    res=DiffResults.value(result)
    jac=DiffResults.jacobian(result)
    h=jac\res
    u-=h
    nm=norm(h)
    push!(history,nm)
    println(it)
    if nm<tol
    return u,history
    end
    it=it+1
    end
    throw("convergence failed")
end

n_x = 20
n_y = 50
n_fine = n_y÷2
n_coarse = n_y - n_fine

grid1y = LinRange(0.0,2.0,n_fine)
grid2y = LinRange(2.0,150.0,n_coarse+1)[2:end]
gridy = reverse(vcat(grid1y,grid2y)) 

gridx = LinRange(0.0,300.0,n_x)

gridX = gridx'.*ones(size(gridy))
gridY = ones(size(gridx))'.*gridy

bc_bott = 0.5
bc_top  = 0.0

#some exotic boundary conditions

#bc_bott = sin.(2*pi/300*gridx)
#bc_bott = 0.5 .- (gridx.-150.0).^2/150.0^2/2
#bc_bott = 0.5 .- gridx/600

#defining system
sys = new_assemble_nonlinear_system(gridX,gridY,bc_bott,bc_top,heating_from_left = false)

#creating initial guess for steady_state_solution
x_T_0 = zeros(Float64,(n_y,n_x))
x_p_0 = zeros(Float64,n_x*n_y)
x_T_0[end,:] .= bc_bott
x_T_0 = reshape(x_T_0,n_x*n_y)
x_0 = vcat(x_T_0,x_p_0)
b = zeros(Float64,n_x*n_y*2)

#calculating steady_state_solution with newtons method
steady_state_solution,res = newton(sys,b,x_0)

Temperature = reshape(steady_state_solution[1:n_x*n_y],(n_y,n_x))
pressure = reshape(steady_state_solution[n_x*n_y+1:end],(n_y,n_x))

plot_results(gridX, gridY,Temperature,pressure,"put_a_nice_filename_here",savepdf = true,projection = "3d")