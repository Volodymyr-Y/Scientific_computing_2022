using SparseArrays
using SparsityDetection
using Cassette
using ForwardDiff
using LinearAlgebra
using DifferentialEquations
using Sundials
using NLsolve
using BenchmarkTools
using SparseDiffTools
using FiniteDiff
using DelimitedFiles
using TimerOutputs

function plot_results(grid,solution,filename)
    #Makie.inline!(false)
    n = length(grid)
    T = solution[1:n]
    P = solution[n+1:end]
    #println(size(solution)) 
    #println(size(P))    
    #println(size(T))
    #println(size(grid))
    fig = Figure()
    ax1 = CairoMakie.Axis(fig[1,1],ylabel = "temperature")
    ax2 = CairoMakie.Axis(fig[2,1],ylabel = "pressure")
end

function heatmap(matrix)
    n = size(matrix)[1]
    x = 1:n
    y = 1:n
    fig = Figure()
    ax = Axis(fig[1,1])
    heatmap!(x,y,matrix)
    display(fig)
end


include("plotting.jl")

"""
FVM discretization of the considered nonlinear system (heat and mass transfer in a porous medium
heated from below). 

**Input**:

    grid: mesh [TBD: specify format]
    bc_bott: {Float} Temperature at the bottom of the domain.
    bc_top:  {Float} Temperature at the top of the domain.

**Output**:
    Discretized system.
"""        
function new_assemble_nonlinear_system(grid,bc_bott::Number,bc_top::Number)
    println("Program initialization: Done")
    println("\nAssembling the FVM model...")

    N       = length(grid)
    λ       = 0.01
    c       = 0.001
    α       = 0.01
    k       = 100.0
    ρ_ref   = 1.0
    ϵ       = 10.0^(-6)
    h       = grid[2:N] - grid[1:N-1]
    h_top   = vcat(h,     [h[end]]) 
    h_bott  = vcat([h[1]] ,h) 

    """
    Computes the gradient of u.
    """
    function gradient(u; side = "left")
        u_prime = zeros(typeof(u[1]), N)

        if side == "left"
            for i in 2:N
                u_prime[i] = u[i-1]-u[i] 
            end
            return u_prime ./ h_bott

        elseif side == "right"
            for i in 1:(N-1)
                u_prime[i] = u[i+1] - u[i]
            end
            return u_prime ./ h_top
        else
            throw(DomainError(side, "The side selection for the gradient (1D) should be left or right"))
        end
    end

    function L(u)
        u1    = zeros(typeof(u[1]), N)
        u1[1] = u[1]*2.0
        
        for i in 2:N
            u1[i] = u[i-1] + u[i] 
        end
        
        return u1 / 2.0
    end

    function R(u)
        u1      = zeros(typeof(u[1]),N)
        u1[end] = u[end]*2.0

        for i in 1:N-1
            u1[i] = u[i+1]+u[i] 
        end

        return u1 / 2.0
    end

    function system!(du,u)
        T = u[1:N]
        P = u[N+1:end]

        F1 = (λ/c)*( gradient(T; side = "left") + gradient(T; side = "right") ) +
             ρ_ref*k* ( L(T).*gradient(P; side = "left") + R(T).*gradient(P; side = "right") ) -
             (ρ_ref^2)*k*(L(T) - R(T)) +
             ρ_ref*k*α*(L(T).^2 - R(T).^2) 

        F1[1]   += (1/ϵ) * T[1]   - (1/ϵ) * bc_bott
        F1[end] += (1/ϵ) * T[end] - (1/ϵ) * bc_top

        F2 = ( gradient(P; side = "left") + gradient(P; side = "right") ) +
             α*(L(T) - R(T))

        F2[end] += (1/ϵ) * P[end] 
        du[:] = vcat(F1,F2)
        nothing
        
    end
    return system!
end


function create_DAE_function(f,n)
    
    E = spdiagm(0 => vcat(ones(Float64,n),zeros(Float64,n)))
    
    function DAE(du,u,p,t)
        return E*du - f(u,p,t)
    end
    
    return DAE
end

function create_constraint_function(f,x)

    n = length(x)

    function g!(dy,y)
        type1   = typeof(y[1])
        z       = vcat(type1.(x),y)
        dz      = zeros(type1,n*2)
        f(dz,z)
        dy[:]   = dz[n+1:end]
        nothing
    end

    return g!
end

function get_jacobian_prototype(n)
    
    A = spdiagm(0 => ones(Float64,n),-1 => ones(Float64,n-1),1 => ones(Float64,n-1))

    return [A A; A A]
end

"""
Returns the sparse jacobian of function f.
"""
function assemble_sparse_jacobian_function(f, n, jacobian_prototype)

    colors = matrix_colors(jacobian_prototype)
    
    function jacobian(x)
        println(typeof(x[1]))
        
        jac = spdiagm(0 => ones(typeof(x[1]),2*n))
        
        forwarddiff_color_jacobian!(jac, x->f(x,0,0), x, colorvec = colors)
        return jac
    end

    return jacobian
end

function get_time_solution(f,T_0,time_interval)
    
    println("Initializing main solver...")
    
    @timeit to "Main solver" begin

        println("\t Computing sparse Jacobian")

        n = length(T_0) 

        @timeit to "Jacobian" begin
            jac_prot   = get_jacobian_prototype(n) 
            g_jac_prot = spdiagm(0 => ones(Float64,n),-1 => ones(Float64,n-1),1 => ones(Float64,n-1))
            g_colors   = matrix_colors(g_jac_prot)

            g!    = create_constraint_function(f,T_0)
            g_jac = copy(g_jac_prot) 

            j!(jac,x)=      forwarddiff_color_jacobian!(jac, g!, x;
                                                        colorvec = g_colors,
                                                        sparsity = g_jac_prot)


            j!(g_jac,rand(n))                                    
        end
        
        g0 = rand(n)
        
        println("\t Computing consistent initial conditions for pressure") 
        
        dg = OnceDifferentiable(g!, j!, g0, g0, g_jac_prot)
        

        @timeit to "NLSolve for IC" P_0 = nlsolve(dg,g0).zero
        
        u0 = vcat(T_0, P_0)
        E  = spdiagm(0 => vcat(ones(Float64,n),zeros(Float64,n)))
        
        # DAE = create_DAE_function(f,n)
        
        println("\t Solving DAE...")

        @timeit to "ODE/DAE solver" begin
            @timeit to "assemling" begin
                FF  = ODEFunction((du,u,p,t) -> f(du,u),mass_matrix=E,jac_prototype=jac_prot)
                prob = ODEProblem(FF,u0,time_interval)
            end

            @timeit to "solving" sol = solve(prob,Rodas5(),dt=0.1,adaptive = false)


        end
    end
     
    return sol
end

function flux(pressure, temperature, grid; side = "top")

    rho_ref = 1.0
    k       = 100
  
    function meshgrid(xin, yin)
      nx   =   length(xin)
      ny   =   length(yin)
      xout =   zeros(ny,nx)
      yout =   zeros(ny,nx)
      
      for jx = 1:nx
          for ix = 1:ny
              xout[ix,jx] = xin[jx]
              yout[ix,jx] = yin[ix]
          end
      end
  
      return (x=xout, y=yout)
    end
  
    # check with Ferhat!
    function gradient(u; side = "left")
      
      ny, nx = N
      u_prime = zeros(typeof(u[1]), N)
  
      if side == "left"
        for i in 2:nx
          for j in 1:ny
              u_prime[j, i] = u[j, i-1] - u[j, i] 
          end  
        end  
        return u_prime ./ h_left
  
      elseif side == "right"
        for i in 1:(nx-1)
          for j in 1:ny
            u_prime[j, i] = u[j, i+1] - u[j, i] 
          end
        end
        return u_prime ./ h_right
  
      elseif side == "bottom"
        for j in 2:ny
          for i in 1:nx
            u_prime[j, i] = u[j-1, i] - u[j, i] 
          end    
        end
        return u_prime ./ h_bottom
  
      elseif side == "top"
        for j in 1:(ny-1)
          for i in 1:nx
            u_prime[j, i] = u[j+1, i]-u[j, i] 
          end
        end
        return u_prime ./ h_top
      else
          throw(DomainError(side, "The side selection for the gradient (1D) should be left or right"))
      end
    end
  
    function rho(T; rho_ref = 1.0, T_ref = 0.0, alpha = 0.01)
      return rho_ref .- alpha .* (T .- T_ref)
    end
  
    n  = length(grid)
    nx = 300
  
    x  = LinRange(0.0, 300.0, nx)
  
    N = (n, nx)
  
    xx, yy = meshgrid(x, grid)
  
    h_x       = xx[:,2:N[2]] - xx[:,1:N[2]-1]     
    h_y       = yy[1:(N[1]-1),:] - yy[2:N[1],:]  
    h_top     = vcat(h_y[1,:]',h_y) 
    h_bottom  = vcat(h_y,h_y[end,:]')  
    h_left    = hcat(h_x[:,1], h_x) 
    h_right   = hcat(h_x, h_x[:,end])  
  
    
    gradP = gradient(pressure; side = side)
  
    if side == "top" || side == "bottom"
      densityT = rho(temperature)
    else
      densityT = zeros(typeof(gradP[1]), size(gradP))
    end

    return -rho_ref .* k .* (gradP - densityT), xx, yy

end


# ================================================================================

to = TimerOutput()

n        = 1000
nx       = 300
n_fine   = n÷5
n_coarse = n - n_fine

grid1 = LinRange(0.0,2.0,n_fine)
grid2 = LinRange(2.0,150.0,n_coarse+1)[2:end]
grid  = vcat(grid1,grid2)

bc_top  = 0.0
bc_bott = 10.0
t_end   = 10.0

x_0     = zeros(Float64, n*2)
x_0[1]  = bc_bott

placeholder = ones(Float64,2*n)

function headingText()
    println("
            ╔═══╦═══╦═══╦═══╦═══╗
            ║╔═╗║╔═╗║╔═╗║╔═╗║╔══╝
            ║║─╚╣║─║║╚══╣╚══╣╚══╗
            ║║─╔╣║─║╠══╗╠══╗║╔══╝
            ║╚═╝║╚═╝║╚═╝║╚═╝║╚══╗
            ╚═══╩═══╩═══╩═══╩═══╝\n")

    println("FVM-based solver of a heat and mass transfer problem in a porous medium\n")
    println("Grid size: ($(nx) x $(n))")
    println("Temporal domain: (0, $(t_end)) s \n")
end

headingText()

@timeit to "FVM assembly" begin
    sys! = new_assemble_nonlinear_system(grid, bc_bott, bc_top)
    sys!(placeholder,rand(2*n))
end

T_0 = LinRange(0.5,0.0,n)
time_solution = get_time_solution(sys!,T_0, (0.0,t_end))


println("\t Done. Calculating the water mass flux.")
# Compute and plot flux
for i in 1:nx
    T[:, i]  = time_solution(1.0)[1:n]
    P[:, i]  = time_solution(1.0)[n+1:end]
end

fluxY, xx, yy = flux(P, T, grid)
fluxX, _, _   = flux(P, T, grid; side="left")


println("\n \t Computation finished. Summary of results:\n")

display(time_solution)
println("\n Performance report:")
display(to)

# plot_flux(xx, yy,fluxY, fluxX; save = true)

# plot_results(grid, solution(0), nx; savepdf = true, projection = "3d")

# Work in progress. In the meanwhile, I'll export the solution as csv and animate the plot in Python (easiest solution)
#
# animate_solution(grid, solution, nx; specie = "temperature", projection = "2d")
writedlm("T_10.csv",  time_solution, ',')