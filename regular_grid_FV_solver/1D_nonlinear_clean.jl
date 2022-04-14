using CairoMakie
using SparseArrays
using ForwardDiff
using LinearAlgebra
using DifferentialEquations
using GLMakie
using PyPlot
using Sundials
using NLsolve
using BenchmarkTools
using TimerOutputs

function steady_solution(f!,N,bc_bott,bc_top)
    jac_prot = get_jacobian_prototype(N) 
    colors = matrix_colors(jac_prot)
    j!(jac,x)=      forwarddiff_color_jacobian!(jac,f!, x;
                                        colorvec = colors,
                                        sparsity = jac_prot)
    
    x_0 = zeros(2*N) # initial guess
    x_0[1] = bc_bott # BC's are needed for better initial guess
    x_0[N] = bc_top
    dg = OnceDifferentiable(f!, j!, x_0, x_0, jac_prot)
    @time sol = nlsolve(dg,x_0).zero
    return sol
end

# To be repeated--------------------------------------------
function plot_results(grid, solution, nx)

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
    
    for i in 1:nx
        T[:, i]  = solution[1:n]
        P[:, i]  = solution[n+1:end]
    end

    PyPlot.clf()
    PyPlot.plot_surface(grid, x, transpose(T), cmap = PyPlot.cm.viridis)
    PyPlot.xlabel("y (-)")
	PyPlot.ylabel("x (-)")
	PyPlot.zlabel("temperature (-)")
    PyPlot.title(" t = 0.26 s")

    PyPlot.gcf()
end

# To be repeated
function create_animatrion(solution,grid,t_end,filename)
    n = length(grid)
    P_max = max(solution(0.0)[n+1],solution(0.0)[end],solution(t_end)[n+1],solution(t_end)[end])
    println("pmax ",P_max)
    T_max = max(solution(0.0)[1],solution(t_end)[n])
    N_frames = 7*30
    fps = 30
    t_sample = LinRange(0,t_end,N_frames)
    fig = Figure(resolution = (1200, 800))
    time = Observable(0.0)
    u = @lift(solution($time))
    T = lift(a -> a[1:n] ,u)
    P = lift(a -> a[n+1:end] ,u)
    
    ax1 = GLMakie.Axis(fig[1,1],ylabel = "temperature",title = @lift("t = $(round(($time), digits = 1))"*" s"))
    ax2 = CairoMakie.Axis(fig[2,1],ylabel = "pressure")
    ylims!(ax2,(0,P_max*1.1))
    lines!(ax1,grid, T)
    lines!(ax2,grid, P)
    t_sample = LinRange(0,t_end,N_frames)
    record(fig,filename, t_sample;framerate = fps) do t
    time[] = t
    end
end
# ----------------------------------------------------------


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
    N = length(grid)
    #println(N)
    λ = 0.01
    c = 0.001
    α = 0.01
    k = 100.0
    ρ_ref = 1.0
    ϵ = 10.0^(-6)
    h = grid[2:N] - grid[1:N-1]
    h_top = vcat(h,[h[end]]) # array of distances between collocation points of size N shifted to the top
    h_bott = vcat([h[1]],h)  # array of distances between collocation points of size N shifted to the bottom

    function ∇_left(u)
        u1 = zeros(typeof(u[1]),N)
        for i in 2:N
            u1[i] = u[i-1]-u[i] 
        end
        return u1 ./ h_bott
    end

    function ∇_right(u)
        u1 = zeros(typeof(u[1]),N)
        for i in 1:N-1
            u1[i] = u[i+1]-u[i] 
        end
        return u1 ./ h_top
    end

    function L(u)
        u1 = zeros(typeof(u[1]),N)
        u1[1] = u[1]*2
        for i in 2:N
            u1[i] = u[i-1]+u[i] 
        end
        return u1 / 2.0
    end

    function R(u)
        u1 = zeros(typeof(u[1]),N)
        u1[end] = u[end]*2
        for i in 1:N-1
            u1[i] = u[i+1]+u[i] 
        end
        return u1 / 2.0
    end

    function system!(du,u)
        T = u[1:N]
        P = u[N+1:end]

        F1 = (λ/c)*(∇_left(T)+∇_right(T)) +
        ρ_ref*k*(L(T).*∇_left(P) + R(T).*∇_right(P))-
        (ρ_ref^2)*k*(L(T)-R(T)) +
        ρ_ref*k*α*(L(T).^2 - R(T).^2) 

        F1[1] += (1/ϵ) * T[1] - (1/ϵ) * bc_bott
        F1[end] += (1/ϵ) * T[end] - (1/ϵ) * bc_top

        F2 = ( ∇_left(P)+∇_right(P) )+
        α*(L(T)-R(T))

        F2[end] += (1/ϵ) * P[end] 
        du[:] = vcat(F1,F2)
        nothing
        
    end
    return system!
end



"""
Newton's method algorithm for finding the roots of differentiable
nonlinear funtions using iterated local linearization of a function to
approximate its roots
"""
function newton(A, b, u0; tol= 1.0e-8, maxit= 100)
    result  =   DiffResults.JacobianResult(u0)
    history =   Float64[]
    u       =   copy(u0)
    it      =   1
    
    while it < maxit
        
        write(stdin.buffer, 0x0C)
        println("\t Newton's method iteration: $(it)")
        
        ForwardDiff.jacobian!(result,(v)->A(v)-b ,u)
        
        res =   DiffResults.value(result)
        jac =   DiffResults.jacobian(result)
        h   =   jac\res
        u  -=   h
        nm  =   norm(h)

        push!(history,nm)
        
        if nm < tol
            println("\t Newton's method converged")
            return u, history
        end
        it  = it+1
    end

    throw("Convergence of the Newton's method error. Solver diverged.")
end


"""
Creates a DAE system f(du,u,p,t) = 0
"""
function create_DAE_function(f,n) 
    E = spdiagm(0 => vcat(ones(Float64,n),zeros(Float64,n)))
    
    function DAE(du,u,p,t)
        return E*du - f(u)
    end
    
    return DAE
end


# what is this doing?
function create_constraint_function(f,x)
    n = length(x)

    function g!(dy,y)
        type1 = typeof(y[1])
        z = vcat(type1.(x),y)
        dz = zeros(type1,n*2)
        f(dz,z)
        dy[:] = dz[n+1:end]
        nothing
    end

    return g!
end

function get_time_solution(f,T_0,time_interval)
    n = length(T_0) # number of collocation points
    
    # jacobian prototype for time solution in 1d
    jac_prot = get_jacobian_prototype(n) 

    # jacobian prototype for constraint function g() in 1d
    g_jac_prot = spdiagm(0 => ones(Float64,n),-1 => ones(Float64,n-1),1 => ones(Float64,n-1))

    # jacobian coloring of g()
    g_colors = matrix_colors(g_jac_prot)

    #
    g! = create_constraint_function(f,T_0)

 
    g_jac = copy(g_jac_prot) # this creates a mutating container for the jacobian of g()

    j!(jac,x)=      forwarddiff_color_jacobian!(jac,g!, x;
                                        colorvec = g_colors,
                                        sparsity = g_jac_prot)


    j!(g_jac,rand(n))                                    
    display(g_jac)
    g0 = rand(n)
    println("computing consistent initial conditions for pressure") 
    dg = OnceDifferentiable(g!, j!, g0, g0, g_jac_prot)
    @time P_0 = nlsolve(dg,g0).zero
    #@time P_0 = nlsolve(g!,j!,g0,autodiff = :forward).zero #
    #println(P_0)
    println("done") 

    u0 = vcat(T_0,P_0)
    E = spdiagm(0 => vcat(ones(Float64,n),zeros(Float64,n)))
    DAE = create_DAE_function(f,n) # create DAE function
    #println(du0)
    println("start solving")
    #prob = DAEProblem(DAE,du0,u0,time_interval)
    #@time sol = solve(prob,IDA(),dt=0.2)
    FF= ODEFunction((du,u,p,t) -> f(du,u),mass_matrix=E,jac_prototype=jac_prot)
    prob = ODEProblem(FF,u0,time_interval)
    @time sol = solve(prob,Rodas5(),dt=0.1,adaptive = false)
    println("done")
    
    return sol
end


# Create the timer object
to = TimerOutput()

n           = 500
nx          = 300
n_fine      = n÷5
n_coarse    = n - n_fine

grid1       = LinRange(0.0,1.0,n_fine)
grid2       = LinRange(1.0,150.0,n_coarse+1)[2:end]
grid        = vcat(grid1,grid2)

bc_bott     = 0.5
bc_top      = 0.0
t_end       = 10.0

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

x_0     = zeros(Float64, n*2)
x_0[1]  = bc_bott

@timeit to "system assembling" sys     = new_assemble_nonlinear_system(grid, bc_bott, bc_top)
steady_sol = steady_solution(sys,n,bc_bott,bc_top)



T_0      = LinRange(0.5,0.0,n)
#solution = get_time_solution(sys, T_0, (0.0,t_end));

#display(solution)
#println("\n Performance report:")
#display(to)

# create_animatrion(solution,grid,t_end,"regular_grid_FV_solver/1D_nonlinear.gif")
println(steady_sol)
plot_results(grid, steady_sol, nx)