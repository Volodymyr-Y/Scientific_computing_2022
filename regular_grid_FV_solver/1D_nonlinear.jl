using SparseArrays
using SparsityDetection
using Cassette
using ForwardDiff
using LinearAlgebra
using DifferentialEquations
using Sundials
using NLsolve
using BenchmarkTools
<<<<<<< HEAD
using SparseDiffTools
using FiniteDiff
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
=======
using TimerOutputs
using DelimitedFiles
>>>>>>> b0ea7cb8457c8d55c042979295cc8f4147258ba4

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
function new_assemble_nonlinear_system(grid, bc_bott::Number, bc_top::Number)
    
    N       = length(grid)
    λ       = 0.01
    c       = 0.001
    α       = 0.01
    k       = 100.0
    ρ_ref   = 1.0
    ϵ       = 10.0^(-6)
    
    h       = grid[2:N] - grid[1:N-1]
    h_top   = vcat(h, [h[end]])          
    h_bott  = vcat([h[1]], h)            

<<<<<<< HEAD
function heatmap(matrix)
    n = size(matrix)[1]
    x = 1:n
    y = 1:n
    fig = Figure()
    ax = Axis(fig[1,1])
    heatmap!(x,y,matrix)
    display(fig)
end


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
=======
    """
    Computes the gradient/flux
    """
    function gradient(u; side = "left")
>>>>>>> b0ea7cb8457c8d55c042979295cc8f4147258ba4

        u_prime = zeros(typeof(u[1]), N)

        if side == "left"
            for i in 2:(N)
                u_prime[i] = u[i-1] - u[i]
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

<<<<<<< HEAD
    function system!(du,u)
        T = u[1:N]
        P = u[N+1:end]
=======
    function system(x)
        
        T = x[1:N]
        P = x[N+1:end]
>>>>>>> b0ea7cb8457c8d55c042979295cc8f4147258ba4

        F1 = (λ/c)*( gradient(T; side = "left") + gradient(T; side = "right") ) +
             ρ_ref*k* ( L(T).*gradient(P; side = "left") + R(T).*gradient(P; side = "right") ) -
             (ρ_ref^2)*k*(L(T) - R(T)) +
             ρ_ref*k*α*(L(T).^2 - R(T).^2) 

        F1[1]   += (1/ϵ) * T[1]   - (1/ϵ) * bc_bott
        F1[end] += (1/ϵ) * T[end] - (1/ϵ) * bc_top

        F2 = ( gradient(P; side = "left") + gradient(P; side = "right") ) +
             α*(L(T) - R(T))

        F2[end] += (1/ϵ) * P[end] 
<<<<<<< HEAD
        du[:] = vcat(F1,F2)
        nothing
        
=======

        return(vcat(F1, F2))
>>>>>>> b0ea7cb8457c8d55c042979295cc8f4147258ba4
    end
    return system!
end


<<<<<<< HEAD
function create_DAE_function(f,n) # create a DAE function f(du,u,p,t) = 0
=======
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
>>>>>>> b0ea7cb8457c8d55c042979295cc8f4147258ba4
    E = spdiagm(0 => vcat(ones(Float64,n),zeros(Float64,n)))
    
    function DAE(du,u,p,t)
        return E*du - f(u,p,t)
    end
    
    return DAE
end

function create_constraint_function(f,x)
    n = length(x)
<<<<<<< HEAD

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

function get_jacobian_prototype(n)
    A = spdiagm(0 => ones(Float64,n),-1 => ones(Float64,n-1),1 => ones(Float64,n-1))
    return [A A ; A A]
end

function assemble_sparse_jacobian_function(f,n,jacobian_prototype)
    # this function returns sparse jacobian for function f
    colors = matrix_colors(jacobian_prototype)
    function jacobian(x)
        println(typeof(x[1]))
        jac = spdiagm(0 => ones(typeof(x[1]),2*n))
        forwarddiff_color_jacobian!(jac,x->f(x,0,0), x, colorvec = colors)
        return jac
    end
    return jacobian
    
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

n = 1000
n_fine = n÷5
n_coarse = n - n_fine
grid1 = LinRange(0.0,2.0,n_fine)
grid2 = LinRange(2.0,150.0,n_coarse+1)[2:end]
grid = vcat(grid1,grid2)
bc_bott = 0.5
bc_top = 0.0
t_end = 10.0
#grid = LinRange(0.0,150.0,n)
println("size of grid: ",length(grid))
x_0 = zeros(Float64,n*2)
x_0[1]= bc_bott
placeholder = ones(Float64,2*n)
sys! = new_assemble_nonlinear_system(grid,bc_bott,bc_top)
sys!(placeholder,rand(2*n))
#println(placeholder)
#jac = assemble_sparse_jacobian_function(sys,n,get_jacobian_prototype(n))
#jac(rand(2*n))
#println(sys(x_0))
#matrix = ForwardDiff.jacobian(x -> sys(x,0,0),rand(Float64,n*2))
#display(get_jacobian_prototype(n))
#display(sparse(matrix))
#dense_to_sparse(matrix)
#heatmap(matrix)
#@code_warntype assemble_nonlinear_system(grid,bc_bott,bc_top)
#@code_warntype sys(x_0)
#sys(x_0)

#display(get_jacobian_prototype(n))
T_0 = LinRange(0.5,0.0,n)
time_solution = get_time_solution(sys!,T_0,(0.0,t_end))

create_animatrion(time_solution,grid,t_end,"regular_grid_FV_solver/1D_nonlinear.gif")
#println("calculating solution")
#@time solution = newton(sys,zeros(Float64,n*2),x_0; tol=1.0e-8, maxit=100)[1]
#solution = nlsolve(sys,x_0,autodiff = :forward).zero
=======
    
    function a(y)
        z = vcat(x,y)
        return f(z)[n+1:end]
    end
    
    return a
end

function get_time_solution(f, T_0, time_interval)

    println("Solver initialized")

    @timeit to "system solution" begin
        n   = length(T_0) 
        DAE = create_DAE_function(f,n)     
        g   = create_constraint_function(f, T_0)

        write(stdin.buffer, 0x0C)
        println("\t Computing consistent initial conditions for pressure") 
        @timeit to "newton's method" P_0 = newton(g,zeros(Float64,n),zeros(Float64,n); tol=1.0e-8, maxit=100)[1]
        u0  = vcat(T_0, P_0)
        du0 = f(u0)

        println("\n\t Solving DAE...")

        @timeit to "DAE solver" begin
            @timeit to "assmebling" prob = DAEProblem(DAE, du0, u0, time_interval)
            @timeit to "solving" sol  = solve(prob, IDA(), dt=0.1; verbose = true)
        end
        
        println("\nSolution available. Summary and solution:\n")

        return sol
    end
end

# Create the timer object
to = TimerOutput()

n           = 1000
nx          = 100
n_fine      = n÷5
n_coarse    = n - n_fine

grid1       = LinRange(0.0,1.0,n_fine)
grid2       = LinRange(1.0,150.0,n_coarse+1)[2:end]
grid        = vcat(grid1,grid2)
>>>>>>> b0ea7cb8457c8d55c042979295cc8f4147258ba4

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

T_0      = LinRange(0.5,0.0,n)
solution = get_time_solution(sys, T_0, (0.0,t_end));

display(solution)
println("\n Performance report:")
display(to)

plot_results(grid, solution(3.3), nx; savepdf = false, projection = "3d")

# Work in progress. In the meanwhile, I'll export the solution as csv and animate the plot in Python (easiest solution)
#
# animate_solution(grid, solution, nx; specie = "temperature", projection = "2d")

writedlm("solution.csv",  solution, ',')