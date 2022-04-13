using SparseArrays
using ForwardDiff
using LinearAlgebra
using DifferentialEquations
using Sundials
using NLsolve
using BenchmarkTools
using TimerOutputs
using DelimitedFiles

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

    """
    Computes the gradient/flux
    """
    function gradient(u; side = "left")

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

    function system(x)
        
        T = x[1:N]
        P = x[N+1:end]

        F1 = (λ/c)*( gradient(T; side = "left") + gradient(T; side = "right") ) +
             ρ_ref*k* ( L(T).*gradient(P; side = "left") + R(T).*gradient(P; side = "right") ) -
             (ρ_ref^2)*k*(L(T) - R(T)) +
             ρ_ref*k*α*(L(T).^2 - R(T).^2) 

        F1[1]   += (1/ϵ) * T[1]   - (1/ϵ) * bc_bott
        F1[end] += (1/ϵ) * T[end] - (1/ϵ) * bc_top

        F2 = ( gradient(P; side = "left") + gradient(P; side = "right") ) +
             α*(L(T) - R(T))

        F2[end] += (1/ϵ) * P[end] 

        return(vcat(F1, F2))
    end
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

function create_constraint_function(f,x)
    n = length(x)
    
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