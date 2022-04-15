using CairoMakie
using SparseArrays
using ForwardDiff
using LinearAlgebra
using DifferentialEquations
using GLMakie
using Sundials
using NLsolve
using BenchmarkTools
using TimerOutputs
using DelimitedFiles

GLMakie.activate!()
#include("plotting.jl")

function new_assemble_nonlinear_system(gridX,gridY,bc_bott::Number,bc_top::Number)
    N = (n_y,n_x)
    λ = 0.01
    c = 0.001
    α = 0.01
    k = 100.0
    ρ_ref = 1.0
    ϵ = 10.0^(-6)
    h_x = gridX[:,2:N[2]] - gridX[:,1:N[2]-1]     
    h_y = gridY[1:(N[1]-1),:] - gridY[2:N[1],:]  
    h_left = vcat(h_y[1,:]',h_y) # array of distances between collocation points of size N shifted to the top
    h_right = vcat(h_y,h_y[end,:]')  # array of distances between collocation points of size N shifted to the bottom
    h_bottom = hcat(h_x[:,1],h_x) # array of distances between collocation points of size N shifted to the left
    h_top = hcat(h_x,h_x[:,end])  # array of distances between collocation points of size N shifted to the right

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

        F2 = (∇_left(P)+∇_right(P)+∇_top(P)+∇_bottom(P) )+
        α*(h_bottom.*B(T)-h_top.*North(T))

        F2[1,:] += (1/ϵ) * P[1,:] 

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

n_x = 10
n_y = 50
n_fine = n_y÷2
n_coarse = n_y - n_fine
grid1y = LinRange(0.0,2.0,n_fine)
grid2y = LinRange(2.0,150.0,n_coarse+1)[2:end]
gridy = reverse(vcat(grid1y,grid2y)) #Uncomment this for nonuniform grid but also comment uniform one (does not work yet)
gridx = LinRange(0,300.0,n_x)
#gridy = reverse(LinRange(0,150.0,n_y)) #uniform grid


gridX = gridx'.*ones(n_y)
gridY = ones(n_x)'.*gridy

bc_bott = 0.5
bc_top  = 0.0
t_end   = 10.0

#defining system
sys = new_assemble_nonlinear_system(gridX,gridY,bc_bott,bc_top)

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

fig = Figure()
println("Solution obtained, start plotting")
ax = Axis3(fig[1,1]; aspect=(1, 1, 1),xlabel = "x",ylabel = "y",zlabel = "pressure")
hm = surface!(ax,gridX, gridY, pressure )
Colorbar(fig[1, 2],hm )
display(fig)

###TIMEDEPENDENT_PROBLEM###
# """
# Creates a DAE system f(du,u,p,t) = 0
# """
# function create_DAE_function(f,n) 
#     E = spdiagm(0 => vcat(ones(Float64,n),zeros(Float64,n)))
    
#     function DAE(du,u,p,t)
#         return E*du - f(u)
#     end
    
#     return DAE
# end

# function create_constraint_function(f,x)
#     n = length(x)
    
#     function a(y)
#         z = vcat(x,y)
#         return f(z)[n+1:end]
#     end
    
#     return a
# end

# function get_time_solution(f, T_0, time_interval)

#     println("Solver initialized")

#     @timeit to "system solution" begin
#         n   = length(T_0) 
#         DAE = create_DAE_function(f,n)     
#         g   = create_constraint_function(f, T_0)

#         write(stdin.buffer, 0x0C)
#         println("\t Computing consistent initial conditions for pressure") 
#         @timeit to "newton's method" P_0 = newton(g,zeros(Float64,n),zeros(Float64,n); tol=1.0e-8, maxit=100)[1]
#         u0  = vcat(T_0, P_0)
#         du0 = f(u0)

#         println("\n\t Solving DAE...")

#         @timeit to "DAE solver" begin
#             @timeit to "assmebling" prob = DAEProblem(DAE, du0, u0, time_interval)
#             @timeit to "solving" sol  = solve(prob, IDA(), dt=0.1; verbose = true)
#         end
        
#         println("\nSolution available. Summary and solution:\n")

#         return sol
#     end
# end

# # Create the timer object
# to = TimerOutput()

# function headingText()
#     println("
#             ╔═══╦═══╦═══╦═══╦═══╗
#             ║╔═╗║╔═╗║╔═╗║╔═╗║╔══╝
#             ║║─╚╣║─║║╚══╣╚══╣╚══╗
#             ║║─╔╣║─║╠══╗╠══╗║╔══╝
#             ║╚═╝║╚═╝║╚═╝║╚═╝║╚══╗
#             ╚═══╩═══╩═══╩═══╩═══╝\n")

#     println("FVM-based solver of a heat and mass transfer problem in a porous medium\n")
#     println("Grid size: ($(n_x) x $(n_y))")
#     println("Temporal domain: (0, $(t_end)) s \n")
# end

# headingText()

# T_0 = x_T_0
# t_end = 1
# solution = get_time_solution(sys, T_0, (0.0,t_end));

# display(solution)
# println("\n Performance report:")
# display(to)

# solution_at_t = solution(0.01)

# Temperature = reshape(solution_at_t[1:n_x*n_y],(n_y,n_x))
# pressure = reshape(solution_at_t[n_x*n_y+1:end],(n_y,n_x))

# fig = Figure()
# println("Solution obtained, start plotting")
# ax = Axis3(fig[1,1]; aspect=(1, 1, 1),xlabel = "x",ylabel = "y",zlabel = "pressure")
# hm = surface!(ax,gridX, gridY, Temperature )
# Colorbar(fig[1, 2],hm )
# display(fig)

# # Work in progress. In the meanwhile, I'll export the solution as csv and animate the plot in Python (easiest solution)
# #
# # animate_solution(grid, solution, nx; specie = "temperature", projection = "2d")

# writedlm("solution.csv",  solution, ',')