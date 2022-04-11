using CairoMakie
using SparseArrays
using ForwardDiff
using LinearAlgebra
using DifferentialEquations
using GLMakie
using Sundials
using NLsolve
using BenchmarkTools
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

    qw = lines!(ax1,grid, T)
    we = lines!(ax2,grid, P)
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

    function system(x)
        T = x[1:N]
        P = x[N+1:end]

        F1 = (λ/c)*(∇_left(T)+∇_right(T)) +
        ρ_ref*k*(L(T).*∇_left(P) + R(T).*∇_right(P))-
        (ρ_ref^2)*k*(L(T)-R(T)) +
        ρ_ref*k*α*(L(T).^2 - R(T).^2) 

        F1[1] += (1/ϵ) * T[1] - (1/ϵ) * bc_bott
        F1[end] += (1/ϵ) * T[end] - (1/ϵ) * bc_top

        F2 = ( ∇_left(P)+∇_right(P) )+
        α*(L(T)-R(T))

        F2[end] += (1/ϵ) * P[end] 

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
    if nm<tol
    return u,history
    end
    it=it+1
    end
    throw("convergence failed")
end

function assemble_nonlinear_system(grid,bc_bott,bc_top)
    N = length(grid)
    #println(N)
    λ = 0.01
    c = 0.001
    α = 0.01
    k = 100.0
    ρ_ref = 1.0
    #bc_bott = 0.5
    #bc_top = 0.0
    ϵ = 10.0^(-6)
    #x = LinRange(0,150,N)
    #x = Float64[0,100,110,120,150]
    h = grid[2:N] - grid[1:N-1]    # array of distances between collocation points of size N-1
    h_top = vcat(h,[h[end]]) # array of distances between collocation points of size N shifted to the top
    h_bott = vcat([h[1]],h)  # array of distances between collocation points of size N shifted to the bottom
    #println(1 ./ h_bott)
    # Let's say u is a species such as pressure or temperature 
    # ∇_left is a matrix such that ∇_left @ u will yeild ∇u (flow) throgh left boundary in each cell
    ∇_left = spdiagm(-1 => ones(Float64,N-1)./h , 0 => -ones(Float64,N)./ h_top)
    ∇_left[1,1] = 0.0
    #display(∇_left)
    # ∇_right is a matrix such that ∇_left @ u will yeild ∇u (flow) throgh right boundary in each cell
    ∇_right = spdiagm(1 => ones(Float64,N-1) ./h , 0 => -ones(Float64,N) ./ h_bott)
    ∇_right[end,end] = 0.0
    #display(∇_right+∇_left)
    L = 0.5*spdiagm(-1 => ones(Float64,N-1) , 0 => ones(Float64,N))
    L[1,1] = 1.0
    R = 0.5*spdiagm(1 => ones(Float64,N-1) , 0 => ones(Float64,N))
    R[end,end] = 1.0
    #display(L-R)
    Z = spzeros(Float64,(N,N))
    penalty_matrix = spzeros(Float64,(2*N,2*N))
    penalty_matrix[1,1] = 1.0/ϵ
    penalty_matrix[N,N] = 1.0/ϵ
    penalty_matrix[end,end] = 1.0/ϵ
    penalty_vector = zeros(Float64,2*N)
    penalty_vector[[1,N,2*N]] = (1/ϵ)*[bc_bott,bc_top,0.0]
    println("assembled matrix")
    println(typeof( [L Z; Z Z]*ones(Float64,2*N).^2 .*  [L Z; L Z]*ones(Float64,2*N)))
    function system!(x)
        F = (λ/c)*([∇_right+∇_left Z;Z Z]*x)+
        ρ_ref*k *( [L Z; Z Z]*x .* [Z ∇_left;Z Z]*x +
                 [R Z; Z Z]*x .* [Z ∇_right;Z Z]*x) -
        (ρ_ref^2)*k*(([L Z; Z Z]-[R Z; Z Z]) *x) +
        ρ_ref*k*α*(([L Z; Z Z]*x).^2 - ([R Z; Z Z]*x).^2) +
        ([Z Z; Z ∇_right+∇_left]*x)+ 
        α*([Z Z; L-R Z])*x +
        penalty_matrix*x - penalty_vector
        
        return F
    end
    return system!
end

function create_DAE_function(f,n) # create a DAE function f(du,u,p,t) = 0
    E = spdiagm(0 => vcat(ones(Float64,n),zeros(Float64,n)))
    function DAE(du,u,p,t)
        return E*du - f(u)
    end
    return DAE
end


function create_constraint_function(f,x)
    n = length(x)
    println(n)
    function a(y)
        z = vcat(x,y)
        return f(z)[n+1:end]
    end
    return a
end

function get_time_solution(f,T_0,time_interval)
    n = length(T_0) # number of collocation points
    DAE = create_DAE_function(f,n) # create DAE function
    g = create_constraint_function(f,T_0)
    println("computing consistent initial conditions for pressure") 
    P_0 = newton(g,zeros(Float64,n),zeros(Float64,n); tol=1.0e-8, maxit=100)[1]
    println("done") 
    u0 = vcat(T_0,P_0)
    du0 = f(u0)
    #println(du0)
    println("start solving")
    prob = DAEProblem(DAE,du0,u0,time_interval)
    @time sol = solve(prob,IDA(),dt=0.1)
    println("done")
    return sol
end

n = 1000
n_fine = n÷5
n_coarse = n - n_fine
grid1 = LinRange(0.0,1.0,n_fine)
grid2 = LinRange(1.0,150.0,n_coarse+1)[2:end]
grid = vcat(grid1,grid2)
bc_bott = 0.5
bc_top = 0.0
t_end = 10.0
#grid = LinRange(0.0,150.0,n)
println("size of grid: ",length(grid))
x_0 = zeros(Float64,n*2)
x_0[1]= bc_bott
sys = new_assemble_nonlinear_system(grid,bc_bott,bc_top)

#@code_warntype assemble_nonlinear_system(grid,bc_bott,bc_top)
#@code_warntype sys(x_0)
#sys(x_0)
T_0 = LinRange(0.5,0.0,n)
#solution = get_time_solution(sys,T_0,(0.0,t_end))

create_animatrion(solution,grid,t_end,"regular_grid_FV_solver/1D_nonlinear.gif")
#println("calculating solution")
#@time solution = newton(sys,zeros(Float64,n*2),x_0; tol=1.0e-8, maxit=100)[1]
#@btime solution = nlsolve(sys,x_0,autodiff = :forward).zero


#println("done")
#println("ploting")
#plot_results(grid,solution,"regular_grid_FV_solver/1D_nonlinear.png")