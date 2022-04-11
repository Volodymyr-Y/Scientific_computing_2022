using CairoMakie
using SparseArrays
using ForwardDiff
using LinearAlgebra
using DifferentialEquations
using WGLMakie
WGLMakie.activate!()
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

function assemble_nonlinear_system(grid,bc_bott::AbstractFloat,bc_top::AbstractFloat)
    N = length(grid)
    println(N)
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
    println(1 ./ h_bott)
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
    Z = spzeros((N,N))
    penalty_matrix = spzeros((2*N,2*N))
    penalty_matrix[1,1] = 1.0/ϵ
    penalty_matrix[N,N] = 1.0/ϵ
    penalty_matrix[end,end] = 1.0/ϵ
    penalty_vector = zeros(Float64,2*N)
    penalty_vector[[1,N,2*N]] = (1/ϵ)*[bc_bott,bc_top,0.0]
    #display(penalty_vector)
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
n = 100
n_fine = n÷5
n_coarse = n - n_fine
grid1 = LinRange(0,2,n_fine)
grid2 = LinRange(2,150,n_coarse+1)[2:end]
grid = vcat(grid1,grid2)
bc_bott = 0.5
bc_top = 0.0
#grid = LinRange(0.0,150.0,n)
println(grid)
x_0 = zeros(Float64,n*2)
x_0[1]= bc_bott
sys = assemble_nonlinear_system(grid,bc_bott,bc_top)



solution = newton(sys,zeros(Float64,n*2),x_0; tol=1.0e-8, maxit=100)[1]
#println(solution)
#println(size(solution))
println("plotting")
plot_results(grid,solution,"regular_grid_FV_solver/1D_nonlinear.png")