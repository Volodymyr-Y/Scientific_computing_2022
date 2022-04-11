using BenchmarkTools
using SparseArrays
using ForwardDiff
using LinearAlgebra

function stitch(A::T,B::T,C::T,D::T) where{T<:AbstractMatrix}
    return hcat(vcat(A,C),vcat(B,D))
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
        for i in 2:N
            u1[i] = u[i-1]+u[i] 
        end
        return u1 / 2.0
    end

    function R(u)
        u1 = zeros(typeof(u[1]),N)
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
        ρ_ref*k*α*(L(T).^2-R(T).^2) 

        F1[1] += (1/ϵ) * T[1] - (1/ϵ) * bc_bott
        F1[end] += (1/ϵ) * T[end] - (1/ϵ) * bc_top

        F2 = ( ∇_left(P)+∇_right(P) )+
        α*(L(T)-R(T))

        F2[end] += (1/ϵ) * P[end] 

        return(vcat(F1,F2))
    end
end


function assemble_nonlinear_system(grid::Array{Float64,1},bc_bott::Float64,bc_top::Float64)
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
    #println("assembled matrix")
    #println(typeof( [L Z; Z Z]*ones(Float64,2*N).^2 .*  [L Z; L Z]*ones(Float64,2*N)))
    function system!(x::Array{Float64})
        F = (λ/c)*([∇_right+∇_left Z;Z Z]*x)+
        ρ_ref*k *( [L Z; Z Z]*x .* [Z ∇_left;Z Z]*x +
                 [R Z; Z Z]*x .* [Z ∇_right;Z Z]*x) -
        (ρ_ref^2)*k*(([L Z; Z Z]-[R Z; Z Z]) *x) +
        ρ_ref*k*α*(([L Z; Z Z]*x).^2 - ([R Z; Z Z]*x).^2) #+
        #([Z Z; Z ∇_right+∇_left]*x)+ 
       #α*([Z Z; L-R Z])*x +
        #penalty_matrix*x - penalty_vector
        return F
    end
    return system!
end

function direct_func(x::Array{Float64})
    grid = LinRange(0.0,150.0,20)
    bc_bott = 0.5
    bc_top = 0.0
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
    #DD = Z*grid
    #ASAS = hcat(Z,Z)
    #[∇_right+∇_left Z;Z Z]
    F = (λ/c)*(stitch(∇_right+∇_left, Z,Z, Z)*x)+
    ρ_ref*k *( ((stitch(L ,Z, Z ,Z)*x) .* (stitch(Z, ∇_left,Z ,Z)*x)) +
               ((stitch(R, Z, Z, Z)*x) .* (stitch(Z ,∇_right,Z ,Z)*x)) ) -
    (ρ_ref^2)*k* (stitch(L-R, Z, Z,Z)*x) +
    ρ_ref*k*α*( ((stitch(L,Z, Z,Z)*x) .^2) - ((stitch(R, Z, Z,Z)*x) .^2) ) #+
    #(stitch(Z, Z, Z, ∇_right+∇_left)*x)+ 
    #α*stitch(Z, Z, L-R, Z)*x +
    #penalty_matrix*x - penalty_vector
    return F
end
n = 50
grid = collect(LinRange(0.0,150.0,n))
x_0 = rand(Float64,2*n)
original_func = assemble_nonlinear_system(grid,0.5,0.0)
new_func = new_assemble_nonlinear_system(grid,0.5,0.0)

#@code_warntype new_func(x_0)


#println(ForwardDiff.jacobian(new_func,x_0))
#println(original_func(x_0) - new_func(x_0))

#@btime original_func(x_0)

#@btime new_func(x_0)

#@code_warntype direct_func(x_0)
#@btime direct_func(x_0)

println("done")