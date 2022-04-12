using CairoMakie
using SparseArrays
using ForwardDiff
using LinearAlgebra
using DifferentialEquations
using GLMakie
using Sundials
using NLsolve
using BenchmarkTools

n = 30
n_fine = n÷5
n_coarse = n - n_fine
grid1y = LinRange(0.0,1.0,n_fine)
grid2y = LinRange(1.0,150.0,n_coarse+1)[2:end]
gridy = vcat(grid1y,grid2y)
gridx = LinRange(0,300.0,n_coarse)


gridX = gridx'.*ones(n_fine+n_coarse)
gridY = ones(n)'.*gridy

bc_bott = 0.5
bc_top = 0.0
t_end = 10.0


function new_assemble_nonlinear_system(gridX,gridY,bc_bott::Number,bc_top::Number)
    N = size(gridX)
    λ = 0.01
    c = 0.001
    α = 0.01
    k = 100.0
    ρ_ref = 1.0
    ϵ = 10.0^(-6)
    h_x = gridX[:,2:N[2]] - gridX[:,1:N[2]-1]     
    h_y = gridY[2:N[1],:] - gridY[1:(N[1]-1),:]  
    h_top = vcat(h_y,h_y[1,:]') # array of distances between collocation points of size N shifted to the top
    h_bott = vcat(h_y[end,:]',h_y)  # array of distances between collocation points of size N shifted to the bottom

    h_left = hcat(h_x[:,1],h_x) # array of distances between collocation points of size N shifted to the left
    h_right = hcat(h_x,h_x[:,end])  # array of distances between collocation points of size N shifted to the right

    function ∇_left(u)
        u1 = zeros(typeof(u[1]),N)
        println(size(u1))
        for i in 2:N[1]
            for j in 1:N[2]
                u1[i,j] = u[i-1,j]-u[i,j] 
            end    
        end
        return u1 ./ h_left
    end

    function ∇_right(u)
        u1 = zeros(typeof(u[1]),N)
        for i in 1:N[1]-1
            for j in 1:N[2]
                u1[i,j] = u[i+1,j]-u[i,j] 
            end
        end
        return u1 ./ h_right
    end

    function ∇_bottom(u)
        u1 = zeros(typeof(u[1]),N)
        for j in 2:N[2]
            for i in 1:N[1]
                u1[i,j] = u[i,j-1]-u[i,j] 
            end    
        end
        return u1 ./ h_bottom
    end

    function ∇_top(u)
        u1 = zeros(typeof(u[1]),N)
        for j in 1:N[2]-1
            for i in 1:N[1]
                u1[i,j] = u[i,j+1]-u[i,j] 
            end
        end
        return u1 ./ h_top
    end

    function L(u)
        u1 = zeros(typeof(u[1]),N)
        u1[:,1] *= 2.0
        for i in 2:N[1]
            for j in 1:N[2]
                u1[i,j] = u[i-1,j]+u[i,j] 
            end
        end
        return u1 / 2.0
    end

    function R(u)
        u1 = zeros(typeof(u[1]),N)
        u1[:,end] *= 2.0
        for i in 1:N[1]-1
            for j in 1:N[2]
                u1[i,j] = u[i+1,j]+u[i,j] 
            end
        end
        return u1 / 2.0
    end

    function B(u)
        u1 = zeros(typeof(u[1]),N)
        u1[1,:] *= 2.0
        for j in 2:N[2]
            for i in 1:N[1]
                u1[i,j] = u[i,j-1]+u[i,j] 
            end
        end
        return u1 / 2.0
    end

    function T(u)
        u1 = zeros(typeof(u[1]),N)
        u1[:,end] *= 2.0
        for j in 1:N[2]-1
            for i in 1:N[1]
                u1[i,j] = u[i,j+1]+u[i,j] 
            end
        end
        return u1 / 2.0
    end

    function system(x)
        T = x[1:N[1]*N[2]]
        P = x[N[1]*N[2]+1:2*N[1]*N[2]]
        T = reshape(T,(N[1],N[2]))
        P = reshape(P,(N[1],N[2]))

        println(size(T))
        println(size(P))

        F1 = (λ/c)*(∇_left(T)+∇_right(T)+∇_top(T)+∇_bottom(T)) +
        ρ_ref*k*(L(T).*∇_left(P) + R(T).*∇_right(P)+T(T).*∇_top(P) + B(T).*∇_bottom(P))-
        (ρ_ref^2)*k*(B(T)-T(T)) +
        ρ_ref*k*α*(B(T).^2 - T(T).^2) 

        F1[end,:] += (1/ϵ) * T[1] - (1/ϵ) * bc_bott
        F1[1,:] += (1/ϵ) * T[end] - (1/ϵ) * bc_top

        F2 = ( ∇_left(P)+∇_right(P)+∇_top(P)+∇_bottom(P) )+
        α*(B(T)-T(T))

        F2[1,:] += (1/ϵ) * P[1,:] 
        reshape(F1,N[1]*N[2])
        reshape(F2,N[1]*N[2])
        return(vcat(F1,F2))
    end
end

sys = new_assemble_nonlinear_system(gridX,gridY,bc_bott,bc_top)
x_0 = zeros(Float64,n*n*2)
#@code_warntype assemble_nonlinear_system(gridX,gridY,bc_bott,bc_top)
#@code_warntype sys(x_0)
sys(x_0)