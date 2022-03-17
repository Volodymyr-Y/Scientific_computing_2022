using LinearAlgebra
using SparseArrays
using Plots
using DifferentialEquations
#using IterativeSolvers

function get_system(N_x,N_y,BC_top,BC_bot)
    λ =    100
    L_x = 300
    L_y = 150
    Δx = L_x/(N_x-1)
    Δy = L_y/(N_y-1)
    
    S_x = diagm(1 => ones(N_x-1), -1 => ones(N_x-1), 0 => -2*ones(N_x))
    S_x[1,2] +=1
    S_x[end,end-1] +=1
    A = kron(I(N_y-2),S_x)/(Δx^2)
    S_y = diagm(1 => ones(N_y-3), -1 => ones(N_y-3), 0 => -2*ones(N_y-2))

    B = kron(S_y,I(N_x))/(Δy^2)
    L_h = -λ*(A+B)
    f_h = zeros(N_x*(N_y-2))
    f_h[1:N_x] = BC_bot
    f_h[N_x*(N_y-2)-N_x+1:end] = BC_top
    f_h = λ*f_h /(Δy^2)
    return L_h, f_h
end

function get_solution_temperature(bc_top, bc_bot,L_h,f_h)
    sol = sparse(L_h) \ f_h
    return vcat(bc_bot,sol,bc_top)
end

function get_time_solution_temperature(N_x,N_y,bc_top,bc_bot,u0,tspan)
    L_h, f_h = get_system(N_x,N_y,bc_top,bc_bot)
    p = [L_h,f_h]
    temperature(u,p,t) = -p[1]*u + p[2]
    prob = ODEProblem(temperature,u0,tspan,p)
    sol = solve(prob, saveat=0.1)[:]
    new_sol = []
    for s in sol
        append!(new_sol,[vcat(bc_bot,s,bc_top)])
    end
    return new_sol
end


function produce_plot(u_h,N_x,N_y) # plot a function of 2 variables on given domain
    x = LinRange(0, 300, N_x)
    x = vec(transpose(x' .* ones(N_y)))
    y = LinRange(0, 150, N_y)
    y = vec(transpose(ones(N_x)' .* y))
    p = plot(x, y, u_h,st = :surface,camera = (60,70))
    display(p)
end




N_x = 80
N_y = 40
#bc_top = 0*ones(N_x)
bc_top = 10*sin.(collect(1:N_x)/2)
#bc_bot = 10*ones(N_x)
bc_bot = 10*sin.(collect(1:N_x)/2)
u0 = vcat(5*ones((N_y-2)*N_x÷2),0*ones((N_y-2)*N_x÷2))
tspan = (0,10)
L_h,f_h = get_system(N_x,N_y,bc_top,bc_bot)
u_h = get_solution_temperature(bc_top, bc_bot,L_h,f_h)
#produce_plot(u_h,N_x,N_y)

#u_t = get_time_solution_temperature(N_x,N_y,bc_top,bc_bot,u0,tspan)
println("sdd")
produce_plot(u_h,N_x,N_y)
"""
@gif for u in u_t[:]

    x = LinRange(0, 300, N_x)
    x = vec(transpose(x' .* ones(N_y)))
    y = LinRange(0, 150, N_y)
    y = vec(transpose(ones(N_x)' .* y))
    #println(size(y))
    plot(x, y, u ,st = :surface,camera = (40,30))
    

end
"""