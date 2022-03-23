include("simple_finite_volume.jl")


# main part

# number of cells
N_x = 70
N_y = 35
λ = 0.01
# x,y coordinates of collocation points
coll_x = LinRange(0, 300, N_x)
coll_y = LinRange(0, 150, N_y)

# lists of x,y coordinates of every collocation point
#println(size(transpose((1:N_x)' .* ones(Int64 ,N_y))))
a= vec(transpose((1:N_x)' .* ones(Int64 ,N_y)))
b = vec(transpose(ones(Int64,N_x)' .* (1:N_y)))
cv_list = [(a[i],b[i]) for i in 1:N_x*N_y] # list of control volumes (contains indices of collocation points coordinates)
mesh_1 = Dict("N_x" => N_x, "N_y" => N_y, "coll_x" => coll_x,"coll_y" => coll_y, "cv_list" =>  cv_list) # put it all in the dictionary
T_bc_bott = 50*ones(N_x)  # temperature BCs bottom
T_bc_top  = 0*ones(N_x) # temperature BCs top

A,rhs = get_T_matrix(mesh_1,T_bc_bott,T_bc_top)
T_steady = get_temperature(A,rhs,T_bc_bott,T_bc_top)
P_steady = get_pressure(mesh_1,T_steady)
println("start computing")

T_0 = vcat(0.0*ones(N_x*(N_y-2)÷2),zeros(N_x*(N_y-2)÷2))

t_end1 = 100000.0
t_end2 = 10000.0
t_end3 = 1000.0

println("calculate temperature")
@time T_sol1 = time_temperature(T_0,A,rhs,(0.0,t_end1),T_bc_bott,T_bc_top,1,0.1)
@time T_sol2 = time_temperature(T_0,A,rhs,(0.0,t_end2),T_bc_bott,T_bc_top,0.1,1)
@time T_sol3 = time_temperature(T_0,A,rhs,(0.0,t_end3),T_bc_bott,T_bc_top,10,0.1)
println("calculating pressure")
#P_sol = time_pressure(mesh_1,T_sol)
println("done")
timestamps1 = LinRange(0.0,t_end1,100)
timestamps2 = LinRange(0.0,t_end2,100)
timestamps3 = LinRange(0.0,t_end3,100)

T_errors_1 = [get_L_2_error(mesh_1,T_sol1(t),T_steady) for t in timestamps1]
println("done 1")
T_errors_2 = [get_L_2_error(mesh_1,T_sol2(t),T_steady) for t in timestamps2]
println("done 2")
T_errors_3 = [get_L_2_error(mesh_1,T_sol3(t),T_steady) for t in timestamps3]
println("done 3")

fig = Figure()
ax  = CairoMakie.Axis(fig[1,1],ylabel = L"log (L_2 error)", xlabel = L"time [s]")
lines!(ax,timestamps1,log.(10,T_errors_1),label = L"\frac{\lambda}{c} = 1")
lines!(ax,timestamps2,log.(10,T_errors_2),label = L"\frac{\lambda}{c} = 10")
lines!(ax,timestamps3,log.(10,T_errors_3),label = L"\frac{\lambda}{c} = 100")
axislegend()
save("log_error.png", fig)
display(fig)

