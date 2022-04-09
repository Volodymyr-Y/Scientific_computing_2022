include("simple_finite_volume.jl")


# main part

# number of cells
N_x = 70
N_y = 35
λ = 0.01
ρ_ref = 1
k = 100
# x,y coordinates of collocation points
coll_x = LinRange(0, 300, N_x)
coll_y = LinRange(0, 150, N_y)

# lists of x,y coordinates of every collocation point
#println(size(transpose((1:N_x)' .* ones(Int64 ,N_y))))
x= vec(transpose((1:N_x)' .* ones(Int64 ,N_y)))
y = vec(transpose(ones(Int64,N_x)' .* (1:N_y)))
cv_list = [(x[i],y[i]) for i in 1:N_x*N_y] # list of control volumes (contains indices of collocation points coordinates)
mesh_1 = Dict("N_x" => N_x, "N_y" => N_y, "coll_x" => coll_x,"coll_y" => coll_y, "cv_list" =>  cv_list) # put it all in the dictionary
T_bc_bott = 0.5*ones(N_x)  # temperature BCs bottom
T_bc_top  = 0*ones(N_x) # temperature BCs top
t_end = 200.0
#T_0 = In_cond_from_png("C:/Users/Volodymyr/Documents/Masters/Scientific_computing/project/Scientific_computing_2022/regular_grid_FV_solver/Julia.png")
#T_0 = In_cond_from_png("regular_grid_FV_solver/Julia.png")[1+N_x:end-N_x]
T_0 = vcat(0.0*ones(N_x*(N_y-2)÷2),zeros(N_x*(N_y-2)÷2))
A,rhs = get_T_matrix(mesh_1,T_bc_bott,T_bc_top)
@time T_sol = time_temperature(T_0,A,rhs,(0.0,t_end),T_bc_bott,T_bc_top,1,0.01)
@time P_sol = time_pressure(mesh_1,T_sol)
#@time q_sol = time_flow(mesh_1,T_sol,P_sol)
#@time q_sol_y = time_flow_y(mesh_1,T_sol,P_sol)
#println(maximum(q_sol(20)[1]),maximum(q_sol(20)[2]))
#T_h = get_temperature(A,rhs,T_bc_bott,T_bc_top)
#P_h = get_pressure(mesh_1,T_h)
#us ,vs = get_flow(mesh_1,T_h,P_h)

#plot_2D(us,N_x,N_y)
# f = Figure(resolution = (800, 400))
# ax = CairoMakie.Axis(f[1, 1], backgroundcolor = "black")
# heatmap!(ax,x,y,us)
# display(f)

create_gif_3D(mesh_1,P_sol,"regular_grid_FV_solver/Pressure.gif",t_end)

# fig = Figure()
# time = Observable(0.0)
# u_h1 = lift(T_sol,time)
# u_h2 = lift(P_sol,time)
# q_x = lift(t -> q_sol(t)[1],time)
# q_y = lift(t -> q_sol(t)[2],time)
# #strength = lift((q_x,q_y) -> sqrt.(q_x .^ 2 .+ q_y .^ 2), q_x, q_y )
# #strength = sqrt.(q .^ 2 .+ q .^ 2)
# ax1 = CairoMakie.Axis(fig[1,1],)#title = @lift("t = $(round(($time), digits = 1))"*" s")
# ax2 = CairoMakie.Axis(fig[1,2])
# ax3 = CairoMakie.Axis(fig[2,1])
# ax4 = CairoMakie.Axis(fig[2,2])
# #arrows!(ax1,x, y, q_x, q_y, lengthscale = 0.01,arrowcolor = strength, linecolor = strength)
# #heatmap!(ax1,x, y, u_h1 )
# #heatmap!(ax2,x, y, u_h2 )
# #heatmap!(ax3,x, y, q_x )
# heatmap!(ax4,x, y, q_y )
# N_frames = 7*30
# fps = 30
# t_sample = LinRange(0,t_end,N_frames)
# record(fig,"letters.gif", t_sample;framerate = fps) do t
# time[] = t
# end



# println(size(T_0))
