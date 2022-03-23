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

T_bc_bott = 50*ones(N_x)  # temperature BCs bottom
T_bc_top  = 0*ones(N_x) # temperature BCs top

T_0 = vcat(0.0*ones(N_x*(N_y-2)÷2),zeros(N_x*(N_y-2)÷2))

L_h,rhs = get_T_matrix(mesh_1,T_bc_bott,T_bc_top)

sol = time_temperature(T_0,L_h,rhs,(0.0,200.0),T_bc_bott,T_bc_top,1,0.1)


