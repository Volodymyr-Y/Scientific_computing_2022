using NLsolve
using GLMakie

#include("simple_finite_volume.jl")

mutable struct Mesh
    λ::Float64
    c::Float64
    α::Float64
    k::Float64
    ρ_ref::Float64
    N_x::Int64
    N_y::Int64
    coll_x::Array{Float64}
    coll_y::Array{Float64}
    x::Array{Int64}
    y::Array{Int64}
    cv_list::Array{Tuple}
    function Mesh(nx,ny)
        cx = LinRange(0, 300, nx)
        cy = LinRange(0, 150, ny)
        x_= vec(transpose((1:nx)' .* ones(Int64 ,ny)))
        y_ = vec(transpose(ones(Int64,nx)' .* (1:ny)))
        cv_lst = [(x_[i],y_[i]) for i in 1:nx*ny]
        new(0.01,0.00001,0.01,100.0,1.0,nx,ny,cx,cy,x_,y_,cv_lst)
    end
end

function plot_3D(u_h,N_x,N_y,label,filename) # plot a function of 2 variables on given domain
    fig = Figure()
    x = LinRange(0, 300, N_x)
    x = vec(transpose(x' .* ones(N_y)))
    y = LinRange(0, 150, N_y)
    y = vec(transpose(ones(N_x)' .* y))
    # change view angles in line below: camera = (⋅,⋅)
    ax = Axis3(fig[1,1]; aspect=(1, 1, 1),xlabel = "x",ylabel = "y",zlabel = label)
    hm = surface!(ax,x, y, u_h )
    Colorbar(fig[1, 2],hm )
    save(filename, fig)
    display(fig)
end
function get_cell_info(i,mesh) 
    """returns: 
    -list of neighbouring cells 
    -distances between collocation points of neighbours and the cell itself
    -length of cell horizontal/vertical boundaries 
    (all info to calculate flux into the cell from neighbours and put it 
    in the correct place in the matrix)
    """
    # unpack mesh object (which is dictionary)
    nx = mesh.N_x
    ny = mesh.N_y
    cx = mesh.coll_x
    cy = mesh.coll_y
    control_volumes = mesh.cv_list

    # find indices of left and right neighbouring cells
    # if there is no nb cell (boundary) put zero 
    if control_volumes[i][1] > 1 && control_volumes[i][1] < nx
        x_neighbours = [i-1,i+1]
    elseif control_volumes[i][1] == 1
        x_neighbours = [0,i+1]
    elseif control_volumes[i][1] == nx
        x_neighbours = [i-1,0]
    end

    # find indices of bottom and top neighbouring cells
    if control_volumes[i][2] > 1 && control_volumes[i][2] < ny
        y_neighbours = [i-nx,i+nx]
    elseif control_volumes[i][2] == 1
        y_neighbours = [0,i+nx]
    elseif control_volumes[i][2] == ny
        y_neighbours = [i-nx,0]
    end
    

    #calculate distances to the left, right, bottom, top neighbouring collocation points  
    if x_neighbours[1] == 0
        h_x_left = 0
    else 
        h_x_left = cx[control_volumes[i][1]] - cx[control_volumes[x_neighbours[1]][1]]
    end
    #
    if x_neighbours[2] == 0
        h_x_right = 0
    else 
        h_x_right = cx[control_volumes[x_neighbours[2]][1]] - cx[control_volumes[i][1]]
    end
    #
    if y_neighbours[1] == 0
        h_y_bott = 0
    else 
        h_y_bott = cy[control_volumes[i][2]] - cy[control_volumes[y_neighbours[1]][2]]
    end
    #
    if y_neighbours[2] == 0
        h_y_top = 0
    else 
        h_y_top = cy[control_volumes[y_neighbours[2]][2]] - cy[control_volumes[i][2]] 
    end
    # calculate cell boundary length
    σ_horizontal = (h_x_left + h_x_right)/2  #length of the cell's edge in x direction
    σ_vertical = (h_y_bott + h_y_top)/2  #length of the cell's edge in y direction

    return vcat(x_neighbours,y_neighbours), [h_x_left,h_x_right,h_y_bott,h_y_top], [σ_vertical,σ_vertical,σ_horizontal,σ_horizontal]
end


function func(x,mesh)
    λ = mesh.λ
    c = mesh.c
    α  =mesh.α
    k = mesh.k
    ρ_ref = mesh.ρ_ref
    N_x = mesh.N_x
    N_y = mesh.N_y
    #coll_x = mesh.coll_x
    #coll_y = mesh.coll_y
    #x = mesh.x
    #y = mesh.y
    #cv_list = mesh.cv_list # list of control volumes (contains indices of collocation points coordinates)
    #mesh = Dict("N_x" => N_x, "N_y" => N_y, "coll_x" => coll_x,"coll_y" => coll_y, "cv_list" =>  cv_list) # put it all in the dictionary
    #T_bc_bott = 3*ones(N_x)  # temperature BCs bottom
    #T_bc_top  = 0*ones(N_x) # temperature BCs top
    T_bc_bott = 3.0*ones(N_x)  # temperature BCs bottom
    T_bc_top = 0.0*ones(N_x)  # temperature BCs bottom
    T = x[1:N_x*(N_y-2)]
    P = x[N_x*(N_y-2)+1:end]
    T_full = vcat(T_bc_bott,T,T_bc_top)
    P_full = vcat(P,zeros(N_x))
    #println(size(T_full))
    F_1 = Float64[]
    F_2 = Float64[]
    #free_deg_freedom = (N_x+1):(N_x*N_y-N_x)
    for i in 1:N_x*N_y 
        nb , h , σ = get_cell_info(i,mesh)
        f1, f2 = 0.0,0.0
        g = [0.0,0.0,1.0,-1.0] # proection of g on left/right/bottom/top cell boundary normal vector   
        for j in 1:4
            if nb[j] == 0
                P_b,T_b,∇T,∇P =  (0.0,0.0,0.0,0.0)
            else
                P_b = (P_full[nb[j]] + P_full[i])/2.0
                T_b = (T_full[nb[j]] + T_full[i])/2.0
                ∇P = (P_full[nb[j]] - P_full[i])/h[j]
                ∇T = (T_full[nb[j]] - T_full[i])/h[j]
            end
            f1 += ((λ/c)*∇T + ρ_ref*k*T_b*∇P- (ρ_ref^2)*k*T_b*g[j] + ρ_ref*k*T_b^2*g[j])*σ[j]
            f2 += (∇P + α*T_b*g[j])*σ[j]
            #println(P_bound) 
        end
        append!(F_1,f1)
        append!(F_2,f2)
    end
    #println("size output")
    #println(size(vcat(F_1[N_x+1:end-N_x],F_2[1:end-N_x])))
    return vcat(F_1[N_x+1:end-N_x],F_2[1:end-N_x])
end


function get_solution(nl_function,N_x,N_y,bc_bott,bc_top)
    x_init = -100.0*ones(2*N_x*N_y -3*N_x)   
    function f!(x, fvec)
        fvec = nl_function(x)
    end
    @time x_out = nlsolve(f!, x_init).zero
    T = x_out[1:N_x*(N_y-2)]
    P = x_out[N_x*(N_y-2)+1:end]
    return vcat(bc_bott,T,bc_top),vcat(P,zeros(N_x))
end

 
#println("size x_init")
#println(size(x_init))
#@time sol = create_nonlinear_equations(N_x,N_y,T_bc_bott,T_bc_top)
#println("time the function itself")
#@time x_out = func!(x_init)
#T,P = get_solution(func!,N_x,N_y,T_bc_bott,T_bc_top)

#plot_3D(T,N_x,N_y,"T","regular_grid_FV_solver/nonlinear_T.png")

m1 = Mesh(4,30)
x_0 = 5*ones(2*m1.N_x*m1.N_y - 3*m1.N_x)
println(func(x_0,m1))
function f!(x,F)
    F = func(x,m1)
end

@time x_out = nlsolve(f!, x_0)