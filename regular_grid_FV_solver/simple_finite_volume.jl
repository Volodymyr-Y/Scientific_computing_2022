""" Simple Finite Volume script. Solves ∇ ⋅ λ∇T = 0  on rectangular domain 
    with Neumann BCs on the sides and Dirichlet on top/bottom.
    Also does time dependent solution with method of lines
    control volumes(cells) are rectangular(dont have to be all same size)
    boundaries of each cell are always in the middle of 2 collocation points (CP).
    If a cell shares a boundary with the domain, CP is placed 
    on that boundary, so it is easy to apply Dirichlet BCs directly to those CPs 
"""

using CSV
using Tables
using LinearAlgebra
using SparseArrays
using IterativeSolvers
using DifferentialEquations
using Images
using FileIO
using GLMakie
using CairoMakie

GLMakie.activate!()

function write_steady_sol(mesh,u_h,filename)
    N_x = mesh["N_x"]
    N_y = mesh["N_y"]
    coll_x = LinRange(0, 300, N_x)
    coll_y = LinRange(0, 150, N_y)
    x = vec(transpose(coll_x' .* ones(Int64 ,N_y)))
    y = vec(transpose(ones(Int64,N_x)' .* coll_y))
    x_y_u = hcat(x,y,u_h)
    CSV.write(filename, Tables.table(x_y_u))
end

function In_cond_from_png(filename)
    img = load(filename)
    img = Gray.(img)
    matrix = reverse(transpose(convert(Array{Float64}, img)),dims=2)
    return 1 .- vec(matrix)
end

function create_gif_3D(mesh,solution,filename,t_end)
    N_x = mesh["N_x"]
    N_y = mesh["N_y"]
    coll_x = LinRange(0, 300, N_x)
    coll_y = LinRange(0, 150, N_y)
    x = vec(transpose(coll_x' .* ones(Int64 ,N_y)))
    y = vec(transpose(ones(Int64,N_x)' .* coll_y))
    fig = Figure()
    time = Observable(0.0)
    u_h = @lift(solution($time))
    ax = Axis3(fig[1,1],title = @lift("t = $(round(($time)/60, digits = 1))"*" min"))
    surface!(ax,x, y, u_h )
    N_frames = 7*30
    fps = 30
    t_sample = LinRange(0,t_end,N_frames)
    record(fig,filename, t_sample;framerate = fps) do t
    time[] = t
    end
end

function plot_3D(u_h,N_x,N_y) # plot a function of 2 variables on given domain
    fig = Figure()
    x = LinRange(0, 300, N_x)
    x = vec(transpose(x' .* ones(N_y)))
    y = LinRange(0, 150, N_y)
    y = vec(transpose(ones(N_x)' .* y))
    # change view angles in line below: camera = (⋅,⋅)
    ax = Axis3(fig[1,1]; aspect=(1, 1, 1))
    hm = surface!(ax,x, y, u_h )
    Colorbar(fig[1, 2], colormap = :viridis,
    flipaxis = false)
    display(fig)
end

#plots a 2d heatmap of the solution
function plot_2D(u_h,N_x,N_y)
    coll_x = LinRange(0, 300, N_x)
    coll_y = LinRange(0, 150, N_y)
    fig = Figure()
    ax = CairoMakie.Axis(fig[1,1])
    a= vec(transpose(coll_x' .* ones(Int64 ,N_y)))
    b = vec(transpose(ones(Int64,N_x)' .* coll_y))
    u= (reshape(u_h,(N_x,N_y)))
    hm = heatmap!(ax, a, b, u)
    Colorbar(fig[1, 2], colormap = :viridis,
    flipaxis = false)
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
    nx = mesh["N_x"]
    ny = mesh["N_y"]
    cx = mesh["coll_x"]
    cy = mesh["coll_y"]
    control_volumes = mesh["cv_list"]

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

 
function get_T_matrix(mesh,bc_bott,bc_top)
    N_x = mesh["N_x"]
    N_y = mesh["N_y"]
    A = zeros(N_x*N_y,N_x*N_y)
    for i in 1:N_x*N_y # loop over control volumes
        nb, h, σ  = get_cell_info(i,mesh)
        #loop over neighbours of the cell
        for k in 1:4
            if !(nb[k] == 0) # if neighbouring cell exists
                j = nb[k] # index of the neigbouring cell
                σ_k = σ[k] # length of the border with neighbouring cell
                h_k = h[k] # distance between the collocation points of given cell and neighbouring cell 
                A[i,j] +=  σ_k/h_k # contribution of flux from neighbouring cells
                A[i,i] += -σ_k/h_k #
            end
        end
    end
    println("full size: ",size(A))
    A = A[N_x+1:end-N_x,N_x+1:end-N_x]
    println("reduced dimension:",size(A))
    #enforce dirichlet BCs using penalty trick (see his slides on fv 2019-2020)
    #bc_values = vcat(0*ones(N_x),10*ones(N_x)) #simple bcs
    bc_values = vcat(bc_bott,bc_top) # bcs
    # first N_x values are bottom BCs next N_x values are top bcs 
    bc_list = vcat(collect(1:N_x),collect(N_x*N_y-N_x+1:N_x*N_y)) # indices of collocation points, which are BCs
    ϵ = 10^-8 # epsilon for penalty trick
    rhs = zeros(N_x*(N_y-2)) # right hand side
    #println("rhs size", size(rhs))
    rhs[1:N_x] -= [get_cell_info(i,mesh)[3][4]/get_cell_info(i,mesh)[2][4] for i in 1:N_x] .* bc_bott
    rhs[end-N_x+1:end] -= [get_cell_info(i,mesh)[3][4]/get_cell_info(i,mesh)[2][3] for i in N_x*N_y-N_x+1:N_x*N_y] .* bc_top
    #println(rhs)
    return A, rhs
end

function get_temperature(L_h,rhs,bc_bott,bc_top)
    #A = sparse(L_h)
    #T_h = cg(A, rhs)
    T_h = inv(L_h)*rhs
    T_h = vcat(bc_bott,T_h,bc_top)
    return T_h
end

function get_pressure(mesh,T_h)
    α = 0.01
    cy = mesh["coll_y"]
    N_x = mesh["N_x"]
    N_y = mesh["N_y"]
    T = reshape(T_h,(N_x,N_y))

    #println(T[:,4])
    dT_dy_1 = (T[:,2] - T[:,1])/(cy[2]-cy[1])
    dT_dy_2 = (T[:,end] - T[:,end-1])/(cy[end]-cy[end-1])
    for i in 2:(N_y-1)
        h_1 = cy[i+1]-cy[i]
        h_2 = cy[i]-cy[i-1]

        dT_dy_i = (h_2*T[:,i+1] - (h_2-h_1)T[:,i] - h_1*T[:,i-1])/(2*h_1*h_2)
        dT_dy_1 = hcat(dT_dy_1,dT_dy_i)
    end
    dT_dy = hcat(dT_dy_1,dT_dy_2)
    dT_dy = vec(dT_dy)
    #print(dT_dy)
    # now assemble matrix
    A = zeros(N_x*N_y,N_x*N_y)
    for i in 1:N_x*N_y # loop over control volumes
        nb, h, σ  = get_cell_info(i,mesh)
        #loop over neighbours of the cell
        for k in 1:4
            if !(nb[k] == 0) # if neighbouring cell exists
                j = nb[k] # index of the neigbouring cell
                σ_k = σ[k] # length of the border with neighbouring cell
                h_k = h[k] # distance between the collocation points of given cell and neighbouring cell 
                A[i,j] += σ_k/h_k # contribution of flux from neighbouring cells
                A[i,i] += -σ_k/h_k #
            end
        end
    end
    bc_values = 0*ones(N_x)
    bc_list = collect(N_x*N_y-N_x+1:N_x*N_y)
    ϵ = 10^-8 # epsilon for penalty trick
    rhs = zeros(N_x*N_y) # right hand side
    for i in 1:N_x*N_y
        σ = get_cell_info(i,mesh)[3]
        #println(σ)
        rhs[i] = dT_dy[i] * α * σ[1] * σ[3] #add forcing function
    end
    for k in 1:N_x # add BCs to rhs
        i = bc_list[k]
        A[i,i] += 1/ϵ
        rhs[i] = 0#(1/ϵ) * bc_values[k]
    end
    A = sparse(A)
    p_h = cg(A, rhs)
    #p_h = inv(A) * rhs
    return p_h
end

function time_temperature(T_0,L_h,rhs,time_interval,bc_bott,bc_top,λ,delta_t)
    
    p = [sparse(L_h),rhs]
    f(u,p,t) = λ*p[1]*u - λ*p[2]
    #f = ODEFunction(f, jac_prototype=Tridiagonal)
    prob = ODEProblem(f,T_0,time_interval,p)
    sol = solve(prob,Euler(),dt=delta_t) #TRBDF2()
    
    full_solution(t) = vcat(bc_bott,sol(t),bc_top)
    return full_solution
end

function time_pressure(mesh,time_temperature)
    P_sol(t) = get_pressure(mesh, time_temperature(t))
    return P_sol
end

function get_flow(mesh,T_h,P_h)
    α = 0.01
    k = 100
    cy = mesh["coll_y"]
    cx = mesh["coll_x"]
    T = reshape(T_h,(N_x,N_y))
    P = reshape(P_h,(N_x,N_y))
    #println("asas")
    dP_dy_1 = (P[:,2] - P[:,1])/(cy[2]-cy[1])
    dP_dy_2 = (P[:,end] - P[:,end-1])/(cy[end]-cy[end-1])

    dP_dx_1 = (P[2,:] - P[1,:])/(cx[2]-cx[1])
    dP_dx_2 = (P[end,:] - P[end-1,:])/(cx[end]-cx[end-1])

    for i in 2:(N_y-1)
        h_1 = cy[i+1]-cy[i]
        h_2 = cy[i]-cy[i-1]

        dP_dy_i = (h_2*P[:,i+1] - (h_2-h_1)*P[:,i] - h_1*P[:,i-1])/(2*h_1*h_2)
        dP_dy_1 = hcat(dP_dy_1,dP_dy_i)
    end

    for i in 2:(N_x-1)
        h_1 = cx[i+1]-cx[i]
        h_2 = cx[i]-cx[i-1]

        dP_dx_i = (h_2*P[i+1,:] - (h_2-h_1)*P[i,:] - h_1*P[i-1,:])/(2*h_1*h_2)
        dP_dx_1 = hcat(dP_dx_1,dP_dx_i)
    end
    dP_dy = hcat(dP_dy_1,dP_dy_2)
    dP_dy = vec(dP_dy)
    dP_dx = transpose(hcat(dP_dx_1,dP_dx_2))
    dP_dx = vec(dP_dx)
    #println(size(dP_dy),size(dP_dx))
    q_x = -k * dP_dx
    q_y = -k * dP_dy .-k + α*k*T_h
    return q_x , q_y
end

function time_flow_x(mesh,time_T,time_P)
    f(t) = get_flow(mesh,time_T(t),time_P(t))[1]
    return f
end

function time_flow_y(mesh,time_T,time_P)
    f(t) = get_flow(mesh,time_T(t),time_P(t))[2]
    return f
end

function time_flow(mesh,time_T,time_P)
    f(t) = get_flow(mesh,time_T(t),time_P(t))
    return f
end

function get_L_2_error(mesh,u_1,u_2)
    nx = mesh["N_x"]
    ny = mesh["N_y"]
    cx = mesh["coll_x"]
    cy = mesh["coll_y"]
    
    areas = [get_cell_info(i,mesh)[3][1] * get_cell_info(i,mesh)[3][3] for i in 1:nx*ny]
    return sqrt(sum(areas.*((u_1-u_2).^2.0)))
end




