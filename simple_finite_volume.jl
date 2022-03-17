""" Simple Finite Volume script. Solves ∇ ⋅ λ∇T = 0  on rectangular domain 
    with Neumann BCs on the sides and Dirichlet on top/bottom.
    control volumes(cells) are rectangular(dont have to be all same size)
    boundaries of each cell are always in the middle of 2 collocation points (CP).
    If a cell shares a boundary with the domain, CP is placed 
    on that boundary, so it is easy to apply Dirichlet BCs directly to those CPs
    
"""

using LinearAlgebra
#using SparseArrays
using Plots


function plot_3D(u_h,N_x,N_y) # plot a function of 2 variables on given domain
    x = LinRange(0, 300, N_x)
    x = vec(transpose(x' .* ones(N_y)))
    y = LinRange(0, 150, N_y)
    y = vec(transpose(ones(N_x)' .* y))
    # change view angles in line below: camera = (⋅,⋅)
    plot(x, y, u_h,st = :surface,camera = (60,70))

end
# function that might be useful for plotting(not used yet)
function merge_series!(sp1::Plots.Subplot, sp2::Plots.Subplot)
    append!(sp1.series_list, sp2.series_list)
    Plots.expand_extrema!(sp1[:xaxis], xlims(sp2))
    Plots.expand_extrema!(sp1[:yaxis], ylims(sp2))
    Plots.expand_extrema!(sp1[:zaxis], zlims(sp2))
    return sp1
end
# function that might be useful for plotting(not used yet)
function merge_series!(plt, plts...)
    for (i, sp) in enumerate(plt.subplots)
        for other_plt in plts
            if i in eachindex(other_plt.subplots)
                merge_series!(sp, other_plt[i])
            end
        end
    end
    return plt
end

#plots a 2d heatmap of the solution
function plot_2D(u_h,N_x,N_y)
    coll_x = LinRange(0, 300, N_x)
    coll_y = LinRange(0, 150, N_y)

    a= vec(transpose(coll_x' .* ones(Int64 ,N_y)))
    b = vec(transpose(ones(Int64,N_x)' .* coll_y))
    p1 = plot(a, b, seriestype = :scatter, title = "My Scatter Plot")
    p2 = heatmap(coll_x, coll_y, transpose(reshape(u_h,(N_x,N_y))))
    
    #p3 = merge_series!(p2,p1)
    display(p2)
    
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

# main part

# number of cells
N_x = 80 
N_y = 40

# x,y coordinates of collocation points
coll_x = LinRange(0, 300, N_x)
coll_y = LinRange(0, 150, N_y)

# lists of x,y coordinates of every collocation point
a= vec(transpose((1:N_x)' .* ones(Int64 ,N_y)))
b = vec(transpose(ones(Int64,N_x)' .* (1:N_y)))
cv_list = [(a[i],b[i]) for i in 1:N_x*N_y] # list of control volumes (contains indices of collocation points coordinates)
mesh_1 = Dict("N_x" => N_x, "N_y" => N_y, "coll_x" => coll_x,"coll_y" => coll_y, "cv_list" =>  cv_list) # put it all in the dictionary 

# assemble matrix
A = zeros(N_x*N_y,N_x*N_y)
for i in 1:N_x*N_y # loop over control volumes
    nb, h, σ  = get_cell_info(i,mesh_1)
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

#enforce dirichlet BCs using penalty trick (see his slides on fv 2019-2020)
#bc_values = vcat(0*ones(N_x),10*ones(N_x)) #simple bcs
bc_values = vcat(10*cos.((1:N_x)/10),10*sin.((1:N_x)/10))# more interesting bcs
# first N_x values are bottom BCs next N_x values are top bcs 
bc_list = vcat(collect(1:N_x),collect(N_x*N_y-N_x+1:N_x*N_y)) # indices of collocation points, which are BCs
ϵ = 10^-8 # epsilon for penalty trick
rhs = zeros(N_x*N_y) # right hand side
for k in 1:2*N_x # asseble rhs
    i = bc_list[k]
    A[i,i] += 1/ϵ
    rhs[i] = (1/ϵ) * bc_values[k]
end

u_h = inv(A) * rhs # direct solver is faster here because of the epsilon in the matrix

#plot_3D(u_h,N_x,N_y)
plot_2D(u_h,N_x,N_y)