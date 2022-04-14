using CairoMakie
using SparseArrays
using ForwardDiff
using LinearAlgebra
using DifferentialEquations
using GLMakie
using Sundials
using NLsolve
using BenchmarkTools

n_x = 20
n_y = 50
n_fine = n_y÷3
n_coarse = n_y - n_fine
grid1y = LinRange(0.0,2.0,n_fine)
grid2y = LinRange(2.0,150.0,n_coarse+1)[2:end]
gridy = reverse(vcat(grid1y,grid2y)) #Uncomment this for nonuniform grid but also comment uniform one
gridx = LinRange(0,300.0,n_x)
gridy = reverse(LinRange(0,150.0,n_y)) #uniform grid

gridX = gridx'.*ones(n_y)
gridY = ones(n_x)'.*gridy
N = (n_y,n_x)

h_x = gridX[:,2:N[2]] - gridX[:,1:(N[2]-1)]     
h_y = gridY[1:(N[1]-1),:] - gridY[2:N[1],:]  
h_top = vcat(h_y[1,:]',h_y) # array of distances between collocation points of size N shifted to the top
h_bottom = vcat(h_y,h_y[end,:]')  # array of distances between collocation points of size N shifted to the bottom
h_left = hcat(h_x[:,1],h_x) # array of distances between collocation points of size N shifted to the left
h_right = hcat(h_x,h_x[:,end])  # array of distances between collocation points of size N shifted to the right