using CairoMakie
using SparseArrays
using SparsityDetection
using Cassette
using ForwardDiff
using LinearAlgebra
using DifferentialEquations
using GLMakie
using Sundials
using NLsolve
using BenchmarkTools
using SparseDiffTools
using FiniteDiff
fcalls = 0
function f(y,x) # in-place
  global fcalls += 1
  for i in 2:length(x)-1
    y[i] = x[i-1] - 2x[i] + x[i+1]
  end
  y[1] = -2x[1] + x[2]
  y[end] = x[end-1] - 2x[end]
  nothing
end

function g(x) # out-of-place
  global fcalls += 1
  y = zero(x)
  for i in 2:length(x)-1
    y[i] = x[i-1] - 2x[i] + x[i+1]
  end
  y[1] = -2x[1] + x[2]
  y[end] = x[end-1] - 2x[end]
  y
end

function construct_sparse_jacobian_function(f,sparsity_pattern)
    # this function returns sparse jacobian for function f
    colors = matrix_colors(sparsity_pattern)
    function jacobian(x)
        jac = copy(sparsity_pattern)
        forwarddiff_color_jacobian!(jac, x->f(x), rand(n), colorvec = colors)
        return jac
    end
    return jacobian
    
end

n = 10
input = rand(n)
output = similar(input)
sparsity_pattern = jacobian_sparsity(f,output,input)
#jac = Float64.(sparse(sparsity_pattern))
jac = spdiagm(0 => ones(Float64,n),-1 => ones(Float64,n-1),1 => ones(Float64,n-1))

println(typeof(jac))
display(jac)
#display(ForwardDiff.jacobian(g,input))
colors = matrix_colors(jac)



#FiniteDiff.finite_difference_jacobian!(jac, f, rand(n), colorvec=colors)
forwarddiff_color_jacobian!(jac, f, rand(n), colorvec = colors)

display(jac)