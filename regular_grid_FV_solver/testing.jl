using BenchmarkTools

function crate_function(x::AbstractFloat,y::AbstractFloat)
    r = x^2+y^2
    function f(a::Int)
        sum = 0.0
        for i in 1:a
            sum += r
        end
        return sum
    end
    return f
end

function original_function(a::Int,x::AbstractFloat,y::AbstractFloat)
    r = x^2+y^2
    sum = 0.0
    for i in 1:a
        sum += r
    end
    return sum
end




created_function = crate_function(12.7,2342.9)



@time original_function(1000000,12.7,2342.9)

@time created_function(1000000)

@time created_function(1000000)
@time created_function(1000000)



