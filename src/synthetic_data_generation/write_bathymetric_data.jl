module SyntheticData

export define_test_function, write_bathymetric_data

"""
    define_test_function(b::Function, end_point, n, sparsity, add_noise, variance) -> x_values_sparse, d_values

Generate synthetic data from the function `b` on `[0, end_point]`.

A uniform grid with `n` points is created, where `sparsity` determine the spacing of what data points are kept. 
If the setting `add_noise` is parsed as true, Gaussian noise with variance `variance` is added to the sampled values.
"""
function define_test_function(b::Function, end_point, n, sparsity, add_noise, variance)
    x_values = Vector(LinRange(0, end_point, n))
    b_values = b.(x_values)

    x_values_sparse = x_values[1:sparsity:end]
    d_values = b_values[1:sparsity:end]

    if add_noise
        sigma = sqrt(variance)
        noise = sigma * randn(length(d_values))
        d_values .= d_values .+ noise
    end

    return x_values_sparse, d_values
end

"""
    write_bathymetric_data(filename::String, x_data::Vector{Float64}, b_data::Vector{Float64})

Write bathymetric data to `filename` in a simple text format.

The file contains:
1) the number of samples,  
2) all `x` values (one per line),  
3) all `b` values (one per line).

`x_data` and `b_data` must have the same length.
"""
function write_bathymetric_data(filename::String, x_data::Vector{Float64}, b_data::Vector{Float64})
    if length(x_data) != length(b_data)
        error("x_data and b_data must be of the same length (got $(length(x_data)) and $(length(b_data))).")
    end

    open(filename, "w") do io
        println(io, "# Number of x values")
        println(io, length(x_data))

        println(io, "# x values")
        for x in x_data
            println(io, x)
        end

        println(io, "# y values")
        for y in b_data
            println(io, y)
        end
    end
end

end # module
