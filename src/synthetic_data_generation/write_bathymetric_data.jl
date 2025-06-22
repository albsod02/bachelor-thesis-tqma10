module SyntheticData

export define_test_function, write_bathymetric_data

"""
define_test_function(f::Function; end_point, n, add_noise, variance)

Generates a synthetic 1D bathymetric dataset by evaluating a given smooth function `f` over a specified interval. Optionally, Gaussian noise can be added to simulate measurement noise.

# Arguments
- `f::Function`: A smooth 1D function to be evaluated.
- `end_point::Real`: The end of the interval `[0, end_point]` over which the function is defined.
- `n::Int`: Number of data points to sample in the interval.
- `add_noise::Bool`: If `true`, Gaussian noise is added to the data.
- `variance::Real`: Variance (σ²) of the Gaussian noise added if `add_noise` is `true`.

# Returns
- `x_values::Vector{Float64}`: The sampled x-values.
- `b_values::Vector{Float64}`: The corresponding (possibly noisy) values of the function.

"""
function define_test_function(f::Function; end_point, n, add_noise, variance)
    x_values = range(0, end_point; length = n)
    b_values = f.(x_values)

    if add_noise
        sigma = sqrt(variance)
        b_values .+= sigma * randn(length(b_values))
    end

    return collect(x_values), collect(b_values)
end 


"""
write_bathymetric_data(filename, x_data, b_data)

Writes 1D bathymetric data to a text file in a format compatible with `TrixiBottomTopography`.

# Arguments
- `filename::String`: Name of the output file (including `.txt` extension).
- `x_data::Vector{Float64}`: The x-coordinate values.
- `b_data::Vector{Float64}`: The corresponding bathymetric (y) values.

# File Format
The output file has the following structure:
1. A comment line indicating the number of x values
2. The number of x values
3. A comment line followed by the list of x values
4. A comment line followed by the list of y (bathymetry) values

# Notes
- Throws an error if `x_data` and `b_data` are not of the same length.
- Output format is compatible with the Trixi.jl plugin `TrixiBottomTopography`.

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

end # end module