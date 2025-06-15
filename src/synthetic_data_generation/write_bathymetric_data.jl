# module for writing bathymetric data to .txt file of specified format compatible with TrixiBottomTopography

module WriteData

export write_data

include("synthetic_test_function.jl")
using .GenerateData

#f(x) = -2x + 0.3sin(10x) + 0.1sin(30x) + 3.0

#plot_function_with_samples(f; n_samples=20, variance=0.001, filename="function_with_noisy_samples.png")

#x_data, b_data = define_test_function(f, add_noise=false, variance=0.01)


"""
write_data(filename, x_data, b_data)

Writes x and bathymetric data to a text file in the specified format:
- Number of x values
- List of x values
- List of corresponding y (bathymetry) values
"""
function write_data(filename::String, x_data::Vector{Float64}, b_data::Vector{Float64})

    if length(x_data) != length(b_data)
        error("x_data and b_data must be of the same length.")
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