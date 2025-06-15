# module for generating bathymetric data from a smooth test function f(x)

module GenerateData

export define_test_function, plot_function_with_samples

using CairoMakie

# Function to define exact or noisy test function data
function define_test_function(f::Function; a=0, b=1, n=200, add_noise=false, variance=0.01)
    x_vals = range(a, b; length=n)
    y_vals = f.(x_vals)

    if add_noise
        sigma = sqrt(variance)
        y_vals .+= sigma * randn(length(y_vals))
    end

    return x_vals, y_vals
end

"""
# Function to plot smooth curve and sampled data
function plot_function_with_samples(f::Function; a=0, b=1, n_samples=20, variance=0.01, filename="sampled_plot.png")
    x_curve, y_curve = define_test_function(f; a=a, b=b, n=300)
    x_data, y_data = define_test_function(f; a=a, b=b, n=n_samples, add_noise=false, variance=variance)

    fig = Figure()
    ax = Axis(fig[1, 1]; xlabel="x", ylabel="f(x)", title="Test Function with Sample Points")
    lines!(ax, x_curve, y_curve, label="f(x)", color=:blue)
    scatter!(ax, x_data, y_data, label="Noisy samples", color=:red, markersize=6)
    axislegend(ax, position=:rt)

    save(filename, fig)
end
"""

end  #end module
