module WriteReport

export write_report, plot_error_graph

include("../spline_interpolation/BSplineInterpolation.jl")
using .SplineInterpolation

using CairoMakie

function write_report(filename, f, interpolation_function, end_point, n)
    path = joinpath(filename * ".txt")

    L1_norm = interpolation_error("1", f, interpolation_function, end_point, n)
    L2_norm = interpolation_error("2", f, interpolation_function, end_point, n)
    L_inf_norm = interpolation_error("inf", f, interpolation_function, end_point, n)

    open(path, "w") do io

        println(io, "# Report\n")

        # Error with respect to different norms
        println(io, "L1-norm: ", L1_norm)
        println(io, "L2-norm: ", L2_norm)
        println(io, "L_inf-norm: ", L_inf_norm)
    end
end

# Plot errors of all norms
function plot_error_graph(path, f, x, b, n, end_point, smoothing_factor, end_condition, alpha)

    fig = Figure()
    ax = Axis(fig[1, 1], yscale=log10)  # Set y-axis to log10 scale

    # Parameter for variance sigma^2 = 10^-1 , ... , 10^-alpha
    alpha_values = collect(1:alpha)

    errors_1 = []
    errors_2 = []
    errors_inf = []

    #add_noise = true
    #add_knots = true

    for a in 1:alpha

        var = 1 / 10 ^ a
        interpolation_function = plot_interpolation_with_original_function(path, f, x, b, n, end_point, smoothing_factor, end_condition; add_noise = true, variance = var, add_knots = true)

        push!(errors_1, interpolation_error("1", f, interpolation_function, end_point, n))
        push!(errors_2, interpolation_error("2", f, interpolation_function, end_point, n))
        push!(errors_inf, interpolation_error("inf", f, interpolation_function, end_point, n))

        println(var)
    end

    println(errors_1)

    lines!(ax, alpha_values, errors_1, color=:blue, linewidth=2, label="Interpolation Error")
    #lines!(ax, alpha_values, errors_2, color=:blue, linewidth=2, label="Interpolation Error")
    #lines!(ax, alpha_values, errors_inf, color=:blue, linewidth=2, label="Interpolation Error")
    scatter!(ax, alpha_values, errors_1, color=:blue, markersize=8)
    #scatter!(ax, alpha_values, errors_2, color=:blue, markersize=8)
    #scatter!(ax, alpha_values, errors_inf, color=:blue, markersize=8)

    axislegend(ax; position=:rt)
    ax.xlabel = "alpha"
    ax.ylabel = "Interpolation Error"   

    save("error_plot.png", fig)
end


end #module 