module SplineInterpolation

export construct_b_spline, plot_interpolation_with_original_function

include("../data_conversion/read_data.jl")
using CairoMakie
using TrixiBottomTopography
using .ConvertData

"""
    construct_b_spline(data, end_condition, smoothing_factor) -> spline_func

Construct a cubic B-spline interpolant from `data` and return it as a callable function `spline_func(x)` using the function CubicBSpliine from the TrixiBottomTopography.jl package.
"""
function construct_b_spline(data, end_condition, smoothing_factor)
    spline_struct = CubicBSpline(data; end_condition, smoothing_factor)
    spline_func(x) = spline_interpolation(spline_struct, x)
    return spline_func
end

"""
    plot_interpolation_with_original_function(path, b, d, x, n, sparsity, end_point, smoothing_factor, end_condition, 
    add_noise, variance, add_knots; outdir::String = "../figures/interpolation_plots", tag::String = "",legend_position = :rt) -> spline_func

Plot the true function `b` and the reconstructed spline on `[0, end_point]`, and save the figure.

The returned value is the reconstructed spline function `spline_func(x)` which is plotted in the figure.
"""
function plot_interpolation_with_original_function(path, b, d, x, n, sparsity, end_point, smoothing_factor, end_condition,
     add_noise, variance, add_knots; outdir::String = "../figures/interpolation_plots", tag::String = "", legend_position = :rt)
    spline_func = construct_b_spline(path, end_condition, smoothing_factor)

    if !isdir(outdir)
        mkpath(outdir)
    end

    fig = Figure()
    ax = Axis(fig[1, 1], title = "Comparison of True and Reconstructed Bathymetry Function", xlabel = "x", ylabel = "b")

    x_eval = range(0, stop=end_point, length=n)
    b_true = b.(x_eval)
    b_interp = spline_func.(x_eval)

    lines!(ax, x_eval, b_true,   color=:"#25539c", linewidth=2, label="Bathymetry b(x)")
    lines!(ax, x_eval, b_interp, color=:red, linestyle=:dash, linewidth=2, label="Reconstructed")

    if add_knots
        scatter!(ax, x_eval, b_true, color=:"#9fc0f5", strokecolor=:black, strokewidth=0.4, markersize=4, marker=:circle, label="Nodal values")
    end

    scatter!(ax, x, d, color=:black, strokewidth=0.01, markersize=6, marker=:x, label="Observed values")

    axislegend(ax; orientation=:vertical, position=legend_position, padding=5, labelsize=9)

    alpha_str = (variance == 0.0) ? "inf" : string(Int(round(-log10(variance))))
    name = tag == "" ? "interpolation" : "interpolation_" * tag
    filename = joinpath(outdir, "$(name)_alpha=$(alpha_str)_sparsity=$(sparsity).png")
    save(filename, fig)

    return spline_func
end

end # module
