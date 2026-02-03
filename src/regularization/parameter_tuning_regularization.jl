module ParameterTuning

export parameter_tuning_lambda

include("../errors.jl")
include("regularization.jl")
include("../spline_interpolation/BSplineInterpolation.jl")

using CairoMakie

using .Regularization
using .SplineInterpolation
using .Errors

"""
    plot_tuning_curve(norm, lambda_values, error_values, sparsity, variance)

Plot the relative error versus `lambda` and save the figure.
"""
function plot_tuning_curve(norm, lambda_values, error_values, sparsity, variance)
    alpha = Int(round(-log10(variance)))

    fig = Figure()
    ax = Axis(fig[1, 1],xlabel = "lambda", ylabel = "Error", title  = "Relative Error for Tikhonov Regularization (alpha = $alpha)")
    lines!(ax, lambda_values, error_values)

    outdir = joinpath(@__DIR__, "..", "..", "figures", "regularization_tuning_plots")
    mkpath(outdir)

    outfile = joinpath(outdir, "error_vs_lambda_$(norm)-norm_alpha=$(alpha)_sparsity=$(sparsity).png")
    save(outfile, fig)
end

"""
    parameter_tuning_lambda(norm, lambda_min, lambda_max, num_points, b, d, m, n, end_point, variance)

Evaluate the relative error for a range of `lambda` values and save an error-versus-`lambda` plot.
"""
function parameter_tuning_lambda(norm, lambda_min, lambda_max, num_points, b, d, m, n, end_point, variance)
    lambda_values = range(lambda_min, stop=lambda_max, length=num_points + 1)
    error_values = Float64[]

    for lambda in lambda_values
        b_regularized = regularize(d, m, n, lambda, variance, end_point)
        error = compute_relative_error(String(norm), b, b_regularized, end_point, n)
        push!(error_values, error)
    end

    sparsity = div(n, m)
    plot_tuning_curve(String(norm), lambda_values, error_values, sparsity, variance)
end

end # module
