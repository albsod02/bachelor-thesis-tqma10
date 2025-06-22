module ConstructBSpline

export construct_b_spline, plot_interpolation_with_original_function, interpolation_error

include("../data_conversion/read_data.jl")
using CairoMakie
using TrixiBottomTopography
using .ConvertData

"""
normalize_data(values)

Normalizes a vector of values to the [0, 1] range.

# Arguments
- `values::AbstractVector{<:Real}`: Data values to normalize.

# Returns
- `normalized_values::Vector{Float64}`: Values scaled to the [0, 1] interval.
"""
function normalize_data(values)
    min = minimum(values)
    max = maximum(values)
    return (values .- min) ./ (max - min)
end

"""
construct_b_spline(data, num_interpolation_points, end_condition, smoothing_factor)

Constructs a cubic B-spline interpolant from input data.

# Arguments
- `data::Vector{Float64}`: The bathymetric data to interpolate.
- `num_interpolation_points::Int`: Number of points for evaluating the interpolation.
- `end_condition::String`: Boundary condition ("free" or "not-a-knot").
- `smoothing_factor::Real`: Spline smoothing parameter.

# Returns
- `x_interpolation_points::Vector{Float64}`: X-values at which the spline is evaluated.
- `b_interpolation_points::Vector{Float64}`: Interpolated values.
- `spline_func::Function`: Callable spline interpolation function.
"""
function construct_b_spline(data, num_interpolation_points, end_condition, smoothing_factor)
    spline_struct = CubicBSpline(data; end_condition, smoothing_factor)
    spline_func(x) = spline_interpolation(spline_struct, x)

    x_interpolation_points = Vector(LinRange(spline_struct.x[1], spline_struct.x[end], num_interpolation_points))
    b_interpolation_points = spline_func.(x_interpolation_points)

    return x_interpolation_points, b_interpolation_points, spline_func
end

"""
plot_interpolation_with_original_function(path, f, x, b, n, end_point, smoothing_factor, end_condition, add_noise, variance, add_knots)

Plots a comparison between the original function and the B-spline interpolation.

# Arguments
- `path::String`: Path to input data file.
- `f::Function`: Original test function.
- `x::Vector{Float64}`: Sampled x-values.
- `b::Vector{Float64}`: Corresponding function values.
- `n::Int`: Number of interpolation points.
- `end_point::Real`: End of the interval [0, end_point].
- `smoothing_factor::Real`: Smoothing parameter for the spline.
- `end_condition::String`: Boundary condition for spline interpolation.
- `add_noise::Bool`: For additional noise. 
- `variance::Real`: Noise variance.
- `add_knots::Bool`: Whether to display original sample points.

# Returns
- `spline_func::Function`: The constructed spline interpolation function for later use.
"""
function plot_interpolation_with_original_function(path, f, x, b, n, end_point, smoothing_factor, end_condition, add_noise, variance, add_knots)
    x_interpolated, b_interpolated, spline_func = construct_b_spline(path, n, end_condition, smoothing_factor)

    fig = Figure()
    ax = Axis(fig[1, 1], title="f(x) vs B-spline", xlabel="x", ylabel="b(x)")

    x_dense = range(0, stop=end_point, length=10 * n)
    b_true = f.(x_dense)

    lines!(ax, normalize_data(x_dense), normalize_data(b_true), color=:blue, label="f(x)")
    lines!(ax, normalize_data(x_interpolated), normalize_data(b_interpolated), color=:red, linestyle=:dash, label="B-spline interpolation")

    if add_knots
        scatter!(ax, normalize_data(x), normalize_data(b), color=:black, label="Sample points")
    end

    axislegend(ax; orientation=:horizontal, position=(0.5, 1.02), halign=:center)

    save("spline_interpolation.png", fig)

    return spline_func
end

"""
interpolation_error(norm, f, interpolation_function, end_point, n)

Computes the relative interpolation error between a true function and its spline interpolant
using the specified norm.

# Arguments
- `norm::String`: Norm type. One of `"1"`, `"2"`, or `"inf"`.
- `f::Function`: The original (true) function.
- `interpolation_function::Function`: The interpolated spline function.
- `end_point::Real`: The end of the interval `[0, end_point]`.
- `n::Int`: Number of original sample points. Error is computed using up to `min(100n, 10000)` points.

# Returns
- `rel_error::Float64`: The relative error between the original and interpolated values under the chosen norm.

# Throws
- `ArgumentError` if an unsupported norm string is passed.
"""
function interpolation_error(norm::String, f::Function, interpolation_function::Function, end_point::Real, n::Int)
    n = min(100 * n, 10000)
    x_dense = range(0, stop=end_point, length=n)
    b_true = f.(x_dense)
    b_interpolation = interpolation_function.(x_dense)
    diff = abs.(b_true .- b_interpolation)

    if norm == "1"
        rel_error = sum(diff) / sum(abs.(b_true))
    elseif norm == "2"
        rel_error = sqrt(sum(diff.^2)) / sqrt(sum(b_true.^2))
    elseif norm == "inf"
        rel_error = maximum(diff) / maximum(abs.(b_true))
    else
        throw(ArgumentError("Unsupported norm: \"$norm\". Choose \"1\", \"2\", or \"inf\"."))
    end

    return rel_error
end

end #module