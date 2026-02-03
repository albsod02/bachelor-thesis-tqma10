using CairoMakie

include("read_data.jl")
using .ConvertData

"""
    normalize_data(values::AbstractVector{<:Real}) -> Vector{Float64}

Normalize `values` linearly to the range `[0, 1]`.
"""
function normalize_data(values::AbstractVector{<:Real})
    min_val = minimum(values)
    max_val = maximum(values)
    return (values .- min_val) ./ (max_val - min_val)
end

"""
    monic_eval(coeffs, x) -> y

Evaluate a polynomial at `x` in the monomial (power) basis with coefficients `coeffs`,
interpreted as `p(x) = coeffs[1] + coeffs[2]x + ...`.
"""
function monic_eval(coeffs, x)
    result = zero(x)
    for i in reverse(eachindex(coeffs))
        result = result * x + coeffs[i]
    end
    return result
end

"""
    newton_eval(coeffs, x, x_data) -> y

Evaluate a polynomial at `x` in Newton form using divided-difference coefficients `coeffs`
defined on nodes `x_data`.
"""
function newton_eval(coeffs, x, x_data)
    result = coeffs[end]
    for i in length(coeffs)-1:-1:1
        result = result * (x - x_data[i]) + coeffs[i]
    end
    return result
end

"""
    lagrange_eval(x_data, y_data, x) -> y

Evaluate the Lagrange interpolating polynomial through `(x_data, y_data)` at the point `x`.
"""
function lagrange_eval(x_data, y_data, x)
    n = length(x_data)
    result = 0.0
    for i in 1:n
        term = y_data[i]
        for j in 1:n
            if j != i
                term *= (x - x_data[j]) / (x_data[i] - x_data[j])
            end
        end
        result += term
    end
    return result
end

"""
    solve_linear_system(A, f) -> x

Solve the linear system `A * x = f`.
"""
function solve_linear_system(A, f)
    return A \ f
end

"""
    divided_differences(x, y) -> a

Compute divided differences for Newton interpolation.

Returns `a` such that the Newton-form interpolant is determined by coefficients `a`
on nodes `x`.
"""
function divided_differences(x, y)
    n = length(x)
    a = copy(y)
    for j in 2:n
        for i in n:-1:j
            a[i] = (a[i] - a[i-1]) / (x[i] - x[i-j+1])
        end
    end
    return a
end

"""
    create_polynomial_interpolation(x, y, basis) -> y_fit

Create interpolated values at the sample points `x` for the chosen polynomial basis.

Supported bases: `"monic"`, `"newton"`, `"lagrange"`.
"""
function create_polynomial_interpolation(x, y, basis)
    n = length(x)

    if basis == "monic"
        A = [x_val^p for x_val in x, p in 0:n-1]
        coeffs = solve_linear_system(A, y)
        return [monic_eval(coeffs, x_val) for x_val in x]

    elseif basis == "newton"
        coeffs = divided_differences(x, y)
        return [newton_eval(coeffs, x_val, x) for x_val in x]

    elseif basis == "lagrange"
        return [lagrange_eval(x, y, x_val) for x_val in x]

    else
        error("Unsupported choice of basis: $basis")
    end
end

"""
    plot_interpolation(filepath::String, basis::String, save_location::String)

Read 1D data from `filepath`, perform polynomial interpolation using `basis`, and save
a comparison plot to `save_location`.
"""
function plot_interpolation(filepath::String, basis::String, save_location::String)
    x, b = read_data_1D(filepath)

    x_norm = normalize_data(x)
    b_norm = normalize_data(b)

    b_fit = normalize_data(create_polynomial_interpolation(x_norm, b_norm, basis))

    fig = Figure()
    ax = Axis(fig[1, 1]; xlabel="x", ylabel="b(x)", title="Polynomial interpolation of b(x)")
    lines!(ax, x_norm, b_norm; label="Original Data", color=:blue)
    lines!(ax, x_norm, b_fit; label="Polynomial Fit", color=:red, linestyle=:dash)
    axislegend(ax, position=:lt)
    save(save_location, fig)
end

if length(ARGS) < 2
    println("Usage: julia process_data.jl <datafile> <basis>")
    exit(1)
end

filepath = ARGS[1]
basis = lowercase(ARGS[2])

if !(basis in ["monic", "newton", "lagrange"])
    println("Invalid basis. Choose one of: monic, newton, lagrange")
    exit(1)
end

output_file = "b_polynomial_fit_" * basis * ".png"
plot_interpolation(filepath, basis, output_file)
