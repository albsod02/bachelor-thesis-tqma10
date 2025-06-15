# this seems to work

# using the monic bases or nexton basis results in bad interpolations for most datasets
# but using lagrange basis results in good interpolation of the data

using CairoMakie

# Load the data conversion module
include("read_data.jl")
using .ConvertData

# Normalize a vector of real values to the interval [0, 1]
function normalize_data(values::AbstractVector{<:Real})
    min_val = minimum(values)
    max_val = maximum(values)
    return (values .- min_val) ./ (max_val - min_val)
end

# Evaluate monic polynomial given coefficients at x
function monic_eval(coeffs, x)
    result = zero(x)
    for i in reverse(eachindex(coeffs))
        result = result * x + coeffs[i]
    end
    return result
end

# Evaluate Newton form polynomial using divided differences
function newton_eval(coeffs, x, x_data)
    result = coeffs[end]
    for i in length(coeffs)-1:-1:1
        result = result * (x - x_data[i]) + coeffs[i]
    end
    return result
end

# Evaluate Lagrange interpolation at point x
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

# Solve linear system Ax = f
function solve_linear_system(A, f)
    return A \ f
end

# Compute divided differences for Newton interpolation
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

# Generate interpolated values for a given basis
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

# Plot the original and interpolated data
function plot_interpolation(filepath::String, basis::String, save_location::String)
    x, b = read_data_1D(filepath)

    x_norm = normalize_data(x)
    b_norm = normalize_data(b)

    # Ensure the same normalized values are used for fitting
    b_fit = normalize_data(create_polynomial_interpolation(x_norm, b_norm, basis))

    fig = Figure()
    ax = Axis(fig[1, 1]; xlabel="x", ylabel="b(x)", title="Polynomial interpolation of b(x)")
    lines!(ax, x_norm, b_norm; label="Original Data", color=:blue)
    lines!(ax, x_norm, b_fit; label="Polynomial Fit", color=:red, linestyle=:dash)
    axislegend(ax, position=:lt)
    save(save_location, fig)
end

# Main entry point from command line
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
