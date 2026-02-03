module Regularization

export regularize, plot_regularized_solution, regularize_discrepancy, choose_lambda_discrepancy

using CairoMakie
using SparseArrays
using LinearAlgebra

"""
    create_finite_difference_matrix(n) -> L

Create the first-order finite difference matrix `L` of size `(n-1) x n`.
"""
function create_finite_difference_matrix(n)
    L = spzeros(n - 1, n)
    @inbounds for i in 1:(n - 1)
        L[i, i] = -1.0
        L[i, i + 1] = 1.0
    end
    return L
end

"""
    create_matrix_F(m, n) -> F

Create a sampling matrix `F` of size `m x n` that selects `m` sparsely uniform nodes from `n`.
"""
function create_matrix_F(m, n)
    F = spzeros(m, n)
    indices = round.(Int, LinRange(1, n, m))
    @inbounds for (i, idx) in enumerate(indices)
        F[i, idx] = 1.0
    end
    return F
end

"""
    tikhonov_regularization(F, d, L, lambda, variance) -> b_reg

Solve the Tikhonov-regularized least squares problem:

`(F'F + variance * lambda * L'L) b = F'd`.
"""
function tikhonov_regularization(F, d, L, lambda, variance)
    A = F' * F + variance * lambda * (L' * L)
    b_reg = A \ (F' * d)
    return b_reg
end

"""
    vector_to_function(b_reg::Vector{Float64}, end_point::Real) -> b_reg_func

Helper function that converts a coefficient vector `b_reg` into a piecewise-constant function on `[0, end_point]`.
"""
function vector_to_function(b_reg::Vector{Float64}, end_point::Real)
    n = length(b_reg)
    x_grid = range(0, stop=end_point, length=n)
    return x -> begin
        idx = round(Int, (x / end_point) * (n - 1)) + 1
        idx = clamp(idx, 1, n)
        b_reg[idx]
    end
end

"""
    bisection(F, d, L, variance, lambda_min, lambda_max, target, tol, maxit) -> lambda

Bisection search (on a log-scale midpoint) to find `lambda` such that the residual norm
`||F b_lambda - d||₂` matches `target` within `tol`.
"""
function bisection(F, d, L, variance, lambda_min, lambda_max, target, tol, maxit)
    iters = 0
    lambda_mid = lambda_max

    while (lambda_max - lambda_min > tol) && (iters < maxit)
        iters += 1

        lambda_mid = sqrt(lambda_min * lambda_max)

        bmid = tikhonov_regularization(F, d, L, lambda_mid, variance)
        rmid = norm(F * bmid - d)

        if rmid > target
            lambda_max = lambda_mid
        else
            lambda_min = lambda_mid
        end
    end

    return lambda_mid
end

"""
    choose_lambda_discrepancy(F, d, L, variance, tau, tol, maxit, lambda_min, lambda_max) -> lambda_DP

Choose `lambda` using the discrepancy principle:

`||F b_lambda - d||₂ = tau * sigma * sqrt(m)`, where `sigma = sqrt(variance)` and `m = size(F, 1)`.
"""
function choose_lambda_discrepancy(F, d, L, variance, tau, tol, maxit, lambda_min, lambda_max)
    m = size(F, 1)
    target = tau * sqrt(variance) * sqrt(m)
    lambda_DP = bisection(F, d, L, variance, lambda_min, lambda_max, target, tol, maxit)

    return lambda_DP
end

"""
    regularize(d, m, n, lambda, variance, end_point) -> b_reg_func

Compute a Tikhonov-regularized reconstruction with user-specified `lambda` and return it
as a function handle on `[0, end_point]`.
"""
function regularize(d, m, n, lambda, variance, end_point)
    L = create_finite_difference_matrix(n)
    F = create_matrix_F(m, n)
    b_reg = tikhonov_regularization(F, d, L, lambda, variance)
    return vector_to_function(b_reg, end_point)
end

"""
    regularize_discrepancy(d, m, n, variance, tau, tol, maxit, lambda_min, lambda_max) -> lambda_DP

Choose `lambda` by the discrepancy principle and return the selected value `lambda_DP`.
"""
function regularize_discrepancy(d, m, n, variance, tau, tol, maxit, lambda_min, lambda_max)
    L = create_finite_difference_matrix(n)
    F = create_matrix_F(m, n)
    lambda_DP = choose_lambda_discrepancy(F, d, L, variance, tau, tol, maxit, lambda_min, lambda_max)
    b_reg = tikhonov_regularization(F, d, L, lambda_DP, variance)
    return lambda_DP
end

"""
    plot_regularized_solution(x_obs, end_point, b, b_reg, d, n, m, lambda, alpha; outdir::String = "../figures/regularization_plots", tag::String = "", legend_position = :rt)

Plot the true function `b` and the regularized reconstruction `b_reg` on `[0, end_point]`,
together with the observed data `d`, and save the figure.
"""
function plot_regularized_solution(x_obs, end_point, b, b_reg, d, n, m, lambda, alpha;
    outdir::String = "../figures/regularization_plots", tag::String = "", legend_position = :rt)

    if !isdir(outdir)
        mkpath(outdir)
    end

    fig = Figure()
    ax = Axis(fig[1, 1],
        title = "Comparison of True and Reconstructed Bathymetry Function",
        xlabel = "x",
        ylabel = "b"
    )

    x_grid = range(0, stop=end_point, length=n)
    sample_x = range(0, stop=end_point, length=m)

    b_true_vals = b.(x_grid)
    b_reg_vals  = b_reg.(x_grid)

    lines!(ax, x_grid, b_true_vals, color=:"#25539c", linewidth=2, label="Bathymetry b(x)")
    lines!(ax, x_grid, b_reg_vals,  color=:red, linestyle=:dash, linewidth=2, label="Reconstruction")
    scatter!(ax, sample_x, d, color=:black, strokewidth=0.01, markersize=6, marker=:x, label="Observed values")

    axislegend(ax; position=legend_position, padding=5, labelsize=9)

    sparsity = div(n, m)
    name = tag == "" ? "regularized_solution" : "regularized_solution_" * tag
    outfile = joinpath(outdir, "$(name)_sparsity=$(sparsity)_alpha=$(alpha)_lambda=$(floor(lambda, digits=1)).png")
    save(outfile, fig)
end

end  # module
