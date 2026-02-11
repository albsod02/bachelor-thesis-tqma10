"""
    example1.2.2.jl

Generate synthetic 1D data, plot `b` and observations `d`, and save a figure to
`project/figures/examples/`.
"""

#!/usr/bin/env julia

using LinearAlgebra
using CairoMakie
using Random

Random.seed!(1234)

# --- Grid / noise ---
n = 80
end_point = 1.0
h = 1 / n

alpha = 3
variance = 1 / 10^alpha
add_knots = true

# --- Gaussian kernel ---
omega = 0.04
A = 1 / (sqrt(2π) * omega)

# --- True bathymetry b(x) and sampled values ---
b(x) = 0.3 + 0.7 / (1 + exp(-25 * (x - 0.30))) - 0.15 * exp(-((x - 0.70) / 0.07)^2)
s = collect(range(h, end_point, length = n))
b_true = b.(s)

# --- Forward operator F (Gaussian convolution) ---
k(s) = A * exp(-(s^2) / (2 * omega^2))

F = Array{Float64}(undef, n, n)
@inbounds for i in 1:n, j in 1:n
    F[i, j] = h * k(h * (i - j))
end

# --- Observations d ---
e = sqrt(variance) .* randn(n)
d = F * b_true .+ e

# --- Plot ---
fig = Figure()
ax = Axis(fig[1, 1], title = "Example 1.2.2", xlabel = "x", ylabel = "b")

lines!(ax, s, b_true, color = "#25539c", linewidth = 2, label = "Bathymetry b(x)")

if add_knots
    scatter!(ax, s, b_true, color = "#9fc0f5", strokecolor = :black, strokewidth = 0.4, markersize = 4, marker = :circle, label = "Nodal values")
end

scatter!(ax, s, d, color = :black, strokewidth = 0.01, markersize = 6, marker = :x, label = "Observed values")

axislegend(ax; orientation = :vertical, position = (0.0, 1.0), padding = 5, rowgap = 2, colgap = 2, labelsize = 9)

# --- Save ---
project_root = normpath(joinpath(@__DIR__, "..", ".."))
outdir = joinpath(project_root, "figures", "examples")
mkpath(outdir)

filename = joinpath(outdir, "example1.2.2.png")
save(filename, fig)

println("Done")
