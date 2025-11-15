#!/usr/bin/env julia

using LinearAlgebra
using SparseArrays
using CairoMakie

n = 80                 
m = 20                  
end_point  = 1.0   

alpha = 3
variance = 1 / 10 ^ alpha
add_knots = true

b(x) = 0.3 + 0.7 / (1 + exp(-25 * (x - 0.30))) - 0.15 * exp(-((x - 0.70) / 0.07)^2)

x = collect(range(0, stop = end_point, length = n))
b_true = b.(x)


obs_idx = round.(Int, range(1, n, length = m))
x_tilde = x[obs_idx]

# F matrix
rows = 1:m
cols = obs_idx
vals = ones(Float64, m)
F = sparse(rows, cols, vals, m, n)

# Add noise
e = sqrt(variance) .* randn(m)
d = F * b_true + e

fig = Figure()
ax = Axis(fig[1, 1], title = "Example 1.2.1", xlabel = "x", ylabel = "b")

lines!(ax, x, b_true, color="#25539c", linewidth=2, label="Bathymetry b(x)")

if add_knots
    scatter!(ax, x, b_true, color="#9fc0f5", strokecolor=:black, strokewidth=0.4, markersize=4, marker=:circle, label="Nodal values")
end

scatter!(ax, x_tilde, d, color=:black, strokewidth=0.01, markersize=6, marker=:x, label="Observed values")

axislegend(ax; orientation=:vertical, position=(0.0, 1.0), padding=5, rowgap=2, colgap=2, labelsize=9)

filename = "figures/example1.2.1.png"
save(filename, fig)

println("Done")
