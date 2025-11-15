#!/usr/bin/env julia

using LinearAlgebra
using CairoMakie

n = 80      
end_point = 1.0   
h = 1 / n # distance between consecutive grid points       

alpha = 3
variance = 1 / 10 ^ alpha
add_knots = true  


# Width and amplitude of Gaussian kernel
omega = 0.04                
A = 1 / (sqrt(2π) * omega)  


b(x) = 0.3 + 0.7 / (1 + exp(-25 * (x - 0.30))) - 0.15 * exp(-((x - 0.70) / 0.07)^2)
s = collect(range(h, end_point, length = n))
b_true = b.(s)

# Gaussian kernel
k(s) = A * exp(-(s^2) / (2 * omega^2))

# F matrix
F = Array{Float64}(undef, n, n)
@inbounds for i in 1:n, j in 1:n
    F[i, j] = h * k(h * (i - j))
end

# Add noise
e = sqrt(variance) .* randn(n)
d = F * b_true .+ e


fig = Figure()
ax = Axis(fig[1, 1], title  = "Example 1.2.2", xlabel = "x", ylabel = "b")

lines!(ax, s, b_true, color="#25539c", linewidth=2, label="Bathymetry b(x)")

if add_knots
    scatter!(ax, s, b_true, color="#9fc0f5", strokecolor=:black, strokewidth=0.4, markersize=4, marker=:circle, label="Nodal values")
end

scatter!(ax, s, d, color=:black, strokewidth=0.01, markersize=6, marker=:x, label="Observed values")

axislegend(ax; orientation=:vertical, position=(0.0, 1.0), padding=5, rowgap=2, colgap=2, labelsize=9)

filename = "figures/example1.2.2.png"
save(filename, fig)

println("Done")