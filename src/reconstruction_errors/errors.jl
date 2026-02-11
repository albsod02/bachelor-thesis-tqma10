module Errors

export compute_relative_error, compute_relative_derivative_error, write_error_table_rows, plot_error_vs_alpha

using Printf
using CairoMakie
using LaTeXStrings


"""
    norm_error(norm::String, diff::AbstractArray, b_true::AbstractArray) -> Float64

Compute a relative error from pointwise absolute differences
using the norm specified by `norm`
"""
function norm_error(norm::String, diff::AbstractArray, b_true::AbstractArray)
    if norm == "1"
        return sum(diff) / sum(abs.(b_true))
    elseif norm == "2"
        return sqrt(sum(diff.^2)) / sqrt(sum(b_true.^2))
    elseif norm == "inf"
        return maximum(diff) / maximum(abs.(b_true))
    else
        error("norm must be \"1\", \"2\", or \"inf\"")
    end
end

"""
    compute_relative_error(norm::String, b::Function, reconstructed_function::Function,
                           end_point::Real, n::Int) -> Float64

Compute the relative error of `reconstructed_function` against the true function `b`
on a uniform grid of size `n` over `[0, end_point]`.
"""
function compute_relative_error(norm::String, b::Function, reconstructed_function::Function, end_point::Real, n::Int)
    x = range(0, stop=end_point, length=n)
    b_true = b.(x)
    b_rec  = reconstructed_function.(x)
    diff = abs.(b_true .- b_rec)
    return norm_error(norm, diff, b_true)
end

"""
    compute_relative_derivative_error(norm::String, b::Function, reconstructed_function::Function,
                                      end_point::Real, n::Int) -> Float64

Compute the relative error of the discrete derivative (forward differences) of
`reconstructed_function` against `b` on a uniform grid of size `n`.

The point with the largest pointwise relative error is removed before computing the norm.
"""
function compute_relative_derivative_error(norm::String, b::Function, reconstructed_function::Function, end_point::Real, n::Int)
    if n < 3
        error("need n >= 3")
    end

    x = range(0, stop=end_point, length=n)
    h = x[2] - x[1]

    b_der = Float64[]
    r_der = Float64[]

    # forward difference on the evaluation grid
    for j = 2:n
        push!(b_der, (b(x[j]) - b(x[j-1])) / h)
        push!(r_der, (reconstructed_function(x[j]) - reconstructed_function(x[j-1])) / h)
    end

    diff = abs.(b_der .- r_der)

    # remove the worst pointwise relative error (near discontinuity)
    rel = diff ./ max.(abs.(b_der), eps())
    k = argmax(rel)
    deleteat!(b_der, k)
    deleteat!(diff, k)

    return norm_error(norm, diff, b_der)
end


"""
    latex_error(x::Float64) -> String

Format a number as mantissa/exponent in LaTeX, e.g. `"1.2 \\cdot 10^{-3}"`.
Handles `0.0` and non-finite values.
"""
function latex_error(x::Float64)
    if x == 0.0
        return "0.0 \\cdot 10^{0}"
    end
    if !isfinite(x)
        return "NaN \\cdot 10^{0}"
    end
    e = floor(Int, log10(abs(x)))
    m = x / (10.0^e)
    return @sprintf("%.1f \\cdot 10^{%d}", m, e)
end

"""
    latex_eoc(x::Float64) -> String

Format an EOC value with two decimals.
"""
latex_eoc(x::Float64) = @sprintf("%.2f", x)

"""
    _print_row(io, alpha_label, e1, eoc1, e2, eoc2, einf, eocinf;
               add_lambda_col=false, lambda_cell="-")

Helper to write one LaTeX table row for a given `alpha_label`.
Optionally appends a λ-column.
"""
function _print_row(io, alpha_label::String, e1::Float64, eoc1::String, e2::Float64, eoc2::String, einf::Float64, eocinf::String; add_lambda_col::Bool=false, lambda_cell::String="-")
    s1   = latex_error(e1)
    s2   = latex_error(e2)
    sinf = latex_error(einf)

    if add_lambda_col
        println(io, "$(alpha_label) & \$$(s1)\$ & $(eoc1) & \$$(s2)\$ & $(eoc2) & \$$(sinf)\$ & $(eocinf) & $(lambda_cell) \\\\")
    else
        println(io, "$(alpha_label) & \$$(s1)\$ & $(eoc1) & \$$(s2)\$ & $(eoc2) & \$$(sinf)\$ & $(eocinf) \\\\")
    end
end


"""
    write_error_table_rows(outpath::String, alpha_max::Int, which::String,
                           make_reconstruction::Function, b::Function,
                           end_point::Real, n_dense::Int;
                           include_infty::Bool=false,
                           make_reconstruction_infty::Union{Nothing,Function}=nothing,
                           add_lambda_col::Bool=false,
                           lambda_getter::Union{Nothing,Function}=nothing,
                           lambda_cell_infty::String="-") -> String

Write LaTeX table rows for `alpha = 1:alpha_max` to `outpath`.

- `which = "b"`  computes errors for `b`.
- `which = "db"` computes errors for `b'` (via forward differences).

If `include_infty=true`, appends an `\\infty` row using `make_reconstruction_infty()`.
Optionally appends a λ-column via `add_lambda_col=true` and `lambda_getter(alpha)`.
"""
function write_error_table_rows(outpath::String, alpha_max::Int, which::String, make_reconstruction::Function, b::Function, end_point::Real, n_dense::Int; include_infty::Bool=false, make_reconstruction_infty::Union{Nothing,Function}=nothing, add_lambda_col::Bool=false, lambda_getter::Union{Nothing,Function}=nothing, lambda_cell_infty::String="-")

    dir = dirname(outpath)
    if dir != "" && !isdir(dir)
        mkpath(dir)
    end

    e1   = zeros(alpha_max)
    e2   = zeros(alpha_max)
    einf = zeros(alpha_max)

    for alpha = 1:alpha_max
        rf = make_reconstruction(alpha)

        if which == "b"
            e1[alpha]   = compute_relative_error("1",   b, rf, end_point, n_dense)
            e2[alpha]   = compute_relative_error("2",   b, rf, end_point, n_dense)
            einf[alpha] = compute_relative_error("inf", b, rf, end_point, n_dense)
        elseif which == "db"
            e1[alpha]   = compute_relative_derivative_error("1",   b, rf, end_point, n_dense)
            e2[alpha]   = compute_relative_derivative_error("2",   b, rf, end_point, n_dense)
            einf[alpha] = compute_relative_derivative_error("inf", b, rf, end_point, n_dense)
        else
            error("which must be \"b\" or \"db\"")
        end
    end

    eoc1   = fill("-", alpha_max)
    eoc2   = fill("-", alpha_max)
    eocinf = fill("-", alpha_max)

    for i = 2:alpha_max
        eoc1[i]   = "\$" * latex_eoc(log10(e1[i-1] / e1[i])) * "\$"
        eoc2[i]   = "\$" * latex_eoc(log10(e2[i-1] / e2[i])) * "\$"
        eocinf[i] = "\$" * latex_eoc(log10(einf[i-1] / einf[i])) * "\$"
    end

    open(outpath, "w") do io
        for alpha = 1:alpha_max
            lam_cell = "-"
            if add_lambda_col
                lam_cell = (lambda_getter === nothing) ? "-" : string(lambda_getter(alpha))
            end

            _print_row(io, string(alpha),
                       e1[alpha], eoc1[alpha],
                       e2[alpha], eoc2[alpha],
                       einf[alpha], eocinf[alpha];
                       add_lambda_col=add_lambda_col,
                       lambda_cell=lam_cell)
        end

        if include_infty
            if make_reconstruction_infty === nothing
                error("include_infty=true but make_reconstruction_infty was not provided. Provide a function that returns the noiseless reconstruction.")
            end

            rf_inf = make_reconstruction_infty()

            if which == "b"
                e1i   = compute_relative_error("1",   b, rf_inf, end_point, n_dense)
                e2i   = compute_relative_error("2",   b, rf_inf, end_point, n_dense)
                einfi = compute_relative_error("inf", b, rf_inf, end_point, n_dense)
            else
                e1i   = compute_relative_derivative_error("1",   b, rf_inf, end_point, n_dense)
                e2i   = compute_relative_derivative_error("2",   b, rf_inf, end_point, n_dense)
                einfi = compute_relative_derivative_error("inf", b, rf_inf, end_point, n_dense)
            end

            _print_row(io, "\$\\infty\$",
                       e1i, "-",
                       e2i, "-",
                       einfi, "-";
                       add_lambda_col=add_lambda_col,
                       lambda_cell=lambda_cell_infty)
        end
    end

    return outpath
end

"""
    write_error_table_rows(outpath::String, alpha_max::Int, make_reconstruction::Function;
                           b::Function, end_point::Real, which::Symbol=:b, n_dense::Int=5000,
                           include_infty::Bool=false, make_reconstruction_infty=nothing,
                           add_lambda_col::Bool=false, lambda_getter=nothing,
                           lambda_cell_infty::String="-") -> String

Wrapper using `which = :b` or `:db`.
"""
function write_error_table_rows(outpath::String, alpha_max::Int, make_reconstruction::Function; b::Function, end_point::Real, which::Symbol = :b, n_dense::Int = 5000, include_infty::Bool = false, make_reconstruction_infty::Union{Nothing,Function}=nothing, add_lambda_col::Bool=false, lambda_getter::Union{Nothing,Function}=nothing, lambda_cell_infty::String="-")
    which_str = (which == :b)  ? "b" :
                (which == :db) ? "db" :
                error("which must be :b or :db")

    return write_error_table_rows(outpath, alpha_max, which_str, make_reconstruction, b, end_point, n_dense;
                                  include_infty=include_infty,
                                  make_reconstruction_infty=make_reconstruction_infty,
                                  add_lambda_col=add_lambda_col,
                                  lambda_getter=lambda_getter,
                                  lambda_cell_infty=lambda_cell_infty)
end


"""
    plot_error_vs_alpha(outpath::String, alpha_max::Int, make_reconstruction::Function;
                        b::Function, end_point::Real, which::Symbol=:b,
                        n_dense::Int=5000, title::String="Error vs alpha")

Plot relative errors versus `alpha = 1:alpha_max` and save to `outpath`.
Returns `(alpha_values, e1, e2, einf, fig)`.
"""
function plot_error_vs_alpha(outpath::String, alpha_max::Int, make_reconstruction::Function;
                             b::Function,
                             end_point::Real,
                             which::Symbol = :b,
                             n_dense::Int = 5000,
                             title::String = "Error vs alpha")

    dir = dirname(outpath)
    if dir != "" && !isdir(dir)
        mkpath(dir)
    end

    alpha_values = collect(1:alpha_max)
    e1   = zeros(alpha_max)
    e2   = zeros(alpha_max)
    einf = zeros(alpha_max)

    for alpha in alpha_values
        rf = make_reconstruction(alpha)
        if which == :b
            e1[alpha]   = compute_relative_error("1",   b, rf, end_point, n_dense)
            e2[alpha]   = compute_relative_error("2",   b, rf, end_point, n_dense)
            einf[alpha] = compute_relative_error("inf", b, rf, end_point, n_dense)
        elseif which == :db
            e1[alpha]   = compute_relative_derivative_error("1",   b, rf, end_point, n_dense)
            e2[alpha]   = compute_relative_derivative_error("2",   b, rf, end_point, n_dense)
            einf[alpha] = compute_relative_derivative_error("inf", b, rf, end_point, n_dense)
        else
            error("which must be :b or :db")
        end
    end

    e1p   = max.(e1, eps())
    e2p   = max.(e2, eps())
    einfP = max.(einf, eps())

    fig = Figure()
    ax = Axis(fig[1, 1];
        title  = title,
        xlabel = "alpha",
        ylabel = "Error",
        yscale = log10,
        xticks = (alpha_values, string.(alpha_values)),
        xgridvisible = true,
        ygridvisible = true,
        xgridcolor = (:gray, 0.25),
        ygridcolor = (:gray, 0.25),
    )

    lines!(ax, alpha_values, e1p;   color=:red,   linestyle=:dash,  linewidth=2, label=L"||\cdot||_1")
    lines!(ax, alpha_values, e2p;   color=:green, linestyle=:solid, linewidth=2, label=L"||\cdot||_2")
    lines!(ax, alpha_values, einfP; color=:blue,  linestyle=:dot,   linewidth=2, label=L"||\cdot||_\infty")

    scatter!(ax, alpha_values, e1p;   marker=:circle,    markersize=7, color=:black, strokecolor=:black, strokewidth=1.0)
    scatter!(ax, alpha_values, e2p;   marker=:rect,      markersize=7, color=:white, strokecolor=:black, strokewidth=1.0)
    scatter!(ax, alpha_values, einfP; marker=:utriangle, markersize=7, color=:white, strokecolor=:black, strokewidth=1.0)

    axislegend(ax; position=:lb, framevisible=true, labelsize=16, labelfont=:bold)

    save(outpath, fig)
    return alpha_values, e1, e2, einf, fig
end

end # module
