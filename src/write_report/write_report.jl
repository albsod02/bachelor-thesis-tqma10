"""
    Script

Main code for the synthetic 1D bathymetry reconstructions.

Parses settings, generates synthetic data, runs interpolation and regularization,
and optionally produces error-vs-alpha plots and LaTeX table rows.
"""
module Script

include("errors.jl")
include("spline_interpolation/BSplineInterpolation.jl")
include("synthetic_data_generation/write_bathymetric_data.jl")
include("write_report/write_report.jl")
include("parse_settings.jl")
include("regularization/regularization.jl")
include("regularization/parameter_tuning_regularization.jl")

using Random
using Printf
using CairoMakie
using TrixiBottomTopography

using .Errors
using .SplineInterpolation
using .SyntheticData
using .WriteReport
using .SettingsParser
using .Regularization
using .ParameterTuning

Random.seed!(1234)

const VAR_EPS = 1e-16


"""
    make_interp_at_alpha(alpha::Int, bfun::Function, tag::String, end_point, n::Int,
                         sparsity::Int, add_noise::Bool, end_condition, smoothing_factor;
                         variance_eps::Float64 = VAR_EPS)

Generate synthetic data at noise level `alpha`, write it to a temporary file,
and return a B-spline interpolation function.
"""
function make_interp_at_alpha(alpha::Int, bfun::Function, tag::String, end_point, n::Int, sparsity::Int, add_noise::Bool, end_condition, smoothing_factor; variance_eps::Float64 = VAR_EPS)
    variance = add_noise ? 10.0^(-alpha) : variance_eps
    x, d = define_test_function(bfun, end_point, n, sparsity, add_noise, variance)
    tmp = joinpath("../data", "synthetic", "1D", "tmp_s$(sparsity)_$(tag)_a$(alpha).txt")
    write_bathymetric_data(tmp, x, d)
    return construct_b_spline(tmp, end_condition, smoothing_factor)
end


"""
    make_reg_at_alpha(alpha::Int, bfun::Function, end_point, m::Int, n::Int,
                      sparsity::Int, add_noise::Bool, tau, tol, maxit,
                      dp_lambda_min, dp_lambda_max, lambda_floor;
                      variance_eps::Float64 = VAR_EPS)

Generate synthetic data at noise level `alpha` and return a regularized
reconstruction. If `add_noise`, `lambda` is chosen according to the discrepancy principle.
"""
function make_reg_at_alpha(alpha::Int, bfun::Function, end_point, m::Int,  n::Int, sparsity::Int, add_noise::Bool, tau, tol, maxit, dp_lambda_min, dp_lambda_max, lambda_floor; variance_eps::Float64 = VAR_EPS)
    variance = add_noise ? 10.0^(-alpha) : variance_eps
    x, d = define_test_function(bfun, end_point, n, sparsity, add_noise, variance)

    lambda = add_noise ?
        regularize_discrepancy(d, m, n, variance, tau, tol, maxit, dp_lambda_min, dp_lambda_max) :
        lambda_floor

    return regularize(d, m, n, lambda, variance, end_point)
end


"""
    main()

Run reconstruction pipeline based on `../config/settings.txt`.
Writes data, figures, reports, and optionally LaTeX tables.
"""
function main()
    settings = parse_settings(joinpath("../config", "settings.txt"))

    filename = settings["filename"]
    path_data = joinpath("../data", "synthetic", "1D", filename * ".txt")

    end_point = settings["end_point"]
    n = settings["n"]
    sparsity = settings["sparsity"]
    m = div(n, sparsity)

    add_noise = settings["add_noise"]
    alpha = settings["alpha"]
    variance = add_noise ? 10.0^(-alpha) : VAR_EPS

    end_condition = settings["end_condition"]
    smoothing_factor = settings["smoothing_factor"]
    add_knots = settings["add_knots"]

    tune_lambda_min = settings["lambda_min"]
    tune_lambda_max = settings["lambda_max"]
    num_points = settings["num_points"]

    generate_latex_tables = settings["generate_latex_tables"]
    generate_error_alpha_plots = settings["generate_error_alpha_plots"]

    alpha_max = settings["alpha_max"]
    n_dense = 5000

    # --- Test functions ---
    b1(x) = (1 - 1 / (1 + exp(-30 * (x - 0.5)))) +
            0.05 * (1 - 1 / (1 + exp(-40 * (x - 0.55)))) * sin(2 * pi * 15 * x)

    b2(x) = x < 0.4 ? 0.6x + 0.05 * sin(10 * pi * x) :
                     1 - 0.4 * (1 - x) + 0.03 * cos(8 * pi * x)

    # Choose bathymetric function to run single reconstruction
    b = b1

    # --- Generate data for the single run ---
    x, d = define_test_function(b, end_point, n, sparsity, add_noise, variance)
    write_bathymetric_data(path_data, x, d)

    legend_position = (b === b2) ? (0.0, 1.0) : (1.0, 1.0)

    # --- Interpolation (single run) ---
    interpolation_function = plot_interpolation_with_original_function(
        path_data, b, d, x, n, sparsity, end_point,
        smoothing_factor, end_condition, add_noise, variance, add_knots;
        legend_position = legend_position
    )

    outpath = joinpath(@__DIR__, "..", "reports", "report_interpolation")
    mkpath(dirname(outpath))
    write_report(outpath, b, interpolation_function, end_point, n; n_dense=n_dense)

    # --- Regularization tuning plots (single run) ---
    parameter_tuning_lambda("1",   tune_lambda_min, tune_lambda_max, num_points, b, d, m, n, end_point, variance)
    parameter_tuning_lambda("2",   tune_lambda_min, tune_lambda_max, num_points, b, d, m, n, end_point, variance)
    parameter_tuning_lambda("inf", tune_lambda_min, tune_lambda_max, num_points, b, d, m, n, end_point, variance)

    # --- Regularization (single run) ---
    tau   = 1.1
    tol   = 1e-8
    maxit = 1e6

    dp_lambda_min = 1e-12
    dp_lambda_max = 1e9

    lambda_floor  = tune_lambda_min

    lambda = add_noise ? regularize_discrepancy(d, m, n, variance, tau, tol, maxit, dp_lambda_min, dp_lambda_max) : lambda_floor
    lambda = settings["lambda"]

    b_regularized = regularize(d, m, n, lambda, variance, end_point)
    plot_regularized_solution(x, end_point, b, b_regularized, d, n, m, lambda, alpha; legend_position = legend_position)

    outpath = joinpath(@__DIR__, "..", "reports", "report_regularized")
    mkpath(dirname(outpath))
    write_report(outpath, b, b_regularized, end_point, n; n_dense=n_dense)

    figdir  = joinpath(@__DIR__, "..", "figures", "error_plots", "s$(sparsity)")
    outdir  = joinpath(@__DIR__, "..", "latex", "tables")
    m_label = string(m)

    if generate_error_alpha_plots
        mkpath(figdir)
    end
    if generate_latex_tables
        mkpath(outdir)
    end

    for (tag_local, b_local) in (("smooth", b1), ("disc", b2))

        make_interp = a -> make_interp_at_alpha(a, b_local, tag_local, end_point, n, sparsity, add_noise, end_condition, smoothing_factor; variance_eps=VAR_EPS)
        make_reg    = a -> make_reg_at_alpha(a, b_local, end_point, m, n, sparsity, add_noise, tau, tol, maxit, dp_lambda_min, dp_lambda_max, lambda_floor; variance_eps=VAR_EPS)

        make_interp_inf = () -> make_interp_at_alpha(1, b_local, tag_local, end_point, n, sparsity, false, end_condition, smoothing_factor; variance_eps=VAR_EPS)
        make_reg_inf    = () -> make_reg_at_alpha(1, b_local, end_point, m, n, sparsity, false, tau, tol, maxit, dp_lambda_min, dp_lambda_max, lambda_floor; variance_eps=VAR_EPS)

        if generate_error_alpha_plots
            Errors.plot_error_vs_alpha(joinpath(figdir, "$(tag_local)_int_error_vs_alpha_s$(sparsity).png"), alpha_max, make_interp; b=b_local, end_point=end_point, which=:b,  n_dense=n_dense, title="Error vs alpha")
            Errors.plot_error_vs_alpha(joinpath(figdir, "$(tag_local)_reg_error_vs_alpha_s$(sparsity).png"), alpha_max, make_reg;   b=b_local, end_point=end_point, which=:b,  n_dense=n_dense, title="Error vs alpha")
        end

        if generate_latex_tables
            Errors.write_error_table_rows(joinpath(outdir, "$(tag_local)_int_m$(m_label).tex"), alpha_max, make_interp; b=b_local, end_point=end_point, which=:b,  n_dense=n_dense, include_infty=true, make_reconstruction_infty=make_interp_inf)
            Errors.write_error_table_rows(joinpath(outdir, "$(tag_local)_der_int_m$(m_label).tex"), alpha_max, make_interp; b=b_local, end_point=end_point, which=:db, n_dense=n_dense, include_infty=true, make_reconstruction_infty=make_interp_inf)

            lambda_getter = a -> begin
                if !add_noise
                    return "-"
                end
                var_a = 10.0^(-a)
                _, d_a = define_test_function(b_local, end_point, n, sparsity, true, var_a)
                lam = regularize_discrepancy(d_a, m, n, var_a, tau, tol, maxit, dp_lambda_min, dp_lambda_max)
                @sprintf("%.1f", lam)
            end

            Errors.write_error_table_rows(joinpath(outdir, "$(tag_local)_reg_m$(m_label).tex"), alpha_max, make_reg; b=b_local, end_point=end_point, which=:b,  n_dense=n_dense, include_infty=true, make_reconstruction_infty=make_reg_inf, add_lambda_col=true, lambda_getter=lambda_getter, lambda_cell_infty="-")
            Errors.write_error_table_rows(joinpath(outdir, "$(tag_local)_reg_der_m$(m_label).tex"), alpha_max, make_reg; b=b_local, end_point=end_point, which=:db, n_dense=n_dense, include_infty=true, make_reconstruction_infty=make_reg_inf, add_lambda_col=true, lambda_getter=lambda_getter, lambda_cell_infty="-")
        end
    end

    println("Done")
    return nothing
end

end # module Script

if abspath(PROGRAM_FILE) == @__FILE__
    Script.main()
end
